#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX;
use English;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/min max/;
use Memoize;
use Storable qw(dclone);

memoize 'get_paternal_from_intercross';

my $crossprefix = "";
my $chromosomes = 0;
my $species     = '';
my $verbose     = '';
my $output      = '';
my $kosambi     = '';

my $options_okay = GetOptions(
    'prefix=s'      => \$crossprefix,
    'chromosomes=s' => \$chromosomes,
    'species=s'     => \$species,
    'output'        => \$output,
    'verbose'       => \$verbose,
    'kosambi'       => \$kosambi,
);

croak
  "Please specify a cross prefix with -p (to match existing files PREFIX_id.err and PREFIX_id.txt)"
  if $crossprefix eq "";
croak "Please specify number of chromosomes with -c" if $chromosomes <= 0;
croak "Please specify a species label with -s" if $species eq '';

croak "$crossprefix\_id.err does not exist" if !-e "$crossprefix\_id.err";
croak "$crossprefix\_id.txt does not exist" if !-e "$crossprefix\_id.txt";

my $start = time;
print STDERR "Start time: " . localtime($start) . "\n";

my @snps       = load_snps($crossprefix);
my @markers    = load_markers($crossprefix);
my @markersnps = load_markersnps($crossprefix);
my %errors     = load_errors($crossprefix);
my %reverses   = load_reverses($crossprefix);

my %paternals;
my %maternals;

for my $chr ( 1 .. $chromosomes ) {
    my %summary = get_markers( \%maternals, $crossprefix, $chr, \@snps, \@markers, \@markersnps, \%errors );

    my %patterns = get_patterns( \%summary, \%paternals, \%maternals, $chr );

    add_new_patterns( \%patterns, \%{ $errors{$chr} } );

    %patterns = reorder_map( \%patterns, defined $reverses{$chr}, $verbose, $kosambi );

    for my $cM ( sort { $a <=> $b } keys %patterns ) {
        $paternals{$chr}{$cM} = $patterns{$cM}{pattern};
    }
}

my $seenall = 1;
for my $chr ( sort { $a <=> $b } keys %errors ) {
    for my $cM ( sort { $a <=> $b } keys %{ $errors{$chr} } ) {
        if ( not defined $errors{$chr}{$cM}{seen} ) {
            $seenall = 0;
            print "Didn't see error $chr $cM $errors{$chr}{$cM}{fix}\n";
        }
    }
}
print "Processed all errors\n" if $seenall;

if ($output) {
    open my $markeroutput, ">", "$crossprefix.markers.tsv"
      or croak "Can't open $crossprefix.markers.tsv for writing\n";

    write_maternals( \%maternals, $markeroutput );

    write_paternals( \%paternals, $markeroutput, 1 + keys %maternals );

    close $markeroutput;
}

my $end = time;
print STDERR "End time: " . localtime($end) . "\n";

my $runtime = $end - $start;
my $hour    = int( $runtime / 3600 );
my $min     = int( ( $runtime - $hour * 3600 ) / 60 );
my $sec     = $runtime - $hour * 3600 - $min * 60;
printf STDERR "Run time: %02d:%02d:%02d\n", $hour, $min, $sec;

sub load_snps {
    my $crossprefix = shift;
    my @snps;
    print STDERR "Loading SNPs file $crossprefix\_snps.txt\n";
    open my $snpfile, "<", "$crossprefix\_snps.txt" or croak "Can't open $crossprefix\_snps.txt\n";
    while ( my $line = <$snpfile> ) {
        next if $line =~ /^#/;
        chomp $line;
        $line =~ s/\*//g;
        push @snps, $line;
    }
    close $snpfile;

    @snps;
}

sub load_markers {
    my $crossprefix = shift;
    my @markers;
    open my $markerfile, "<", "$crossprefix\_id.err" or croak "Can't open $crossprefix\_id.err\n";
    print STDERR "Loading marker file $crossprefix\_id.err\n";
    while ( my $line = <$markerfile> ) {
        next if $line !~ /^\d/;
        chomp $line;
        my ( $id, $type, $count, $paternal, $maternal ) = split "\t", $line;
        $markers[$id] = {
            type     => $type,
            count    => $count,
            paternal => $paternal,
            maternal => $maternal
        };
    }
    close $markerfile;
    @markers;
}

sub load_markersnps {
    my $crossprefix = shift;
    my @markersnps;
    print STDERR "Loading marker SNPs file $crossprefix\_id.txt\n";
    open my $markersnpfile, "<", "$crossprefix\_id.txt" or croak "Can't open $crossprefix\_id.txt\n";
    while ( my $line = <$markersnpfile> ) {
        next if $line =~ /^#/;
        chomp $line;
        push @markersnps, $line;
    }
    close $markersnpfile;

    @markersnps;
}

sub load_errors {
    my $crossprefix = shift;
    my %errors;
    if ( -e "$crossprefix.errors.tsv" ) {
        print STDERR "Loading errors file $crossprefix.errors.tsv\n";
        open my $errorfile, '<', "$crossprefix.errors.tsv" or croak "Can't open $crossprefix.errors.tsv\n";
        while ( my $line = <$errorfile> ) {
            chomp $line;
            my ( $chromosome, $cM, $destination ) = split /\t/, $line;
            $errors{$chromosome}{$cM}{fix} = $destination;
        }
    }
    else {
        print STDERR "No errors file found\n";
    }
    %errors;
}

sub load_reverses {
    my $crossprefix = shift;
    my %reverses;
    if ( -e "$crossprefix.reverse.txt" ) {
        print STDERR "Loading reverses from $crossprefix.reverse.txt\n";
        open my $reversefile, '<', "$crossprefix.reverse.txt" or croak "Can't open $crossprefix.reverse.txt\n";
        while ( my $line = <$reversefile> ) {
            chomp $line;
            $reverses{$line}++;
        }
    }
    else {
        print STDERR "No reverses file found\n";
    }
    %reverses;
}

sub get_markers {
    my ( $maternals, $crossprefix, $chr, $snps, $markers, $markersnps, $errors ) = @_;
    my %summary;
    open my $orderfile, "<", "$crossprefix.order$chr.txt" or croak "can't open $crossprefix.order$chr.txt\n";
    my $phases    = '';
    my $unique_id = -1;
    while ( my $line = <$orderfile> ) {
        chomp $line;
        my @f = split /\t/, $line;

        my $marker_number = $f[0];
        next if $marker_number !~ /^\d+$/;

        my $paternal_cm = $f[1];
        if ( defined $errors->{$chr}{$paternal_cm} and $errors->{$chr}{$paternal_cm}{fix} !~ /^[01]+$/ ) {
            $errors->{$chr}{$paternal_cm}{seen}++;
            next if $errors->{$chr}{$paternal_cm}{fix} eq '-';
            $paternal_cm = $errors->{$chr}{$paternal_cm}{fix};
        }

        my $error = "-";
        if ( $f[3] =~ /\( (.+) \)/ ) {
            $error = $1;
        }

        my $marker = $markers->[ $markersnps->[ $marker_number - 1 ] ];
        if ( $line !~ "duplicate" ) {
            $phases = $f[-1];
            if ( $marker->{type} == 2 ) {
                $maternals->{$chr}{ $marker->{maternal} } = { count => $marker->{count} }
                  if $marker->{type} == 2;
            }
            else {
                $unique_id++;
                $summary{$paternal_cm}{ $marker->{type} }{$unique_id} = {
                    count    => $marker->{count},
                    error    => $error,
                    phase    => $phases,
                    maternal => $marker->{maternal},
                    paternal => $marker->{paternal}
                };
            }
        }

        push @{ $summary{$paternal_cm}{ $marker->{type} }{$unique_id}{pos} }, $snps->[ $marker_number - 1 ]
          if $marker->{type} != 2;
    }
    close $orderfile;

    my $maternal =
      ( sort { $maternals{$chr}{$b}{count} <=> $maternals{$chr}{$a}{count} } keys %{ $maternals{$chr} } )[0];
    $maternals{$chr} = $maternal;

    %summary;
}

sub get_id {
    my $markertype = shift;
    return (
        sort {
                 $markertype->{$a}{error} <=> $markertype->{$b}{error}
              or $markertype->{$b}{count} <=> $markertype->{$a}{count}
        } keys %{$markertype}
    )[0];
}

sub get_patterns {
    my ( $summary, $paternals, $maternals, $chr ) = @_;
    my %patterns;

    for my $cM ( sort { $a <=> $b } keys %{$summary} ) {
        if ( not defined $summary->{$cM}{1} and defined $summary->{$cM}{3} ) {
            my $id = get_id( $summary->{$cM}{3} );
            $summary->{$cM}{1}{$id} = dclone $summary->{$cM}{3}{$id};
            $summary->{$cM}{1}{$id}{maternal} = '?' x length( $summary->{$cM}{1}{$id}{maternal} );
            $summary->{$cM}{1}{$id}{paternal} =
              get_paternal_from_intercross( $summary->{$cM}{1}{$id}{paternal}, $maternals->{$chr} );
        }
        next if not defined $summary->{$cM}{1};
        my $id     = get_id( $summary->{$cM}{1} );
        my $marker = $summary->{$cM}{1}{$id};

        next if $marker->{error} != 0 or $marker->{paternal} =~ /\?/;
        $patterns{$cM}{pattern} =
          ( substr( $marker->{phase}, 0, 1 ) eq '0' ) ? $marker->{paternal} : phase( $marker->{paternal} );
        $patterns{$cM}{count} = $marker->{count};
        $patterns{$cM}{pos}   = $marker->{pos};
    }
    %patterns;
}

sub add_new_patterns {
    my ( $patterns, $errors ) = @_;
    for my $cM ( keys %{$errors} ) {
        if ( $errors->{$cM}{fix} =~ /^[01]+$/ ) {
            $patterns->{$cM}{pattern} = $errors->{$cM}{fix};
            $patterns->{$cM}{count}   = 0;
            $patterns->{$cM}{pos}     = [];
            $errors->{$cM}{seen}++;
        }
    }
}

sub reorder_map {
    my ( $patterns, $reverse, $verbose, $kosambi ) = @_;

    my @cMs = sort { $a <=> $b } keys %{$patterns};

    my ( $diffs, $length ) = get_diffs( \@cMs, $patterns, $verbose );
    my @newcMs = ("0.000");
    for my $diff ( @{$diffs} ) {
        my $numdiff   = $diff =~ tr/X//;
        my $r         = $numdiff / $length;
        if ($kosambi) {
            $r   = log( ( 1 + 2 * $r ) / ( 1 - 2 * $r ) ) / 4;
        }
        my $newcMdist = $r * 100;
        push @newcMs, sprintf "%.3f", $newcMs[-1] + $newcMdist;
    }
    if ($reverse) {
        my $maxcm = max @newcMs;
        @newcMs = map { sprintf "%.3f", $maxcm - $_ } @newcMs;
    }
    my %newpatterns;
    my @order = 0 .. $#newcMs;
    for my $i ( $reverse ? reverse(@order) : @order ) {
        $newpatterns{ $newcMs[$i] } = $patterns->{ $cMs[$i] };
        my $diffout = ( $verbose and defined $diffs->[$i] ) ? "\t\t$diffs->[$i]\n" : '';
        print $diffout if $reverse;
        print "$newcMs[$i]\t$cMs[$i]\t$patterns->{$cMs[$i]}{pattern}\t$patterns->{$cMs[$i]}{count}\t"
          . get_range( $patterns->{ $cMs[$i] }{pos} ) . "\n"
          if $verbose;
        print $diffout if not $reverse;
    }
    %newpatterns;
}

sub get_diffs {
    my ( $cMs, $patterns, $verbose ) = @_;
    my @diffs;
    my @output;
    my $length = 0;
    for my $i ( 1 .. $#{$cMs} - 1 ) {
        $length = length $patterns->{ $cMs->[$i] }{pattern};
        my @last = split //, $patterns->{ $cMs->[ $i - 1 ] }{pattern};
        my @this = split //, $patterns->{ $cMs->[$i] }{pattern};
        my @next = split //, $patterns->{ $cMs->[ $i + 1 ] }{pattern};
        my $last_to_this = @diffs ? $diffs[-1] : join '', map { $last[$_] eq $this[$_] ? ' ' : 'X' } 0 .. $#this;
        if ( not @diffs ) {
            push @diffs, $last_to_this;
        }

        my $this_to_next = join '', map { $this[$_] eq $next[$_] ? ' ' : 'X' } 0 .. $#this;
        push @diffs, $this_to_next;

        my @last_diffs = split //, $diffs[-2];
        my @this_diffs = split //, $diffs[-1];
        for my $j ( 0 .. $#last_diffs ) {
            if ( $last_diffs[$j] ne ' ' and $this_diffs[$j] ne ' ' ) {
                $last_diffs[$j] = '*';
                $this_diffs[$j] = '*';
            }
        }
        $diffs[-2] = join '', @last_diffs;
        $diffs[-1] = join '', @this_diffs;

    }
    \@diffs, $length;
}

sub get_range {
    my ($poslist) = @_;
    my %pos;
    for my $pos ( @{$poslist} ) {
        my ( $scaffold, $position ) = split /\t/, $pos;
        $pos{$scaffold}{$position}++;
    }
    my @pos;
    for my $scaffold ( sort { substr( $a, 5 ) <=> substr( $b, 5 ) } keys %pos ) {
        my @scfpos = sort { $a <=> $b } keys %{ $pos{$scaffold} };
        push @pos, sprintf "%s:%8d-%8d", $scaffold, $scfpos[0], $scfpos[-1];
    }
    join "\t", @pos;
}

sub write_maternals {
    my ( $maternals, $markeroutput ) = @_;
    for my $chr ( sort { $a <=> $b } keys %{$maternals} ) {
        print $markeroutput "$chr\t$chr\t-1\t2\t$maternals->{$chr}\n";
    }
}

sub write_paternals {
    my ( $paternals, $markeroutput, $marker_id ) = @_;

    for my $chr ( sort { $a <=> $b } keys %{$paternals} ) {
        for my $cM ( sort { $a <=> $b } keys %{ $paternals->{$chr} } ) {
            print $markeroutput "$marker_id\t$chr\t$cM\t1\t$paternals->{$chr}{$cM}\n";
            $marker_id++;
        }
    }
}

sub hamming {
    return ( $_[0] ^ $_[1] ) =~ tr/\001-\255//;
}

sub phase {
    my ($marker) = @_;
    $marker = $marker;
    $marker =~ tr/01/10/;
    return $marker;
}

sub get_paternal_from_intercross {
    my $print    = shift;
    my $maternal = shift;

    my @intgts = split //, $print;

    my $mirror = $maternal;
    $mirror =~ tr/01/10/;

    my $matpat = make_paternal( \@intgts, $maternal );
    my $mirpat = make_paternal( \@intgts, $mirror );

    $matpat =~ tr/?// < $mirpat =~ tr/?// ? $matpat : $mirpat;
}

sub make_paternal {
    my ( $intgts, $maternal ) = @_;
    my @matgts = split //, $maternal;
    my $paternal;
    my @patgts;
    for my $i ( 0 .. $#{$intgts} ) {
        my $pgt = ( $intgts->[$i] eq '?' ) ? '?' : $intgts->[$i] - $matgts[$i];
        $pgt = ( $pgt eq '?' or $pgt != 0 and $pgt != 1 ) ? '?' : $pgt;
        push @patgts, $pgt;
    }
    join '', @patgts;
}
