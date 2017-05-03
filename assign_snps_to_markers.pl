#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX;
use English;
use Data::Dumper;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use List::Util qw/sum/;

my $markers       = '';
my $posteriors    = '';
my $rewrite       = '';
my $loddifference = 3;
my $species       = '';
my $cross         = '';
my $interrun      = '';
my $jsdir         = '';
my $options_okay  = GetOptions(
    'markers=s'       => \$markers,
    'posteriors=s'    => \$posteriors,
    'rewrite'         => \$rewrite,
    'loddifference=i' => \$loddifference,
    'species=s'       => \$species,
    'cross=s'         => \$cross,
    'interrun'        => \$interrun,
    'jsdir'           => \$jsdir,
);

croak "Specify a marker file with -m"         if $markers eq '';
croak "Specify a posteriors file with -p"     if $posteriors eq '';
croak "Specify a species with -s"             if $species eq '';
croak "Specify a cross with -c"               if $cross eq '';
croak "$markers does not exist"               if !-e $markers;
croak "$posteriors does not exist"            if !-e $posteriors;
croack "JoinSingles directory does not exist" if !-e $jsdir;

my $start = time;
print STDERR "Start time: " . localtime($start) . "\n";

my $outputstub = $cross;
$outputstub .= ".intermediates" if $interrun;

# Load sample information and SNP positions from posteriors file
print STDERR "Loading samples and SNPs\n";
my ( $snps, $samples ) = get_cross_info($posteriors);
my $snpnum = @{$snps};
print "$snpnum lines in posteriors file\n";

# Load markers from marker file and find intermediates if requested
print STDERR "Loading markers";
print STDERR " and finding intermediates" if $interrun;
print STDERR "\n";
my ( $markerhash, $maternals, $marker_groups ) = load_markers( $markers, $samples, $outputstub );

run_JoinSingles2( $markerhash, $posteriors, $snpnum, $rewrite, $outputstub );

choose_intermediates( $markerhash, $marker_groups, $outputstub );

write_new_map( $markerhash, $snps, $snpnum, $outputstub );

write_new_markers( $markerhash, $outputstub );

my $end = time;
print STDERR "End time: " . localtime($end) . "\n";

my $runtime = $end - $start;
my $hour    = int( $runtime / 3600 );
my $min     = int( ( $runtime - $hour * 3600 ) / 60 );
my $sec     = $runtime - $hour * 3600 - $min * 60;
printf STDERR "Run time: %02d:%02d:%02d\n", $hour, $min, $sec;

sub run_JoinSingles2 {
    my ( $markerhash, $posteriors, $snpnum, $rewrite, $outputstub ) = @_;

    my @marker_ids = sort { $a <=> $b } keys %{$markerhash};

    open my $markerposteriors, '>', "$outputstub.marker.posteriors.tmp"
      or croak "Can't open marker posteriors file for writing\n";
    map { print $markerposteriors "$markerhash->{$_}{posterior}\n" } @marker_ids;
    close $markerposteriors;

    open my $markeridfile, '>', "$outputstub.marker_ids.tmp" or croak "Can't open marker IDs file for writing\n";
    print $markeridfile "#\n";
    map { print $markeridfile "0\n" } 1 .. $snpnum;
    map { print $markeridfile "$_\n" } @marker_ids;
    close $markeridfile;

    if ( !-e "$outputstub.JoinSingles2.out" or $rewrite ) {
        print STDERR "Running JoinSingles\n";
        system
"zcat $posteriors | cat - $outputstub.marker.posteriors.tmp | java -cp $jsdir JoinSingles2 $outputstub.marker_ids.tmp lod3Mode=3 data=- > $outputstub.JoinSingles2.out 2> $outputstub.JoinSingles2.err";
    }
    else {
        print STDERR "JoinSingles output exists, skipping\n";
    }

    unlink "$outputstub.marker.posteriors.tmp", "$outputstub.marker_ids.tmp";
}

sub choose_intermediates {
    my ( $markerhash, $marker_groups, $outputstub ) = @_;
    print STDERR "Filtering markers\n";
    open my $jsout, '<', "$outputstub.JoinSingles2.out" or croak "Can't open JoinSingles output!\n";

    my %reverse_lookup;
    for my $marker (keys %{$markerhash}) {
        $reverse_lookup{$markerhash->{$marker}{lookup}}{$marker}++;
    }
    
    my %intermarkers;
    for my $group_id ( keys %{$marker_groups} ) {
        for my $perm ( keys %{ $marker_groups->{$group_id} } ) {
            map { $intermarkers{$_}++ } @{ $marker_groups->{$group_id}{$perm} };
        }
    }

    my %markerprobs;
    while ( my $line = <$jsout> ) {
        next if $line =~ /^#/;
        chomp $line;
        if ( $line =~ /\t/ ) {
            my %f = split /\t/, $line;
            my @ids = sort { $f{$b} <=> $f{$a} } keys %f;
            my @markers = map { { id => $_, lod => $f{$_} } } grep { $_ == $ids[0] or defined $intermarkers{$_} } @ids;

            my $maxlod = $markers[0]{lod};
            my @probs  = map { 10**( $_->{lod} - $maxlod ) } @markers;
            my $sum    = sum @probs;
            for my $m ( 0 .. $#markers ) {
                $markerprobs{ $markers[$m]{id} }{prob} += $probs[$m] / $sum;
                $markerprobs{ $markers[$m]{id} }{snps}++;
            }
        }
    }
    close $jsout;

    map { $markerprobs{$_}{ratio} = $markerprobs{$_}{prob} / $markerprobs{$_}{snps} } keys %markerprobs;

    for my $group_id ( sort { $a <=> $b } keys %{$marker_groups} ) {
        my %groupmarkers;
        my %perm;
        for my $perm ( sort keys %{ $marker_groups->{$group_id} } ) {
            my $permprob = 0;
            my $permsnps = 0;
            for my $marker ( @{ $marker_groups->{$group_id}{$perm} } ) {
                $groupmarkers{$marker}++;
                $permprob += $markerprobs{$marker}{prob};
                $permsnps += $markerprobs{$marker}{snps};
            }
            $perm{$perm}{ratio} = $permprob / $permsnps;
            $perm{$perm}{markers} = join ',', @{ $marker_groups->{$group_id}{$perm} };
        }
        my $maxperm = ( sort { $perm{$b}{ratio} <=> $perm{$a}{ratio} } keys %perm )[0];

        map { delete $groupmarkers{$_} } split /,/, $perm{$maxperm}{markers};
        for my $gm (keys %groupmarkers) {
            next if not defined $reverse_lookup{$gm};
            for my $lookup (keys %{$reverse_lookup{$gm}}) {
                delete $markerhash->{$lookup};
            }
        }
    }

}

sub write_new_map {
    my ( $markerhash, $snps, $snpnum, $outputstub ) = @_;
    print STDERR "Writing map\n";
    open my $jsout,        '<', "$outputstub.JoinSingles2.out" or croak "Can't open JoinSingles output!\n";
    open my $mapout,       '>', "$outputstub.map.tsv"          or croak "Can't open map output for writing!\n";
    open my $rejectsout,   '>', "$outputstub.rejects.tsv"      or croak "Can't open rejects output for writing!\n";
    open my $candidateout, '>', "$outputstub.candidates.tsv"   or croak "Can't open candidates file for writing!\n";
    open my $markersout,   '>', "$outputstub.markerpositions.tsv"
      or croak "Can't open marker positions output for writing!\n";

    my $header = "Species\tCross\tScaffold\tPosition\tChromosome\tcM\tMarker\tMarkerType\tLOD\tLODDiff\n";
    print $mapout $header;
    print $rejectsout $header;
    print $candidateout $header;

    my $snpcount = 0;

    while ( my $line = <$jsout> ) {
        next if $line =~ /^#/;
        chomp $line;
        my $mapline       = '';
        my $markerline    = '';
        my $rejectline    = "-1\t-1\t-1\t-1";
        my $markerlod     = 0;
        my $markerloddiff = 0;

        my $snppos = join "\t", @{ $snps->[$snpcount] };

        my @markers;
        if ( $line =~ /\t/ ) {
            my %f = split /\t/, $line;
            @markers =
              map { { id => $_, lod => $f{$_} } } sort { $f{$b} <=> $f{$a} } grep { defined $markerhash->{$_} } keys %f;
        }

        if (@markers) {
            if ( @markers > 1 ) {
                $markerloddiff = $markers[0]{lod} - $markers[1]{lod};
                if ( $markerloddiff < $loddifference ) {
                    unshift @markers, { id => 0, lod => $markerloddiff };
                }
            }
            else {
                $markerloddiff = $markers[0]{lod};
            }

            $markerloddiff = sprintf( "%.3f", $markerloddiff );
            my $choice = $markers[0]{id};
            if (    $choice ne '0'
                and defined $markerhash->{$choice}
                and defined $markerhash->{ $markerhash->{$choice}{lookup} } )
            {
                my $lookup = $markerhash->{$choice}{lookup};
                $markerhash->{$lookup}{onmap}++;
                $markerline =
"$markerhash->{$lookup}{chromosome}\t$markerhash->{$lookup}{cM}\t$lookup\t$markerhash->{$lookup}{type}";
                $mapline = $markerline;
                $markerlod = sprintf( "%.3f", $markers[0]{lod} );
            }
            else {
                $markerline .= "0\t0\t0\t0";
                $rejectline = "0\t0\t0\t0";
            }
            for my $marker (@markers) {
                if (    defined $markerhash->{ $marker->{id} }
                    and defined $markerhash->{ $markerhash->{ $marker->{id} }{lookup} } )
                {
                    my $lookup = $markerhash->{ $marker->{id} }{lookup};
                    $markerline .=
"\t$marker->{id}\t$markerhash->{$lookup}{chromosome}:$markerhash->{$lookup}{cM}\t$lookup\t$marker->{lod}";
                    if ( $choice eq '0' and $marker->{lod} > $markers[1]{lod} - $loddifference ) {
                        my $candidatelod     = sprintf( "%.3f", $marker->{lod} );
                        my $candidateloddiff = sprintf( "%.3f", $markers[1]{lod} - $marker->{lod} );
                        print $candidateout
"$species\t$cross\t$snppos\t$markerhash->{$lookup}{chromosome}\t$markerhash->{$lookup}{cM}\t$lookup\t$markerhash->{$lookup}{type}\t$candidatelod\t$candidateloddiff\n";
                    }
                }
                else {
                    $markerline .= "\t$marker->{id}\t$marker->{lod}";
                }
            }
        }
        print $markersout "$snppos\t$markerline\n";
        print $mapout "$species\t$cross\t$snppos\t$mapline\t$markerlod\t$markerloddiff\n" if $mapline;
        print $rejectsout "$species\t$cross\t$snppos\t$rejectline\t0\t0\n" if not $mapline;
        $snpcount++;
        last if $snpcount == $snpnum;
    }
    close $markersout;
    close $mapout;
    close $rejectsout;
    close $candidateout;
    close $jsout;
}

sub write_new_markers {
    my ( $markers, $outputstub ) = @_;

    open my $newmarkersout, '>', "$outputstub.markers.updated.tsv"
      or croak "Can't open new markers file for writing!\n";
    my @newmarker_ids = sort {
             $markers->{$b}{type} <=> $markers->{$a}{type}
          or $markers->{$a}{chromosome} <=> $markers->{$b}{chromosome}
          or $markers->{$a}{cM} <=> $markers->{$b}{cM}
    } grep { defined $markers->{$_}{onmap} } keys %{$markers};

    for my $nmi (@newmarker_ids) {
        my $m = $markers->{$nmi};
        print $newmarkersout "$nmi\t$m->{chromosome}\t$m->{cM}\t$m->{type}\t$m->{pattern}\n";
    }
    close $newmarkersout;
}

sub load_markers {
    my ( $markerfilename, $samples, $outputstub ) = @_;
    open my $markerfile, '<', $markerfilename or croak "Can't open $markerfilename\n";
    my %markers;
    my %maternals;
    my %marker_groups;
    my %intermediates;
    my $final_marker_id = 0;
    while ( my $line = <$markerfile> ) {
        chomp $line;
        my ( $marker_id, $chromosome, $cM, $type, $pattern ) = split "\t", $line;

        if ( $type == 2 ) {
            $maternals{$marker_id} = $pattern;
        }

        if (    $interrun
            and $marker_id > 1
            and $chromosome == $markers{ $marker_id - 1 }{chromosome}
            and $type == $markers{ $marker_id - 1 }{type} )
        {
            get_intermediates(
                $markers{ $marker_id - 1 }{pattern},
                $pattern, $marker_id, $chromosome, $markers{ $marker_id - 1 }{cM},
                $cM, $type, \%intermediates
            );
        }

        my $posteriors = make_posteriors( $pattern, $type, $samples );

        $markers{$marker_id} = {
            chromosome => $chromosome,
            cM         => $cM,
            type       => $type,
            pattern    => $pattern,
            posterior  => "$marker_id\t-\t$posteriors",
            lookup     => $marker_id
        };

        if ( $type == 1 ) {
            $markers{$marker_id}{intercross} = make_intercross( $pattern, $maternals{$chromosome}, $samples );
        }

        $final_marker_id = $marker_id > $final_marker_id ? $marker_id : $final_marker_id;
    }
    close $markerfile;

    if ($interrun) {
        $final_marker_id =
          add_intermediates( \%markers, \%marker_groups, \%intermediates, $samples, $final_marker_id, \%maternals,
            $outputstub );
    }

    add_intercross( \%markers, $final_marker_id );

    \%markers, \%maternals, \%marker_groups;
}

sub get_intermediates {
    my ( $last, $this, $marker_id, $chromosome, $lastcM, $thiscM, $type, $intermediates ) = @_;
    my %newints;
    croak "Patterns $this and $last are different lengths" if length $last != length $this;
    my @last = split //, $last;
    my @this = split //, $this;
    my @mismatches = get_mismatches( \@last, \@this );
    if ( @mismatches > @this / 2 ) {
        @this = split //, phase($this);
        @mismatches = get_mismatches( \@last, \@this );
    }
    return if @mismatches < 2 or @mismatches > 4;
    my $cMdiff = $thiscM - $lastcM;

    for my $perm ( permute(@mismatches) ) {
        my @new = split //, $last;
        my $permlist = join ',', @{$perm};
        for my $i ( 0 .. $#{$perm} - 1 ) {    # -1 because last entry is the last marker so doesn't need adding
            my $mismatch = $perm->[$i];
            $new[$mismatch] = $new[$mismatch] == 0 ? 1 : 0;
            my $newint = join '', @new;
            if ( not defined $newints{$newint} ) {
                $newints{$newint} = {
                    chromosome => $chromosome,
                    cM         => sprintf( "%.3f", $lastcM + $cMdiff / @mismatches * ( $i + 1 ) ),
                    type       => $type,
                    next_id    => $marker_id,
                };
            }
            push @{ $newints{$newint}{perms} }, $permlist;
        }
    }
    for my $newint ( keys %newints ) {
        if ( defined $intermediates->{$newint} ) {
            print "CLASH! $newint\n";
            print Dumper $intermediates->{$newint};
            print Dumper $newints{$newint};
        }
        else {
            $intermediates->{$newint} = $newints{$newint};
        }
    }
}

sub get_mismatches {
    my ( $last, $this ) = @_;
    return grep { $_ ne '' } map { $_ if $this->[$_] ne $last->[$_] } 0 .. $#{$this};
}

sub add_intermediates {
    my ( $markers, $marker_groups, $intermediates, $samples, $final_marker_id, $maternals, $outputstub ) = @_;
    open my $intout, '>', "$outputstub.tsv" or croak "Can't open intermediate output file for writing\n";
    for my $int (
        sort {
                 $intermediates->{$a}{chromosome} <=> $intermediates->{$b}{chromosome}
              or $intermediates->{$a}{cM} <=> $intermediates->{$b}{cM}
              or $a cmp $b
        } keys %{$intermediates}
      )
    {
        $final_marker_id++;
        my $posteriors = make_posteriors( $int, $intermediates->{$int}{type}, $samples );

        my $intercross = make_intercross( $int, $maternals->{ $intermediates->{$int}{chromosome} }, $samples );

        $markers->{$final_marker_id} = {
            chromosome => $intermediates->{$int}{chromosome},
            cM         => $intermediates->{$int}{cM},
            type       => $intermediates->{$int}{type},
            pattern    => $int,
            posterior  => "$final_marker_id\t-\t$posteriors",
            lookup     => $final_marker_id,
            intercross => $intercross
        };
        my $next_id = $intermediates->{$int}{next_id};
        my $prev_id = $next_id - 1;
        my $nextm   = $markers->{$next_id};
        my $prevm   = $markers->{$prev_id};
        print $intout "$prev_id\t$prevm->{pattern}\t$prevm->{chromosome}\t$prevm->{cM}\n";
        print $intout
          "$final_marker_id\t$int\t$markers->{$final_marker_id}{chromosome}\t$markers->{$final_marker_id}{cM}\n";
        print $intout "$next_id\t$nextm->{pattern}\t$nextm->{chromosome}\t$nextm->{cM}\n\n";

        for my $perm ( @{ $intermediates->{$int}{perms} } ) {
            push @{ $marker_groups->{$next_id}{$perm} }, $final_marker_id;
        }
    }
    close $intout;

    return $final_marker_id;
}

sub add_intercross {
    my ( $markers, $final_marker_id ) = @_;
    my @marker_ids = sort { $a <=> $b } keys %{$markers};
    for my $marker_id (@marker_ids) {
        if ( defined $markers->{$marker_id}{intercross} ) {
            $final_marker_id++;
            $markers->{$final_marker_id} = {
                chromosome => $markers->{$marker_id}{chromosome},
                cM         => $markers->{$marker_id}{cM},
                type       => $markers->{$marker_id}{type},
                pattern    => $markers->{$marker_id}{intercross},
                posterior  => "$final_marker_id\t-\t$markers->{$marker_id}{intercross}",
                lookup     => $marker_id
            };
        }
    }
}

sub make_posteriors {
    my ( $pattern, $type, $samples ) = @_;
    my @gts = split //, $pattern;
    my @posteriors;

    for my $sample ( @{$samples} ) {
        my $posterior;
        if ( $sample->[6] eq 'Father' ) {
            $posterior =
              $type == 1 ? "0 1 0 0 0 0 0 0 0 0" : $type == 2 ? "1 0 0 0 0 0 0 0 0 0" : "1 1 1 1 1 1 1 1 1 1";
        }
        elsif ( $sample->[6] eq 'Mother' ) {
            $posterior =
              $type == 1 ? "1 0 0 0 0 0 0 0 0 0" : $type == 2 ? "0 1 0 0 0 0 0 0 0 0" : "1 1 1 1 1 1 1 1 1 1";
        }
        else {
            my $gt = shift @gts;
            $posterior =
              $gt eq '0' ? "1 0 0 0 0 0 0 0 0 0" : $gt eq '1' ? "0 1 0 0 0 0 0 0 0 0" : "1 1 1 1 1 1 1 1 1 1";
        }
        push @posteriors, $posterior;
    }
    return join "\t", @posteriors;
}

sub make_intercross {
    my ( $paternal, $maternal, $samples ) = @_;
    my @pgts = split //, $paternal;
    my @mgts = split //, $maternal;
    my @intercross =
      map { ( $pgts[$_] == 0 && $mgts[$_] == 0 ) ? 0 : ( $pgts[$_] == 1 && $mgts[$_] == 1 ) ? 2 : 1 } 0 .. $#pgts;

    my @posteriors;
    for my $sample ( @{$samples} ) {
        my $posterior;
        if ( $sample->[6] eq 'Father' ) {
            $posterior = "0 1 0 0 0 0 0 0 0 0";
        }
        elsif ( $sample->[6] eq 'Mother' ) {
            $posterior = "0 1 0 0 0 0 0 0 0 0";
        }
        else {
            my $gt = shift @intercross;
            $posterior =
                $gt eq '0' ? "1 0 0 0 0 0 0 0 0 0"
              : $gt eq '1' ? "0 1 0 0 0 0 0 0 0 0"
              : $gt eq '2' ? "0 0 0 0 1 0 0 0 0 0"
              :              "1 1 1 1 1 1 1 1 1 1";
        }
        push @posteriors, $posterior;
    }
    return join "\t", @posteriors;
}

sub get_cross_info {
    my ($posteriorfilename) = @_;
    my @snps;
    my @samples;
    my $posteriorfile = new IO::Uncompress::Gunzip $posteriorfilename
      or croak "Failed to unzip $posteriorfilename: $GunzipError\n";
    while ( my $line = $posteriorfile->getline() ) {
        next if $line =~ /^#/;
        chomp $line;
        my @f = split /\t/, $line;
        if ( $line =~ /^CHR/ ) {
            shift @f;    # CHR
            shift @f;    # POS
            for my $i ( 0 .. $#f ) {
                push @{ $samples[$i] }, $f[$i];
            }
        }
        else {
            $f[1] =~ s/\*//;
            push @snps, [ $f[0], $f[1] ];
        }
    }

    for my $sample (@samples) {
        if ( $sample->[2] eq '0' and $sample->[3] eq '0' ) {
            if ( $sample->[4] eq '1' ) {
                $sample->[6] = "Father";
            }
            elsif ( $sample->[4] eq '2' ) {
                $sample->[6] = "Mother";
            }
            else {
                croak "Can't identify type of sample @{$sample}\n";
            }
        }
        else {
            $sample->[6] = "Offspring";
        }
    }
    close $posteriorfile;
    \@snps, \@samples;
}

sub phase {
    my ($marker) = @_;
    $marker =~ tr/01/10/;
    return $marker;
}

sub permute {
    return ( [] ) unless (@_);
    return map {
        my @cdr = @_;
        my $car = splice @cdr, $_, 1;
        map { [ $car, @$_ ]; } &permute(@cdr);
    } 0 .. $#_;
}
