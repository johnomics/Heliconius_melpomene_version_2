#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Memoize;
use DBD::SQLite;
use Parallel::ForkManager;
use Term::ExtendedColor qw/:all/;

use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

memoize "mirror";
memoize "phase";
memoize "linked";
memoize "calc_LR";

my %chromosomes = (
    "ABBBBBBABABAABAAAABBBBBAABABBABABAABAAABABBBABBABBABBAABBABBAABBAABBB" => "1",
    "ABBBABABBAABBBBBBABBAABBBBBBABBBABBABABAABBAAABBABBABBAAAAAAAABBAABAB" => "2",
    "ABABABAABAABBBBBBABBBAABBBBAABAABAABAABAAAAABAABBAAAAAAAABBBBABBABBBA" => "3",
    "ABBABBBBBABBABAAAAABAABBBBBBAAAAABABBAABABABAABABBAABABBAAABBBAABBABB" => "4",
    "AABABBABBBABABABBBBABBBBBAABBABABAABABBABAAAAABBBAAABBAABBBBABBBBBABA" => "5",
    "AABABBBBBABBAABABAABAABABABAAAAABAAAABAABBBBAABAABBBBAAABABAABAAAABAB" => "6",
    "AABBAABBABBBAAABBBBBBBBABAAABBAAABAAAAAAABBAAABAAAABAABBBBABABABAABBB" => "7",
    "AABBABBAABBBAAABAABABAAAAABBABBAAAABBABAABABABBBAABAABAAABABAAAAAAAAB" => "8",
    "AABBBABBABAABBABBABBBAABABBBABBABAAAAAABABABAAAAABBBBBAABBABBAABABBAA" => "9",
    "BAHHBABHBBBHABBAABBABHBABAABABABAHBABAAAABAABAAAAAAAABHBBHBBHAABHBAAA" => "10",
    "ABAABBBABBAABAABBAAAAAABABBAABAAAAABBAAAAAAABABBABAAABBABAAAABABABBAB" => "11",
    "AAABBAABBBBABBBABABABBABBAAAAABAAABBABAABBBABAABAAABAAABBBBBBAABBBBAA" => "12",
    "ABAABBBBAAAABABAABBBAAAAAABBBBBBBBAAAAABAAABABABBBABABAABAAABBAABBBAA" => "13",
    "AAABAHHBABHBAHHABBHHBBHAAAABBABBAHBBAAHBAAHBAHBBBHBAABHABBHBBBAABHAHB" => "14",
    "ABBBBBABBAAAABBBBBBABBABBBAAABAABBABAAAAAAAABBABAABABBABBAAAAABAABBAA" => "15",
    "ABBAABBABBBABBBBBBAAAABABBAAABBBBAABABAABBAAABAAABBBBBBBBBABBBAAABABB" => "16",
    "ABBBAABAABBABAABABABABBAABAABBABBAABBAAAABAAABBBAAAABAABBBAAABBBAABAB" => "17",
    "ABAAAABAAAABAABAAAAABAABBAABABBAABAABBAABBBAAAABBAABABBBABBABBBAAABBB" => "18",
    "ABBAABBAABBBAABBBBBAAABBBBAABBBAABAABBBABABBAAABAABBAABBAAAAAABBAABAA" => "19",
    "ABBBABAAABBBBABAABBBBAABBABBABBABBABBBBBBBBBAAABBAAAABBBBBBBABBBBAAAA" => "20",
    "ABBABABBBBAABBBBBBBAAAAABABAAABBBBABABABBBABBABABBAAAAAABABBAABBAABAB" => "21",
);

my %args;
$args{input}       = "";
$args{output}      = "test";
$args{known}       = "";
$args{corrections} = "";
$args{verbose}     = "";

my $options_okay = GetOptions(
    'input=s'       => \$args{input},
    'output=s'      => \$args{output},
    'known=s'       => \$args{known},
    'corrections=s' => \$args{corrections},
    'verbose'       => \$args{verbose}
);

croak "Please supply an input database with -i" if $args{input} eq "";

croak "Known markers file $args{known} does not exist! $OS_ERROR\n" unless $args{known} eq "" or -e $args{known};

croak "Corrections file $args{corrections} does not exist! $OS_ERROR\n"
  unless $args{corrections} eq ""
  or -e $args{corrections};

my $metadata = { verbose => $args{verbose} };

print STDERR "INFER MARKERS\n";
my $blocklist = refine_markers( $args{input}, $args{corrections}, $metadata );

print STDERR "MAKING LINKAGE MAPS\n";
my ( $linkage_groups, $marker_blocks, $maternal_blocks ) = make_linkage_maps( $blocklist, $args{output}, $args{known} );

output_marker_blocks( $linkage_groups, $marker_blocks, $maternal_blocks, $blocklist, $args{input}, $args{output} );

print STDERR "Done\n";

exit;

## INFER MARKERS
sub refine_markers {
    my ( $input, $corrections, $metadata ) = @_;

    print STDERR "Loading blocks...\n";
    my $blocklist = load_blocks( $input, $corrections, $metadata );

    correct_paternal( $blocklist, $metadata );

    collapse( $blocklist, $metadata );

    output_blocks( $blocklist, $metadata );

    $blocklist;
}

sub load_blocks {
    my ( $input, $corrections, $metadata ) = @_;

    my $blocklist = [];
    my $header;
    my $types;

    my %correct;
    if ( $corrections ne '' ) {
        open my $correctfile, '<', $corrections or croak "Can't open corrections file $corrections! $OS_ERROR\n";
        while ( my $correctline = <$correctfile> ) {
            chomp $correctline;
            my ( $wrong, $right, $comment ) = split /\t/, $correctline;
            $correct{$wrong} = $right;
        }
        close $correctfile;
    }

    print STDERR "Load cleanblocks from database...\n" if $metadata->{verbose};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    ( $header, $types ) = get_header_types($dbh);

    my $sth = $dbh->prepare("SELECT * FROM cleanblocks ORDER BY scaffold, start");
    my $fileblocklist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    push @{$blocklist}, @{$fileblocklist};
    $sth->finish;

    $dbh->disconnect;

    # Make empty pattern by checking first Maternal pattern (could be any pattern)
    my $empty;
    my $samplenum;
    for my $block ( @{$blocklist} ) {
        if ( $block->{Maternal} ne "" ) {
            $samplenum = length $block->{Maternal};
            $empty     = ' ' x $samplenum;
            last;
        }
    }

    print STDERR "Correcting patterns...\n";
    for my $block ( @{$blocklist} ) {
        for my $field ( keys %{$block} ) {
            if ( defined $correct{ $block->{$field} } ) {
                $block->{$field} = $correct{ $block->{$field} };
            }
        }
    }

    $metadata->{input}   = $input;
    $metadata->{header}  = $header;
    $metadata->{samples} = $samplenum;
    $metadata->{types}   = $types;
    $metadata->{empty}   = $empty;

    get_block_stats( "After loading blocks", $blocklist, $metadata );

    $blocklist;
}


sub correct_paternal {
    my ( $blocklist, $metadata ) = @_;

    my $patterns = get_paternal_patterns($blocklist, $metadata);

    my $pattern_num = scalar keys %{$patterns};
    for my $pattern ( sort { $patterns->{$b}{length} <=> $patterns->{$a}{length} } keys %{$patterns} ) {

        for my $upgrade (
            sort { $patterns->{$pattern}{validity}{$b}{length} <=> $patterns->{$pattern}{validity}{$a}{length} }
            keys %{ $patterns->{$pattern}{validity} }
          )
        {
            for my $b ( @{ $patterns->{$pattern}{validity}{$upgrade}{blocks} } ) {
                $blocklist->[$b]{validity} = "   MVPVV";
            }
        }
    }

    update_block_stats( "After correcting paternal", $blocklist, $metadata );

    return;

}

sub get_paternal_patterns {
    my ($blocklist, $metadata) = @_;
    my %patterns;
    for my $b ( 0 .. $#{$blocklist} ) {
        my $block    = $blocklist->[$b];
        if ($block->{validity} =~ /I[ed]/ or $block->{validity} =~ /[X\-]/) {
            $block->{'Maternal'} = $metadata->{empty};
            $block->{'Paternal'} = $metadata->{empty};
        };

        my $paternal = $block->{'Paternal'};
        $patterns{$paternal}{length} += $block->{length};
        $patterns{$paternal}{blocks}++;
        $patterns{$paternal}{validity}{ $block->{validity} }{length} +=
          $block->{length};
        push @{ $patterns{$paternal}{validity}{ $block->{validity} }{blocks} }, $b;
        $patterns{$paternal}{maternal}{ $block->{'Maternal'} }{length} += $block->{length};
        $patterns{$paternal}{maternal}{ $block->{'Maternal'} }{blocks}++;
    }

    # Collapse mirrors
    for my $paternal ( sort { $patterns{$a}{length} <=> $patterns{$b}{length} } keys %patterns ) {
        next if !defined $patterns{$paternal};
        my $mirror = mirror($paternal);
        if ( defined $patterns{$mirror} ) {
            $patterns{$mirror}{length} += $patterns{$paternal}{length};
            $patterns{$mirror}{blocks} += $patterns{$paternal}{blocks};
            for my $validity ( keys %{ $patterns{$paternal}{validity} } ) {
                $patterns{$mirror}{validity}{$validity}{length} += $patterns{$paternal}{validity}{$validity}{length};
                for my $b ( @{ $patterns{$paternal}{validity}{$validity}{blocks} } ) {
                    $blocklist->[$b]{'Paternal'} = $mirror;
                }
                push @{ $patterns{$mirror}{validity}{$validity}{blocks} },
                  @{ $patterns{$paternal}{validity}{$validity}{blocks} };
            }
            for my $maternal ( keys %{ $patterns{$paternal}{maternal} } ) {
                $patterns{$mirror}{maternal}{$maternal}{length} += $patterns{$paternal}{maternal}{$maternal}{length};
                $patterns{$mirror}{maternal}{$maternal}{blocks} += $patterns{$paternal}{maternal}{$maternal}{blocks};
            }
            delete $patterns{$mirror};
        }
    }

    \%patterns;
}


sub collapse {
    my ( $blocklist, $metadata ) = @_;
    for my $block ( @{$blocklist} ) {
        delete $block->{'validity'};
        delete $block->{'IntercrossHet'};
        delete $block->{'IntercrossHom'};
    }

    my $i = 0;
    my $j = 1;
    while ( $j <= $#{$blocklist} ) {
        my $same = 1;
        $same = 0
          if ( $blocklist->[$i]{'scaffold'} ne $blocklist->[$j]{'scaffold'} );

        map {
            if ( $blocklist->[$i]{$_} ne $blocklist->[$j]{$_} ) {
                my $phase_j = phase( $blocklist->[$j]{$_} );
                if ( $blocklist->[$i]{$_} eq $phase_j ) {
                    $blocklist->[$j]{$_} = $phase_j;
                }
                else {
                    $same = 0;
                }
            }
        } ( 'Maternal', 'Paternal' );

        if ($same) {
            $blocklist->[$i]{'end'} = $blocklist->[$j]{'end'};
            $blocklist->[$i]{'length'} =
              $blocklist->[$i]{'length'} + $blocklist->[$j]{'length'};
            splice( @{$blocklist}, $j, 1 );
        }
        else {
            $i = $j;
            $j++;
        }
    }

    return;
}


sub output_blocks {
    my ( $blocklist, $metadata ) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$metadata->{input}", "", "", { AutoCommit => 0 } );

    $dbh->do("DROP TABLE IF EXISTS mapblocks");

    $dbh->do(
"CREATE TABLE mapblocks (scaffold text, start integer, end integer, length integer, Maternal text, Paternal text)"
    );

    $dbh->commit;

    my $insert = $dbh->prepare_cached("INSERT INTO mapblocks VALUES (?,?,?,?,?,?)");

    for my $block ( @{$blocklist} ) {
        $insert->execute(
            $block->{scaffold}, $block->{start},    $block->{end},
            $block->{length},   $block->{Maternal}, $block->{Paternal}
        );
    }

    $dbh->commit;
    $dbh->disconnect;

}



## MAKE LINKAGE MAPS

sub make_linkage_maps {
    my ( $blocklist, $output, $known ) = @_;

    my ( $linkage_groups, $marker_blocks, $maternal_blocks ) = build_linkage_groups( $blocklist, $known, $metadata );

    my $pm = new Parallel::ForkManager( scalar keys %{$linkage_groups} );

    for my $chromosome_print ( sort { $chromosomes{$a} <=> $chromosomes{$b} } keys %{$linkage_groups} ) {
        my $chromosome = $chromosomes{$chromosome_print} // '-';
        $pm->start() and next;
        my $lg         = $linkage_groups->{$chromosome_print};
        my $marker_num = keys %{ $linkage_groups->{$chromosome_print} };
        printf STDERR
          "Building map for print $chromosome_print, chromosome %3s.\t%4d paternal markers to process\n",
          $chromosome, $marker_num;

        my ( $revised_lg, $rejects ) = build_initial_chromosome_map( $lg, $chromosome_print, $output );

        revise_chromosome_map( $lg, $revised_lg, $rejects, $chromosome_print, $output );

        output_rejected_markers( $chromosome_print, $rejects, $lg, $output );

        output_map_summary( $chromosome, $lg, $rejects );

        make_linkage_map( $chromosome_print, $lg, $output );

        $pm->finish();
    }
    $pm->wait_all_children;

    # Collapse rejected markers into one file
    system('cat test.[ABH]*.rejected.tsv | sort -k1,1n -k4,4 -k2,2 > test.rejected.tsv');

    ( $linkage_groups, $marker_blocks, $maternal_blocks );
}


sub build_linkage_groups {
    my ( $blocklist, $known, $metadata ) = @_;

    my %linkage_groups;
    my %marker_blocks;
    my %maternal_blocks;
    
    for my $b ( 0 .. $#{$blocklist} ) {
        my $block = $blocklist->[$b];

        my $maternal = $block->{Maternal};
        my $paternal = $block->{Paternal};

        next
          if ( $maternal eq $metadata->{empty} );
        
        if ($paternal eq $metadata->{empty}) {
            $maternal_blocks{$maternal}{$block->{scaffold}}{$block->{start}} = $block->{end};
            next;
        }

        $linkage_groups{$maternal}{$paternal}{length} += $block->{length};
        $linkage_groups{$maternal}{$paternal}{blocks}++;
        $linkage_groups{$maternal}{$paternal}{output} = $paternal;
        if ( defined $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } ) {
            if ( $block->{end} > $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } ) {
                $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } = $block->{end};
            }
        }
        else {
            $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } = $block->{end};
        }
    }
    
    if ( $known ne '' ) {
        open my $knownfile, '<', $known or croak "Can't open known file $known: $OS_ERROR\n";
        while ( my $knownline = <$knownfile> ) {
            chomp $knownline;
            my ( $maternal, $paternal, $comment ) = split /\t/, $knownline;
            if ( !defined $linkage_groups{$maternal}{$paternal} ) {
                $linkage_groups{$maternal}{$paternal}{length} = 0;
                $linkage_groups{$maternal}{$paternal}{blocks} = 0;
                $linkage_groups{$maternal}{$paternal}{output} = $paternal;
            }
            $linkage_groups{$maternal}{$paternal}{known}++;
        }
        close $knownfile;
    }

    ( \%linkage_groups, \%marker_blocks, \%maternal_blocks );
}


sub build_initial_chromosome_map {
    my ( $lg, $chromosome_print, $output ) = @_;

    # First pass through markers in order of length
    my %revised_lg;
    my $marker_count;
    my %rejects;
    for my $marker (
        sort { defined $lg->{$a}{known} ? -1 : defined $lg->{$b}{known} ? 1 : $lg->{$b}{length} <=> $lg->{$a}{length} }
        keys %{$lg}
      )
    {
        $marker_count++;
        my ( $fixed_marker, $outcome ) =
          integrate_marker( $marker, $lg->{$marker}, \%revised_lg, $chromosome_print, $output );
        $lg->{$marker}{output} = $fixed_marker if $outcome eq 'Valid';
        $rejects{$outcome}{$marker} = $fixed_marker;
    }
    ( \%revised_lg, \%rejects );
}

sub revise_chromosome_map {
    my ( $lg, $revised_lg, $rejects, $chromosome_print, $output ) = @_;

    # Attempt to integrate rejected markers until no more are added on one pass
    my $more_added = 1;

    while ($more_added) {
        $more_added = 0;
        for my $outcome ( sort keys %{$rejects} ) {
            next if $outcome eq 'Valid';
            for my $marker (
                sort {
                        defined $lg->{$a}{known} ? -1
                      : defined $lg->{$b}{known} ? 1
                      : $lg->{$b}{length} <=> $lg->{$a}{length}
                } keys %{ $rejects->{$outcome} }
              )
            {
                my ( $fixed_marker, $new_outcome ) =
                  integrate_marker( $marker, $lg->{$marker}, $revised_lg, $chromosome_print, $output );
                if ( $new_outcome eq 'Valid' ) {
                    $more_added++;
                    $lg->{$marker}{output} = $fixed_marker;
                    delete $rejects->{$outcome}{$marker};
                    delete $rejects->{$outcome} if ( keys %{ $rejects->{$outcome} } == 0 );
                    $rejects->{'Valid'}{$marker} = $fixed_marker;
                }
            }
        }
    }
    return;
}

sub integrate_marker {
    my ( $marker, $marker_ref, $rlg, $chromosome_print, $output, $verbose ) = @_;
    $verbose = $verbose // 0;

    $rlg->{$marker}{output} = $marker;
    $rlg->{$marker}{length} = $marker_ref->{length};
    return ( $marker, 'Valid' ) if keys %{$rlg} < 3 or defined $marker_ref->{known} or $marker_ref->{length} >= 200000;

    # Build map with marker in it
    my $marker_id = run_mstmap( $chromosome_print, $rlg, $output );
    my $map = load_map( $output, $chromosome_print );

    my $chromosome_lgs = keys %{$map};
    if ( $chromosome_lgs >= 2 ) {
        phase_linkage_groups( $map, $rlg, $marker_id );
        $marker_id = run_mstmap( $chromosome_print, $rlg, $output );
        $map = load_map( $output, $chromosome_print );
        my $chromosome_lgs = keys %{$map};
    }

    my ( $fixed_marker, $outcome ) = check_marker_on_map( $marker, $map, $rlg, $marker_id, $verbose );

    delete $rlg->{$marker} if $outcome eq 'Disorder' or $outcome eq 'End' or $outcome eq 'Missing';

    if ( $outcome eq 'Valid' and $fixed_marker ne $marker ) {
        if ( !defined $rlg->{$fixed_marker} ) {
            $rlg->{$fixed_marker}{output} = $fixed_marker;
            $rlg->{$fixed_marker}{length} = $rlg->{$marker}{length};
        }
        else {
            $rlg->{$fixed_marker}{length} += $rlg->{$marker}{length};
        }
        delete $rlg->{$marker};
    }
    else {
    }

    ( $fixed_marker, $outcome );
}

sub check_marker_on_map {
    my ( $marker, $map, $markers, $marker_id, $verbose ) = @_;

    my $original_marker;
    my $new_marker    = '';
    my $marker_status = '';

    print "\n" if $verbose;
    for my $lg ( sort keys %{$map} ) {
        my @cMs = sort { $a <=> $b } keys %{ $map->{$lg} };
        for my $i ( 0 .. $#cMs ) {
            my $new          = 0;
            my $cM           = $cMs[$i];
            my $cM_marker_id = ( keys %{ $map->{$lg}{$cM} } )[0];
            my $cM_marker    = $marker_id->{$cM_marker_id};

            print "\n$lg\t$cM\t$markers->{$cM_marker}{output}\t$markers->{$cM_marker}{length}" if $verbose;

            if ( $cM_marker eq $marker ) {
                $new++;

                print "\tNew" if $verbose;
            }

            if ( $i == 0 or $i == $#cMs ) {

                print "\tEnd" if $verbose;
                $marker_status = 'End' if $new;
                next;
            }

            $original_marker = $markers->{$cM_marker}{output};

            my $this = $cM_marker;

            my $prev = $marker_id->{ ( keys %{ $map->{$lg}{ $cMs[ $i - 1 ] } } )[0] };
            my $next = $marker_id->{ ( keys %{ $map->{$lg}{ $cMs[ $i + 1 ] } } )[0] };

            my @prev = split //, $markers->{$prev}{output};
            my @this = split //, $markers->{$this}{output};
            my @next = split //, $markers->{$next}{output};

            my $fixed_marker = '';
            for my $i ( 0 .. $#this ) {
                if ( $prev[$i] eq $next[$i] and $prev[$i] ne $this[$i] ) {
                    $this[$i] = $prev[$i];
                }
            }
            $fixed_marker = join '', @this;
            $new_marker = $fixed_marker if $new;

            if ( $fixed_marker ne $original_marker ) {
                if ($new) {

                    print "\tFixed" if $verbose;
                }
                else {
                    print "\tDisordered" if $verbose;
                    $marker_status = 'Disorder';
                }
            }
        }
    }
    $marker_status = $marker_status ne '' ? $marker_status : $new_marker eq '' ? 'Missing' : 'Valid';

    print "\t$marker_status" if $verbose;
    ( $new_marker, $marker_status );
}

sub make_linkage_map {
    my ( $chromosome_print, $markers, $output ) = @_;

    # Make map for this maternal pattern
    my $marker_id = run_mstmap( $chromosome_print, $markers, $output );

    my $map = load_map( $output, $chromosome_print );
    my $chromosome_lgs = keys %{$map};

    #    print STDERR "\tBuilt $chromosome_lgs linkage groups";

    # If more than one linkage map returned for this chromosome,
    # attempt to rephase markers and remake the map
    my $phased_markers;
    if ( $chromosome_lgs >= 2 ) {
        phase_linkage_groups( $map, $markers, $marker_id );
        $marker_id = run_mstmap( $chromosome_print, $markers, $output );
        $map = load_map( $output, $chromosome_print );
        my $chromosome_lgs = keys %{$map};

        #        print STDERR "\tAfter phasing, $chromosome_lgs linkage groups";
    }

    #    $marker_id = run_mstmap( $chromosome_print, $markers, $output );
    #    $map = load_map( $output, $chromosome_print );

    my $map_markers = write_map_markers( $map, $markers, $marker_id, $output, $chromosome_print );

    #    print STDERR "\n";

    $map_markers;
}

sub run_mstmap {

    my ( $chromosome_print, $markers, $output ) = @_;

    my $mstmap_input = write_mstmap_header( $output, $chromosome_print, $markers );

    my $id = 1;
    my %marker_id;
    for my $marker (
        sort { $markers->{$a}{output} cmp $markers->{$b}{output} || $markers->{$b}{length} <=> $markers->{$a}{length} }
        keys %{$markers}
      )
    {
        $marker_id{$id} = $marker;
        $markers->{$marker}{output} =
          check_mirror( $markers->{$marker}{output}, $markers );
        $id = write_marker( $id, $markers->{$marker}, $mstmap_input );
    }

    close $mstmap_input;

    system(
"MSTMap.exe $output.$chromosome_print.mstmap.input $output.$chromosome_print.mstmap.map > $output.$chromosome_print.mstmap.log"
    );

    \%marker_id;
}

sub write_mstmap_header {
    my ( $output, $chromosome_print, $markers ) = @_;

    my %mst_header = (
        population_type              => "RIL2",
        population_name              => "HeliconiusWGS",
        distance_function            => "kosambi",
        cut_off_p_value              => "0.000001",
        no_map_dist                  => "0",
        no_map_size                  => "0",
        missing_threshold            => "1",
        estimation_before_clustering => "yes",
        detect_bad_data              => "yes",
        objective_function           => "ML",
    );

    my @mst_header = (
        "population_type", "population_name",    "distance_function", "cut_off_p_value",
        "no_map_dist",     "no_map_size",        "missing_threshold", "estimation_before_clustering",
        "detect_bad_data", "objective_function", "number_of_loci",    "number_of_individual",
    );

    open my $mstmap_input, ">", "$output.$chromosome_print.mstmap.input"
      or croak "Can't open $output.$chromosome_print.mstmap.input: $OS_ERROR\n";

    my $samples = split //, ( keys %{$markers} )[0];
    $mst_header{"number_of_loci"}       = keys %{$markers};
    $mst_header{"number_of_individual"} = $samples;
    map { print $mstmap_input "$_ $mst_header{$_}\n"; } @mst_header;

    print $mstmap_input "locus_name";
    map { print $mstmap_input "\t$_" } ( 1 .. $samples );
    print $mstmap_input "\n";

    $mstmap_input;
}

sub write_marker {
    my ( $id, $marker, $mstmap_input, $marker_lookup ) = @_;

    print $mstmap_input "$id";
    my @gt = split //, $marker->{output};
    map { print $mstmap_input "\t"; print $mstmap_input $_ eq 'H' ? 'X' : $_; } @gt;
    print $mstmap_input "\n";
    $id++;

    $id;
}

sub load_map {
    my ( $output, $chromosome_print ) = @_;

    open my $mstmap_output, '<', "$output.$chromosome_print.mstmap.map"
      or croak "Can't open map for $chromosome_print! $OS_ERROR\n";

    my %lg;
    my $lg_num;
    my $in_group = 0;
    while ( my $mst_line = <$mstmap_output> ) {
        chomp $mst_line;

        $lg_num = $1 if ( $mst_line =~ /^group (.+)$/ );
        if ( $mst_line eq ";BEGINOFGROUP" ) {
            $in_group = 1;
            next;
        }
        $in_group = 0 if ( $mst_line eq ";ENDOFGROUP" );
        if ($in_group) {
            if ( $mst_line =~ /^(.+)\t([\d\.]+)$/ ) {
                my $id = $1;
                my $cM = $2;
                $lg{$lg_num}{$cM}{$id}++;
            }
        }
    }

    # Delete linkage groups with a single marker
    map { delete $lg{$_} if keys %{ $lg{$_} } == 1 } keys %lg;
    close $mstmap_output;

    \%lg;
}

sub check_mirror {
    my ( $marker, $markers ) = @_;
    my $mirror     = mirror($marker);
    my @markerlist = map { $markers->{$_}{output} } keys %{$markers};
    my $markerh    = sum_hamming( $marker, \@markerlist );
    my $mirrorh    = sum_hamming( $mirror, \@markerlist );

    $markerh < $mirrorh ? $marker : $mirror;
}

sub sum_hamming {
    my ( $marker, $list ) = @_;
    my $h = 0;
    for my $l ( @{$list} ) {
        next if $l eq $marker;
        $h += hamming( $l, $marker );
    }

    $h;
}

sub phase_linkage_groups {
    my ( $map, $markers, $marker_id ) = @_;

    my $first_lg = ( keys %{$map} )[0];
    for my $lg ( keys %{$map} ) {
        for my $cM ( keys %{ $map->{$lg} } ) {
            for my $id ( keys %{ $map->{$lg}{$cM} } ) {
                my $marker = $marker_id->{$id};
                my $output = $markers->{$marker}{output};
                my $phased = $lg eq $first_lg ? mirror($output) : $output;
                $markers->{$marker}{output} = $phased;
            }
        }
    }
}


sub write_map_markers {
    my ( $map, $markers, $marker_id, $output, $chromosome_print ) = @_;

    my %chromosome;
    open my $map_marker_file, ">", "$output.$chromosome_print.mstmap.markers"
      or croak "Can't open $output.$chromosome_print.mstmap.markers: $OS_ERROR\n";

    print $map_marker_file "ID\tLG\tCM\tOriginal\tOutput\tLength\n";
    for my $lg ( sort keys %{$map} ) {
        for my $cM ( sort { $a <=> $b } keys %{ $map->{$lg} } ) {
            for my $id ( sort { $a <=> $b } keys %{ $map->{$lg}{$cM} } ) {
                my $marker = $marker_id->{$id};
                $chromosome{$lg}{$cM}{$marker}{length} =
                  $markers->{$marker}{length};
                $chromosome{$lg}{$cM}{$marker}{output} =
                  $markers->{$marker}{output};
                print $map_marker_file "$id\t$lg\t$cM";
                print $map_marker_file "\t$marker\t$markers->{$marker}{output}";
                print $map_marker_file "\t$markers->{$marker}{length}";
                print $map_marker_file "\n";
            }
        }
    }
    close $map_marker_file;

    \%chromosome;
}

sub output_map_summary {
    my ( $chromosome, $lg, $rejects ) = @_;
    my $summary       = sprintf "%2d: ", $chromosome;
    my $total_markers = 0;
    my $total_bases   = 0;
    my %total_uniques;
    for my $outcome ( "Disorder", "End", "Missing", "Valid" ) {
        my %unique_markers;
        my $marker_num = 0;
        my $bases      = 0;
        if ( defined $rejects->{$outcome} ) {
            $marker_num = keys %{ $rejects->{$outcome} };
            for my $marker ( keys %{ $rejects->{$outcome} } ) {
                $unique_markers{ $lg->{$marker}{output} }++;
                $total_uniques{ $lg->{$marker}{output} }++;
                $bases += $lg->{$marker}{length};
                delete $lg->{$marker} if $outcome ne 'Valid';
            }
        }
        $total_markers += $marker_num;
        $total_bases   += $bases;
        $summary .= sprintf "%8s:%4d:%4d:%9d ", $outcome, $marker_num, scalar keys %unique_markers, $bases;
    }
    $summary .= sprintf "   Total:%4d:%4d:%9d\n", $total_markers, scalar keys %total_uniques, $total_bases;
    print STDERR $summary;
    return;
}

sub output_rejected_markers {
    my ( $chromosome_print, $rejects, $lg, $output ) = @_;

    open my $rejected, '>', "$output.$chromosome_print.rejected.tsv"
      or croak "Can't open rejected markers file! $OS_ERROR\n";

    for my $outcome ( sort keys %{$rejects} ) {
        next if $outcome eq 'Valid';
        for my $marker ( keys %{ $rejects->{$outcome} } ) {
            print $rejected
              "$chromosomes{$chromosome_print}\t$lg->{$marker}{output}\t$lg->{$marker}{length}\t$outcome\n";
        }
    }
    close $rejected;
}

## MAP GENOME

sub output_marker_blocks {
    my ( $linkage_groups, $marker_blocks, $maternal_blocks, $blocklist, $input, $output ) = @_;

    my ( $dbh, $insert_chromosome, $insert_scaffold ) = create_map_tables($input);

    my $chromosome_numbers = get_chromosome_numbers($linkage_groups);

    my %filled_blocks;
    my %genome;
    for my $chromosome ( sort { $a <=> $b } keys %{$chromosome_numbers} ) {
        my $chromosome_print = $chromosome_numbers->{$chromosome};
        my $map = load_chromosome_map( $chromosome_print, \%genome, $output );

        for my $lg ( sort keys %{$map} ) {
            for my $cM ( sort { $a <=> $b } keys %{ $map->{$lg} } ) {
                for my $marker ( keys %{ $map->{$lg}{$cM} } ) {

                    $insert_chromosome->execute(
                        $chromosome, $chromosome_print, $cM, $marker,
                        colour_clean_marker( $map->{$lg}{$cM}{$marker}{output} ),
                        $map->{$lg}{$cM}{$marker}{length}
                    );

                    insert_scaffold_blocks( $marker_blocks->{$marker}, $chromosome, $cM, $insert_scaffold, \%filled_blocks );
                }
            }
        }
        insert_scaffold_blocks($maternal_blocks->{$chromosome_print}, $chromosome, -1, $insert_scaffold, \%filled_blocks );
    }
    
    for my $block ( @{$blocklist} ) {
        next if defined $filled_blocks{$block->{scaffold}}{$block->{start}};

        $insert_scaffold->execute( 0, -1, $block->{scaffold}, $block->{start}, $block->{end}, $block->{length} );
    }
    
    $dbh->commit;
    $dbh->disconnect;
}

sub create_map_tables {
    my ($input) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "", { AutoCommit => 0 } );

    $dbh->do("DROP TABLE IF EXISTS chromosome_map") or die $dbh->errstr;
    $dbh->do(
"CREATE TABLE chromosome_map (chromosome integer, print text, cm real, original text, clean text, length integer)"
    ) or die $dbh->errstr;
    $dbh->do("DROP TABLE IF EXISTS scaffold_map") or die $dbh->errstr;
    $dbh->do(
"CREATE TABLE scaffold_map (chromosome integer, cm real, scaffold text, start integer, end integer, length integer)"
    ) or die $dbh->errstr;

    $dbh->commit;

    my $insert_chromosome = $dbh->prepare_cached("INSERT INTO chromosome_map VALUES (?,?,?,?,?,?)");
    my $insert_scaffold   = $dbh->prepare_cached("INSERT INTO scaffold_map VALUES (?,?,?,?,?,?)");

    ( $dbh, $insert_chromosome, $insert_scaffold );
}

sub get_chromosome_numbers {
    my ($linkage_groups) = @_;

    my %chromosome_numbers;
    my $new_chromosome_number = 100;
    for my $print ( keys %{$linkage_groups} ) {
        my $chromosome = $chromosomes{$print} // $new_chromosome_number;
        $new_chromosome_number++ if $chromosome == $new_chromosome_number;
        $chromosome_numbers{$chromosome} = $print;
    }

    \%chromosome_numbers;
}

sub load_chromosome_map {
    my ( $chromosome_print, $genome, $output ) = @_;
    open my $chrom_file, '<', "$output.$chromosome_print.mstmap.markers"
      or croak "Can't open file for chromosome $chromosomes{$chromosome_print}! $OS_ERROR\n";
    my $header = <$chrom_file>;
    while ( my $chrom_line = <$chrom_file> ) {
        chomp $chrom_line;
        my ( $id, $lg, $cM, $original, $output, $length ) = split /\t/, $chrom_line;
        $genome->{$chromosome_print}{$lg}{$cM}{$original}{output} = $output;
        $genome->{$chromosome_print}{$lg}{$cM}{$original}{length} = $length;
    }
    close $chrom_file;

    croak "More than one linkage group found for chromosome print $chromosome_print"
      if keys %{ $genome->{$chromosome_print} } > 1;

    $genome->{$chromosome_print};
}

sub colour_clean_marker {
    my ($marker) = @_;

    my $colour_marker = "";
    for my $call ( split //, $marker ) {
        my $colour = $call eq 'A' ? 'red1' : $call eq 'B' ? 'royalblue1' : $call eq 'H' ? 'sandybrown' : 'bold';
        $colour_marker .= fg $colour, $call;
    }

    $colour_marker;
}

sub insert_scaffold_blocks {
    my ( $blocks, $chromosome, $cM, $insert_scaffold, $filled_blocks ) = @_;
    for my $scaffold ( sort keys %{$blocks} ) {
        for my $start ( sort { $a <=> $b } keys %{ $blocks->{$scaffold} } ) {
            my $end = $blocks->{$scaffold}{$start};
            my $length = $end - $start + 1;
            $insert_scaffold->execute( $chromosome, $cM, $scaffold, $start, $end, $length );
            $filled_blocks->{$scaffold}{$start} = 1
        }
    }
    return;
}

## STATISTICS

sub get_block_stats {
    my ( $step, $blocklist, $metadata ) = @_;

    my $stats = get_stats( $blocklist, $metadata );

    $metadata->{stats} = $stats;
    return if !$metadata->{verbose};

    for my $type ( @{ $metadata->{header} } ) {
        next if !defined $metadata->{types}{$type};
        output_type_stats( $step, $type, $stats ) if $metadata->{verbose};
    }

    return;
}

sub update_block_stats {
    my ( $step, $blocklist, $metadata ) = @_;

    my $stats = get_stats( $blocklist, $metadata );

    my $all_equal = 1;

    for my $type ( @{ $metadata->{header} } ) {
        next if !defined $metadata->{types}{$type};

        next if stats_equal( $metadata->{stats}{$type}, $stats->{$type} );

        $all_equal = 0;

        output_type_stats( $step, $type, $stats ) if $metadata->{verbose};
    }

    if ( $all_equal and $metadata->{verbose} ) {
        printf STDERR "%-60s", $step;
        print STDERR "No change\n";
    }

    $metadata->{stats} = $stats;

    return;
}

sub output_type_stats {
    my ( $step, $type, $stats ) = @_;

    printf STDERR "%-60s", $step;
    printf STDERR "%-25s", $type;
    for my $full ('ok') {
        printf STDERR "%12d", $stats->{$type}{$full}{patterns};
        printf STDERR "%12d", $stats->{$type}{$full}{blocks};
        printf STDERR "%12d", $stats->{$type}{$full}{bases};
    }
    print STDERR "\n";

    return;
}

sub get_stats {
    my ( $blocklist, $metadata ) = @_;

    my %stats;
    my $empty = $metadata->{empty};

    my %validities;
    for my $block ( @{$blocklist} ) {
        for my $type ( sort keys %{ $metadata->{types} } ) {

            next if $block->{$type} eq $empty;

            my $full = $block->{$type} =~ /\-/ ? 'no' : 'ok';

            $stats{$type}{$full}{blocks}++;
            $stats{$type}{$full}{bases} += $block->{'length'};
            $stats{$type}{$full}{patterns}{ $block->{$type} }++;
        }
        $validities{ $block->{validity} }{blocks}++;
        $validities{ $block->{validity} }{bases} += $block->{'length'};
        $validities{ $block->{validity} }{patterns}{ $block->{'Paternal'} }++;
    }
    for my $validity (
        sort { $validities{$b}{bases} <=> $validities{$a}{bases} }
        keys %validities
      )
    {
        printf STDERR "%s\t%5d\t%5d\t%10d\n", $validity, $validities{$validity}{blocks},
          scalar keys %{ $validities{$validity}{patterns} },
          $validities{$validity}{bases};
    }
    for my $type ( keys %{ $metadata->{types} } ) {
        for my $full ( 'ok', 'no' ) {
            $stats{$type}{$full}{blocks} = $stats{$type}{$full}{blocks} // 0;
            $stats{$type}{$full}{bases}  = $stats{$type}{$full}{bases}  // 0;
            $stats{$type}{$full}{patterns} =
              defined $stats{$type}{$full}{patterns}
              ? scalar keys %{ $stats{$type}{$full}{patterns} }
              : 0;
        }
    }

    \%stats;
}

sub stats_equal {
    my ( $oldstats, $newstats ) = @_;
    my $stats_equal = 1;

    for my $full ( 'ok', 'no' ) {
        $stats_equal = 0
          if $oldstats->{$full}{patterns} ne $newstats->{$full}{patterns};
        $stats_equal = 0
          if $oldstats->{$full}{blocks} ne $newstats->{$full}{blocks};
        $stats_equal = 0
          if $oldstats->{$full}{bases} ne $newstats->{$full}{bases};
    }

    $stats_equal;
}

## LIBRARY FUNCTIONS


sub phase {
    my $pat = shift;
    my @checkbases;
    if ( $pat =~ 'A' && $pat =~ 'B' && $pat =~ 'H' ) {
        @checkbases = ( 'A', 'B' );
    }
    elsif ( $pat =~ 'A' && $pat =~ 'B' && $pat !~ 'H' ) {
        @checkbases = ( 'A', 'B' );
    }
    elsif ( $pat =~ 'B' && $pat =~ 'H' && $pat !~ 'A' ) {
        @checkbases = ( 'B', 'H' );
    }
    return " " x length($pat) if !@checkbases;

    my @p = split //, $pat;
    my $out = "";
    my %trans;

    for my $gt (@p) {
        my $added = 0;
        for my $bi ( 0, 1 ) {
            if ( $gt eq $checkbases[$bi] ) {
                if ( !defined $trans{$gt} ) {
                    $trans{ $checkbases[$bi] } = 'A';
                    my $other = $bi == 0 ? 1 : 0;
                    $trans{ $checkbases[$other] } = 'B';
                }
                $out .= $trans{$gt};
                $added++;
            }
        }
        $out .= $gt =~ /[ ~\.]/ ? '-' : $gt if !$added;
    }

    $out;
}

sub mirror {
    my $pat = shift;
    if (   ( $pat =~ 'A' and $pat =~ 'B' and $pat =~ 'H' )
        or ( $pat =~ 'A' and $pat =~ 'B' and $pat !~ 'H' ) )
    {
        $pat =~ tr/AB/BA/;
    }
    if ( $pat =~ 'A' and $pat =~ 'H' and $pat !~ 'B' ) {
        $pat =~ tr /AH/HA/;
    }
    if ( $pat =~ 'B' and $pat =~ 'H' and $pat !~ 'A' ) {
        $pat =~ tr /HB/BH/;
    }

    $pat;
}

sub hamming {
    return ( $_[0] ^ $_[1] ) =~ tr/\001-\255//;
}

sub get_header_types {
    my $dbh = shift;

    my $sth = $dbh->prepare("SELECT * FROM blocks WHERE scaffold = 'getcolumnnames'");
    $sth->execute;
    my @header = @{ $sth->{NAME} };

    my %types;
    map { $types{$_} = "" } @header;
    delete $types{"scaffold"};
    delete $types{"start"};
    delete $types{"end"};
    delete $types{"length"};

    ( \@header, \%types );
}
