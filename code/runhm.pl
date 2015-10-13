#!/usr/bin/env perl

# runhm.pl
# Run HaploMerger

# John Davey
# johnomics@gmail.com

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use File::Basename 'fileparse';
use File::Copy 'copy';
use File::Copy::Recursive 'dircopy';
use File::Path 'rmtree';
use IO::Uncompress::Gunzip qw($GunzipError);
use Cwd;

# Autoflush output so reporting on progress works
$| = 1;

my $input           = "";
my $configdir       = "";
my $outputdir       = "";
my $prefix          = "";
my $scaffold_prefix = "";
my $badmerges       = "";
my $annotation      = "";
my $breakgenes;
my $g;

my $options_okay = GetOptions(
    'input=s'           => \$input,
    'configdir=s'       => \$configdir,
    'outputdir=s'       => \$outputdir,
    'prefix=s'          => \$prefix,
    'scaffold_prefix=s' => \$scaffold_prefix,
    'badmerges=s'       => \$badmerges,
    'annotation=s'      => \$annotation,
    'x'                 => \$breakgenes,
    'g'                 => \$g,
);

croak "No output directory! Please specify -o $OS_ERROR\n" if $outputdir eq "";

if ( !-d $outputdir ) {
    croak "No config directory! Please specify -c $OS_ERROR\n" if $configdir eq "";
    croak "Config directory $configdir does not exist!\n" if !-d $configdir;
    dircopy $configdir, $outputdir;
}

if ( !-e "$outputdir/genome.fa" ) {
    croak "No FASTA file! Please specify -f $OS_ERROR\n" if $input eq "";
    croak "FASTA file $input does not exist!\n" if !-e $input;
    copy $input, "$outputdir/genome.fa";
}

my %badmerges;
if ($badmerges) {
    if ( -e $badmerges ) {
        open my $bmfh, '<', $badmerges or croak "Can't open bad merges file $badmerges $OS_ERROR\n";
        while ( my $badmergeline = <$bmfh> ) {
            chomp $badmergeline;
            my ( $a, $b ) = split "\t", $badmergeline;
            $badmerges{$a}{$b} = 1;
            $badmerges{$b}{$a} = 1;
        }
        close $bmfh;
    }
    else {
        croak "Can't find badmerges file $badmerges! $OS_ERROR\n";
    }
}

my %genes;
if ( $annotation and -e $annotation ) {
    open my $gff, '<', $annotation or croak "Can't open annotation file $annotation $OS_ERROR\n";
    while ( my $featureline = <$gff> ) {
        chomp $featureline;
        my ( $scaffold, $source, $featuretype, $start, $end, $score, $strand, $phase, $attributes ) = split "\t",
          $featureline;
        if ( $featuretype eq 'gene' ) {
            $genes{$scaffold}{$start}{end}        = $end;
            $genes{$scaffold}{$start}{attributes} = $attributes;
        }
    }
    close $gff;
}

my ( $filename, $dirs, $suffix ) = fileparse( $input, qr/\.[^.]*/ );
chdir $outputdir;

my $log;
if ( !-e "runhm.log" ) {
    open $log, '>', "runhm.log";
}
else {
    open $log, '>>', "runhm.log";
}

print $log "HaploMerger run log\n";
print $log "Genome: $input\n";
print $log "Config directory: $configdir\n";
print $log "Output directory: $outputdir\n";
print $log "Prefix: '$prefix'\n";
print $log "Scaffold prefix: '$scaffold_prefix'\n";
print $log "Run G: ";
print $log ( $g ? "Yes" : "No" ), "\n";

my $start_time = localtime;
my $start      = time;
print $log "Start time: $start_time\n";

my $optimized_file = "optiNewScaffolds.fa.gz";
my $unpaired_file  = "unpaired.fa.gz";

( $optimized_file, $unpaired_file, $log, my $next_start ) = run_abc( $log, $start );

( $optimized_file, $unpaired_file, $log, $next_start ) = run_e( $log, $next_start, \%badmerges, \%genes, $breakgenes )
  if $badmerges ne "" or %genes;

( $optimized_file, $unpaired_file, $log, $next_start ) = run_f( $log, $next_start, $scaffold_prefix, \%genes )
  if $scaffold_prefix ne "" or %genes;

output_final_genome( $outputdir, $optimized_file, $unpaired_file, $prefix );

if ($g) {
    ( $optimized_file, $unpaired_file, $log, $next_start ) = run_g( $log, $next_start );
    output_final_genome( $outputdir, $optimized_file, $unpaired_file, $prefix, "refined" );
}

open $log, '>>', "runhm.log";

if ( -d 'genome.genomex.result/raw.axt' ) {
    printf $log "Removing raw.axt folder\n";
    rmtree ['genome.genomex.result/raw.axt'];
}

my $end_time = localtime;
print $log "End time: $end_time\n";

output_duration( "Total", $log, $start );
print $log "Done\n";
close $log;

sub run_abc {
    my ( $log, $start ) = @_;
    my $c_start = $start;
    if ( !-e "genome.genomex.result/mafFiltered.net.maf.tar.gz" ) {
        system "./hm.batchA.initiation_and_all_lastz genome > runhm.out 2>&1";
        my $b_start = output_duration( "A", $log, $start );
        system "./hm.batchB.chainNet_and_netToMaf genome >> runhm.out 2>&1";
        $c_start = output_duration( "B", $log, $b_start );
    }

    my $next_start = $c_start;
    if ( !-e "genome.genomex.result/optiNewScaffolds.fa.gz" ) {
        system "./hm.batchC.haplomerger genome >> runhm.out 2>&1";
        $next_start = output_duration( "C", $log, $c_start );
    }

    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start );
}

sub run_e {
    my ( $log, $start, $badmerges, $genes, $breakgenes ) = @_;
    my $next_start = $start;

    if ( !-e "genome.genomex.result/hm.nodes_edited" ) {
        my $gene_boundaries = find_broken_genes($genes);
        edit_nodes($badmerges, $genes, $gene_boundaries, $breakgenes);
        chdir "genome.genomex.result";
        system "ln -s hm.sc_portions hm.sc_portions_edited"         if !-l "hm.sc_portions_edited";
        system "ln -s hm.assembly_errors hm.assembly_errors_edited" if !-l "hm.assembly_errors_edited";
        chdir "..";

        if ( -e "genome.genomex.result/hm.nodes_edited" ) {
            system "./hm.batchE.refine_haplomerger_updating_nodes_and_portions genome >> runhm.out 2>&1";
            $next_start = output_duration( "E", $log, $next_start );
        }
    }
    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start );
}

sub run_f {
    my ( $log, $start, $scaffold_prefix, $genes ) = @_;
    my $next_start = $start;

    if ( !-e "genome.genomex.result/hm.new_scaffolds_edited" ) {
        edit_new_scaffolds( $scaffold_prefix, $genes );
        if ( -e "genome.genomex.result/hm.new_scaffolds_edited" ) {
            system "./hm.batchF.refine_haplomerger_connections_and_Ngap_fillings genome >> runhm.out 2>&1";
            $next_start = output_duration( "F", $log, $next_start );
        }
    }

    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start );
}

sub find_broken_genes {
    my ($genes) = @_;

    open my $new_scaffolds, '<', "_hm.haploMerger.log"
      or croak "Can't open hm.haploMerger.log!\n";

    my %gene_boundaries;
    while ( my $line = <$new_scaffolds> ) {
        chomp $line;
        next if $line =~ /^#/ or $line =~ /^$/;
        my @f = split "\t", $line;
        next if @f != 6;
        
        my ($scaffold, $size, $hapid, $start, $end, $strand) = @f;

        for my $genestart (keys %{$genes->{$scaffold}}) {
            my $gene = $genes->{$scaffold}{$genestart};
            if ($genestart < $start and $gene->{end} >= $start) {
                $gene_boundaries{$scaffold}{$start} = $gene;
            }
            if ($genestart <= $end and $gene->{end} > $end) {
                $gene_boundaries{$scaffold}{$end} = $gene;
            }
        }
    }
    return \%gene_boundaries;
}

sub edit_nodes {
    my ($badmerges, $genes, $gene_boundaries, $breakgenes) = @_;

    my $nodes_file = "genome.genomex.result/hm.nodes";
    open my $nodes, '<', $nodes_file or croak "Can't open nodes file!\n";
    open my $edited, '>', "genome.genomex.result/hm.nodes_edited" or croak "Can't open edited nodes file!\n";

    while ( my $portion = <$nodes> ) {
        if ( $portion =~ /^#/ or $portion =~ /^$/ ) {
            print $edited $portion;
            next;
        }
        chomp $portion;
        my @f = split "\t", $portion;
        my $portion1 = get_node_portion(\@f, 1);
        my $portion2 = get_node_portion(\@f, 2);
        if (   ( defined $badmerges->{ $portion1->{scaffold} } and defined $badmerges->{ $portion1->{scaffold} }{ $portion2->{scaffold} } )
            or ( defined $badmerges->{ $portion2->{scaffold} } and defined $badmerges->{ $portion2->{scaffold} }{ $portion1->{scaffold} } ) )
        {
            $f[-1] = 1;
        }

        my $curated1 = check_curated($portion1, $genes);
        my $curated2 = check_curated($portion2, $genes);
        if ($curated1 and $curated2 and $portion1->{scaffold} ne $portion2->{scaffold}) {
            $f[-1] = 1;
        }
        else {
            my ($broken1, $broken_curated1) = check_broken($portion1, $gene_boundaries);
            my ($broken2, $broken_curated2) = check_broken($portion2, $gene_boundaries);
            # Gene is broken; fix it by removing this node if it is curated or if x option (breakgenes) is off
            if ($broken_curated1 or (not $breakgenes and $broken1) or $broken_curated2 or (not $breakgenes and $broken2)) {
                $f[-1] = 1;
            }            
        }

        my $edit = join "\t", @f;
        $edit .= "\n";

        print $edited $edit;
    }
    close $nodes;
    close $edited;
}

sub get_node_portion {
    my ($f, $p) = @_;
    my %region;
    if ($p == 1) {
        $region{scaffold} = $f->[0];
        $region{start}    = $f->[5];
        $region{length}   = $f->[6];
        $region{end}      = $f->[7];
    }
    elsif ($p == 2) {
        $region{scaffold} = $f->[1];
        $region{start}    = $f->[10];
        $region{length}   = $f->[11];
        $region{end}      = $f->[12];
    }
    return \%region;
}

sub check_broken {
    my ($portion, $gene_boundaries) = @_;
    my $scfgenes = $gene_boundaries->{$portion->{scaffold}};
    my $broken = 0;
    my $broken_curated = 0;
    if (defined $scfgenes) {
        for my $breakpoint (keys %{$scfgenes}) {
            my $gene = $scfgenes->{$breakpoint};
            if ($portion->{start} == $breakpoint or $portion->{end} == $breakpoint) {
                $broken = 1;
                $broken_curated = 1 if $gene->{attributes} =~ /Description/;
            }
        }
    }
    return $broken, $broken_curated;
}

sub edit_new_scaffolds {
    my ( $scaffold_prefix, $genes ) = @_;

    my $new_scaffolds_file = "genome.genomex.result/hm.new_scaffolds";
    if ( -e "genome.genomex.result/hm.new_scaffolds_updated" ) {    # If batchE has been run
        $new_scaffolds_file = "genome.genomex.result/hm.new_scaffolds_updated";
    }

    open my $new_scaffolds, '<', $new_scaffolds_file
      or croak "Can't open new scaffolds file!\n";
    open my $edited, '>', "genome.genomex.result/hm.new_scaffolds_edited"
      or croak "Can't open edited new scaffolds file!\n";

    my $scfnum = 0;
    my $curid  = -1;

    my %portions;
    while ( my $portion = <$new_scaffolds> ) {
        chomp $portion;
        my @f = ( $portion =~ /^#/ or $portion =~ /^$/ ) ? ($portion) : split "\t", $portion;
        if ( @f == 1 ) {
            print $edited "$f[0]\n";
            next;
        }
        $curid = $f[0] if $curid == -1;
        if ( $f[0] != $curid ) {
            $scfnum++;
            $curid = $f[0];
        }
        push @{ $portions{$scfnum} }, \@f;
    }

    for my $scfnum ( sort { $a <=> $b } keys %portions ) {

        for my $f ( @{ $portions{$scfnum} } ) {

            my $scaffold1 = $f->[5];
            my $scaffold2 = $f->[12];

            my $active_portion = 0;
            if ( $scaffold_prefix ne "" and $scaffold1 ne '0' and $scaffold2 ne '0' ) {
                if ( $scaffold1 =~ /^$scaffold_prefix/ and $scaffold2 !~ /^$scaffold_prefix/ ) {
                    $active_portion = 2;
                }
                elsif ( $scaffold1 !~ /^$scaffold_prefix/ and $scaffold2 =~ /^$scaffold_prefix/ ) {
                    $active_portion = 1;
                }
            }

            $f->[-2] = $active_portion if $active_portion;
        }
        if ($genes) {
            preserve_genes( \@{ $portions{$scfnum} }, $genes );
        }

        for my $f ( @{ $portions{$scfnum} } ) {
            my $edit = join "\t", @{$f};
            print $edited "$edit\n";
        }
    }
    close $new_scaffolds;
    close $edited;
}

sub preserve_genes {
    my ( $portions, $genes ) = @_;
    for my $f ( @{$portions} ) {
        my $portion1 = get_ns_portion($f, 1);
        my $portion2 = get_ns_portion($f, 2);
        my $curated1 = check_curated($portion1, $genes);
        my $curated2 = check_curated($portion2, $genes);
        if ($curated1 or $curated2) {
            if ($curated1 and $portion2->{active} or $curated2 and $portion1->{active}) { # Curated portion is haplotype; switch portions
                $f->[-2] = $portion1->{active} ? 2 : 1;
            }
        }
        else {
            my $broken1 = check_ns_broken($portion1, $genes);
            my $broken2 = check_ns_broken($portion2, $genes);
            $f->[-2] = 1 if $broken1 and !$broken2;
            $f->[-2] = 2 if $broken2 and !$broken1;
        }
    }
}

sub check_curated {
    my ($portion, $genes) = @_;
    my $scfgenes = $genes->{$portion->{scaffold}};
    my %foundgenes;
    if ( defined $scfgenes ) {
        for $start ( keys %{ $scfgenes } ) {
            next if $scfgenes->{$start}{attributes} !~ /Description/;
            my $end = $scfgenes->{$start}{end};
            if ($portion->{start} <= $start and $end <= $portion->{end}) {
                $foundgenes{$start}++;
            }
            elsif ($start < $portion->{start} and $portion->{start} <= $end or $start <= $portion->{end} and $portion->{end} < $end) {
                $foundgenes{$start}++;
            }
        }
    }
    return scalar keys %foundgenes;
}

sub check_ns_broken {
    my ($portion, $genes) = @_;
    my $scfgenes = $genes->{$portion->{scaffold}};
    my $broken = 0;
    if (defined $scfgenes) {
        for my $start (keys %{$scfgenes}) {
            my $end = $scfgenes->{$start}{end};
            if ($start < $portion->{start} and $portion->{start} <= $end or $start <= $portion->{end} or $portion->{end} < $end) {
                $broken = 1;
            }
        }
    }
    return $broken;
}

sub get_ns_portion {
    my ( $f, $portion ) = @_;
    my $active_portion = $f->[-2] ? $f->[-2] : $f->[4];
    my %region;
    $region{active} = $active_portion == $portion ? 1 : 0;
    if ( $portion == 1 ) {
        $region{scaffold} = $f->[5];
        $region{start}    = $f->[8] + 1;
        $region{end}      = $f->[9];
        $region{dir}      = $f->[10];
    }
    elsif ( $portion == 2 ) {
        $region{scaffold} = $f->[12];
        $region{start}    = $f->[15] + 1;
        $region{end}      = $f->[16];
        $region{dir}      = $f->[17];
    }
    return \%region;
}



sub run_g {
    my ( $log, $start ) = @_;
    my $next_start = $start;

    if ( !-e "genome.genomex.result/unpaired_refined.fa.gz" ) {
        system "./hm.batchG.refine_unpaired_sequences genome >> runhm.out 2>&1";
        $next_start = output_duration( "G", $log, $start );
    }
    ( "optiNewScaffolds.fa.gz", "unpaired_refined.fa.gz", $log, $next_start );
}

sub output_final_genome {
    my ( $outputdir, $optimized_file, $unpaired_file, $prefix, $suffix ) = @_;

    chdir "genome.genomex.result";

    my $outname = $suffix ? "$outputdir\_$suffix" : "$outputdir";

    if ( !-e "$outname.fa" ) {
        open my $finalgenome, '>', "$outname.fa";

        output_optimized_genome( $optimized_file, $finalgenome, $prefix, $outname );

        output_unpaired_genome( $unpaired_file, $finalgenome, $prefix, $outname );

        close $finalgenome;

        system "summarizeAssembly.py $outname.fa > $outname.summary";

        system "agp_from_fasta.py -f $outname.fa";
    }
    chdir "..";
}

sub output_optimized_genome {
    my ( $optimized_file, $finalgenome, $prefix, $outname ) = @_;

    if ( !-e $optimized_file ) {
        print $log "$optimized_file does not exist, abandon writing to final genome\n";
        return;
    }

    print $log "Writing $optimized_file to $outname.fa\n";
    my $opti = IO::Uncompress::Gunzip->new($optimized_file)
      or die "IO::Uncompress::Gunzip failed to open $optimized_file: $GunzipError\n";

    while ( my $fastaline = <$opti> ) {
        if ( $fastaline =~ /^>(.+) (.+)$/ ) {
            print $finalgenome ">$prefix$1 $2\n";
        }
        else {
            print $finalgenome $fastaline;
        }
    }

    close $opti;

}

sub output_unpaired_genome {
    my ( $unpaired_file, $finalgenome, $prefix, $outname ) = @_;

    if ( !-e $unpaired_file ) {
        print $log "$unpaired_file does not exist, abandon writing to final genome\n";
        return;
    }

    print $log "Writing $unpaired_file to $outname.fa\n";
    my $unpaired = IO::Uncompress::Gunzip->new($unpaired_file)
      or die "IO::Uncompress::Gunzip failed to open $unpaired_file: $GunzipError\n";

    while ( my $fastaline = <$unpaired> ) {
        if ( $fastaline =~ />(.+) old_(.+?);(.+)/ ) {
            print $finalgenome ">$prefix$1 old_$2;$3\n";
        }
        else {
            print $finalgenome $fastaline;
        }
    }

    close $unpaired;

}

sub output_duration {
    my ( $stage, $file, $start ) = @_;
    my $end      = time;
    my $duration = $end - $start;

    printf $file "$stage run time: %02d:%02d:%02d\n", int( $duration / 3600 ), int( ( $duration % 3600 ) / 60 ),
      int( $duration % 60 );

    $end;
}
