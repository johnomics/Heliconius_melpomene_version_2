#!/usr/bin/env perl

# vcf-parallel.pl
# John Davey johnomics@gmail.com

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Parallel::ForkManager;

$OUTPUT_AUTOFLUSH = 1;

my %args;

$args{vcf_filename}     = "";
$args{vcf_subs}         = "";
$args{threads}          = 1;
$args{output_prefix}    = "vcf-parallel-output";
$args{lengths_filename} = "";
$args{byscf}            = 0;
$args{maxscf}           = 0;
$args{uniquescf}        = "";
$args{extraargs}        = "";

my $options_okay = GetOptions(
    'vcf=s'       => \$args{vcf_filename},
    'script=s'    => \$args{vcf_subs},
    'threads=i'   => \$args{threads},
    'output=s'    => \$args{output_prefix},
    'lengths=s'   => \$args{lengths_filename},
    'byscaffold'  => \$args{byscf},
    'maxscf=i'    => \$args{maxscf},
    'uniquescf=s' => \$args{uniquescf},
    'extraargs=s' => \$args{extraargs},
);
croak "Can't process general options: $OS_ERROR\n" if !$options_okay;
croak "No VCF file! Please specify -v $OS_ERROR\n"
  if ( $args{vcf_filename} eq "" );

our $setup;
our $process;
our $merge;
our $output;

if ( $args{vcf_subs} eq "" ) {
    print STDERR "No script file given with -s, so will count lines in file\n";

    $setup = sub {
        my $args = shift;

        open my $vcf_file, '<', $args->{vcf_filename}
          or croak "Can't open VCF file $args->{vcf_filename}! $OS_ERROR\n";

        my %samples;

        my $vcf_line;
        while ( $vcf_line = <$vcf_file> ) {
            if ( $vcf_line =~ /^#CHROM/ ) {
                last;
            }
        }
        close $vcf_file;

        chomp $vcf_line;
        my @sample_names = split /\t/, $vcf_line;
        map {
            $samples{offspring}{lookup}{ $sample_names[$_] } = $_;
            push @{ $samples{offspring}{order} }, $sample_names[$_];
        } 9 .. $#sample_names;

        \%samples;
    };

    $process = sub {
        my ( $scf, $lines, $samples, $scfdata, $userdata ) = @_;
        if ( ref($lines) eq "ARRAY" ) {
            $scfdata->{lines} += @{$lines};
        }
        else {
            $scfdata->{lines}++;
        }
    };

    $merge = sub {
        my ( $part, $all ) = @_;
        return if !defined $part;
        $all->{lines} += $part->{lines};
    };

    $output = sub {
        my ( $data, $genome, $outfix ) = @_;

        print "$data->{lines} lines in VCF file\n";
    };
}
else {
    require( $args{vcf_subs} );
}

my $genome = load_genome( \%args );

$genome->{scfp} = get_partitions( $genome, $args{threads} );
$args{threads} = keys %{ $genome->{scfp} };    # last thread may not be used;
                                               # each part is always larger than minimum size, so scaffolds for final
                                               # thread may be spread around the other threads

print_partitions( $genome->{scfp}, $genome->{scfl} );

my $userdata = $setup->( \%args );

$output->( parse_vcf( $genome, $userdata, \%args ), $genome, $args{output_prefix}, $userdata );

sub parse_vcf {
    my ( $genome, $userdata, $argref ) = @_;
    my %data;

    # If only one thread specified, do not fork thread
    # (This allows code to be profiled)
    if ( $argref->{threads} == 1 ) {
        my ( $partdata, $scfi ) = run_part( 1, $genome, $userdata, $argref );

        print STDERR "1:Done, processed $scfi scaffolds of " . keys( %{ $genome->{scfp}{1} } ) . "\n";
        $merge->( $partdata, \%data );
        print STDERR "1:Merged\n";
        return \%data;
    }

    my $part_pm = new Parallel::ForkManager( $argref->{threads} );
    $part_pm->set_max_procs( $argref->{threads} );

    $part_pm->run_on_finish(
        sub {
            my $pid       = shift;
            my $exit      = shift;
            my $part      = shift;
            my $childdata = pop;
            $merge->( $childdata, \%data );
            print STDERR "$part:Merged\n";
        }
    );

    foreach my $part ( 1 .. $argref->{threads} ) {
        $part_pm->start($part) and next;

        my ( $partdata, $scfi ) = run_part( $part, $genome, $userdata, $argref );

        print STDERR "$part:Done, processed $scfi scaffolds of " . keys( %{ $genome->{scfp}{$part} } ) . "\n";
        $part_pm->finish( 0, $partdata );
    }
    $part_pm->wait_all_children;

    \%data;
}

sub run_part {
    my ( $part, $genome, $userdata, $argref ) = @_;

    open my $vcf_file, '<', $argref->{vcf_filename}
      or croak "Can't open $argref->{vcf_filename} $OS_ERROR!\n";

    # Skip header
    while (<$vcf_file>) {
        last if (/^#CHROM/);
    }

    my $curscf    = "";
    my $scfi      = 1;
    my $foundpart = 0;
    my %data;
    my @scf_vcf;
    my $scf;
    my $process_this_scf;
    while ( my $vcf_line = <$vcf_file> ) {
        last if ( $argref->{maxscf} && $scfi >= $argref->{maxscf} );
        $scf = $vcf_line =~ /^(.+?)\t/ ? $1 : "";

        $curscf = $scf if $curscf eq "";    # Fill curscf at the start

        $process_this_scf = ( $argref->{uniquescf} eq "" or $curscf eq $argref->{uniquescf} ) ? 1 : 0;

        if ( defined $genome->{scfp}{$part}{$scf} ) {
            $foundpart = 1;

            if ( $curscf ne $scf ) {
                if ( $argref->{byscf} && @scf_vcf ) {
                    $process->( $curscf, \@scf_vcf, \%data, $userdata, $genome->{scfl} ) if $process_this_scf;
                    @scf_vcf = ();
                }
                print_part_progress( $scfi, $part, $genome->{scfp} )
                  if ( $scfi % 10 == 0 );
                $scfi++;

                last if $curscf eq $argref->{uniquescf};
                $curscf = $scf;
                $process_this_scf = ( $argref->{uniquescf} eq "" or $curscf eq $argref->{uniquescf} ) ? 1 : 0;
            }
            $argref->{byscf}
              ? push @scf_vcf, $vcf_line
              : $process->( $curscf, $vcf_line, \%data, $userdata, $genome->{scfl} )
              if $process_this_scf;
        }
        else {
            if ($foundpart) {
                $process->( $curscf, \@scf_vcf, \%data, $userdata, $genome->{scfl} )
                  if $argref->{byscf}
                  && $process_this_scf;
                @scf_vcf = ();
                $scfi--;    # Last iteration will add an extra scaffold to the count,
                            #  so remove it here, after part has been processed
                last;
            }
        }
    }
    close $vcf_file;
    $process->( $curscf, \@scf_vcf, \%data, $userdata, $genome->{scfl} )
      if $argref->{byscf}
      && @scf_vcf
      && $process_this_scf;

    return ( \%data, $scfi );
}

sub print_part_progress {
    my ( $scfi, $part, $scfpref ) = @_;

    printf STDERR "%3d:%4d scaffolds processed, %4d remaining\n", $part, $scfi, keys( %{ $scfpref->{$part} } ) - $scfi;

}

sub print_partitions {
    my ( $scfpref, $scflref ) = @_;

    my $genscf   = 0;
    my $genpartl = 0;
    print STDERR "Part\tScaffolds\tLength\n";
    foreach my $part ( sort { $a <=> $b } keys %{$scfpref} ) {
        my $partl  = 0;
        my $numscf = keys %{ $scfpref->{$part} };
        foreach my $scf ( keys %{ $scfpref->{$part} } ) {
            $partl += $scflref->{$scf};
        }
        print STDERR "$part\t$numscf\t$partl\n";
        $genscf   += $numscf;
        $genpartl += $partl;
    }

    print STDERR "Genome\t$genscf\t$genpartl\n";
    return;
}

sub get_partitions {
    my ( $genome, $threads ) = @_;

    my %scfp;
    my $part      = 1;
    my $threshold = $genome->{length} / $threads;
    my $part_size = 0;
    for my $scf ( @{ $genome->{scf} } ) {
        $scfp{$part}{$scf}++;
        $part_size += $genome->{scfl}{$scf};
        if ( $part_size > $threshold ) {
            $part_size = 0;
            $part++;
        }
    }

    \%scfp;
}

sub load_genome {
    my $argref = shift;

    my @scflines =
         parse_vcf_header( $argref->{vcf_filename} )
      or load_lengths( $argref->{lengths_filename} )
      or croak "No scaffold lengths in VCF header or scaffold lengths file (-l)\n";
    return load_scflen(@scflines);
}

sub load_scflen {
    my @scflines = @_;

    my %scfl;
    my @scf;
    my $genl = 0;
    foreach my $scfline (@scflines) {
        my ( $scf, $len ) = split /\t/, $scfline;
        $scfl{$scf} = $len;
        push @scf, $scf;
        $genl += $len;
    }
    { scfl => \%scfl, scf => \@scf, length => $genl };
}

sub parse_vcf_header {
    my $vcf_filename = shift;
    open my $vcf_file, '<', $vcf_filename
      or croak "Can't open VCF file $vcf_filename! $OS_ERROR\n";

    my @scflines;
    my $vcf_line = <$vcf_file>;
    while ( $vcf_line !~ /^#CHROM/ ) {
        if ( $vcf_line =~ /^##contig=<ID=(.+),length=(\d+)>$/ ) {
            push @scflines, "$1\t$2";
        }
        $vcf_line = <$vcf_file>;
    }
    close $vcf_file;

    @scflines;
}

sub load_lengths {
    my $lengths_filename = shift;

    my @scflines;

    return @scflines if ( $lengths_filename eq "" );

    open my $lengths_file, "<", $lengths_filename
      or croak "Can't open scaffold lengths file $lengths_filename! $OS_ERROR\n";

    while ( my $scf_line = <$lengths_file> ) {
        if ( $scf_line =~ /^(.+)\t(.+)$/ ) {
            push @scflines, "$1\t$2";
        }
    }
    close $lengths_file;

    @scflines;
}

