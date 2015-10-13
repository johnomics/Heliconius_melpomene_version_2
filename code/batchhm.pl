#!/usr/bin/env perl

# batchhm.pl
# Run HaploMerger iteratively until pure haploid genome is produced

# John Davey
# johnomics@gmail.com

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
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
my $tsv             = "";
my $rounds          = 0;
my $g;
my $gff = "";

my $options_okay = GetOptions(
    'input=s'           => \$input,
    'configdir=s'       => \$configdir,
    'outputdir=s'       => \$outputdir,
    'prefix=s'          => \$prefix,
    'scaffold_prefix=s' => \$scaffold_prefix,
    'badmerges=s'       => \$badmerges,
    'tsv=s'             => \$tsv,
    'rounds=i'          => \$rounds,
    'g'                 => \$g,
    'annotation=s'      => \$gff,
);

croak "Can't open badmerges file $badmerges!" if $badmerges and !-e $badmerges;
croak "No TSV file!" if $tsv eq "";

my $basename = fileparse $tsv, ".tsv";
if ($gff) {
    my $transfer_features_command = "transfer_features.py -m $tsv -o $basename -g $gff";
    print "$transfer_features_command\n";
    system $transfer_features_command;
}

my $i       = 1;
my $command = "runhm.pl -i $input -c $configdir -o $outputdir\_$i -p $prefix$i -g";
$command .= " -s $scaffold_prefix" if $scaffold_prefix;
$command .= " -b $badmerges"       if $badmerges;
$command .= " -a $basename.gff"             if $gff;


print "$command\n";
system $command;
my $map_merge_command = "map_merge.py -m $outputdir\_$i -p $prefix$i";
print "$map_merge_command\n";
system $map_merge_command;
my $transfer_merge_command = "transfer_merge.py -d $tsv -n $outputdir\_$i/genome.genomex.result/$outputdir\_$i\_old.tsv -o $outputdir\_$i/$outputdir\_$i\_orig.tsv";
print "$transfer_merge_command\n";
system $transfer_merge_command;
my $transfer_features_command = "transfer_features.py -m $outputdir\_$i/$outputdir\_$i\_orig.tsv -o $outputdir\_$i/$outputdir\_$i";
$transfer_features_command .= " -g $gff" if $gff;
print "$transfer_features_command\n";
system $transfer_features_command;

while ( -e "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa"
    and -s "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa" > 0 )
{
    my $j        = $i + 1;
    my $filename = "$outputdir\_$i/genome.genomex.result/$outputdir\_$i";
    $filename .= $g ? "_refined.fa" : ".fa";

    if ($badmerges) {
        my $newbadmerges = "$outputdir\_$i/bad_merges.txt";
        my $transfer_bad_merges_command = 
"transfer_bad_merges.pl -i $badmerges -o $newbadmerges -t $outputdir\_$i/genome.genomex.result/$outputdir\_$i\_old.tsv";
        print "$transfer_bad_merges_command\n";
        system $transfer_bad_merges_command;
        $badmerges = $newbadmerges;
    }

    my $command = "runhm.pl -i $filename -c $configdir -o $outputdir\_$j -p $prefix$j -g";
    $command .= " -s $scaffold_prefix" if $scaffold_prefix;
    $command .= " -b $badmerges"       if $badmerges;
    $command .= " -a $outputdir\_$i/$outputdir\_$i.gff" if $gff;
    print "$command\n";
    system $command;
    my $map_merge_command = "map_merge.py -m $outputdir\_$j -p $prefix$j";
    print "$map_merge_command\n";
    system $map_merge_command;
    my $transfer_merge_command = "transfer_merge.py -d $outputdir\_$i/$outputdir\_$i\_orig.tsv -n $outputdir\_$j/genome.genomex.result/$outputdir\_$j\_old.tsv -o $outputdir\_$j/$outputdir\_$j\_orig.tsv";
    print "$transfer_merge_command\n";
    system $transfer_merge_command;
    my $transfer_command = "transfer_features.py -m $outputdir\_$j/$outputdir\_$j\_orig.tsv -o $outputdir\_$j/$outputdir\_$j";
    $transfer_command .= " -g $gff" if $gff;
    print "$transfer_command\n";
    system $transfer_command;

    $i++;

    last if $rounds > 0 and $rounds < $i;
}

