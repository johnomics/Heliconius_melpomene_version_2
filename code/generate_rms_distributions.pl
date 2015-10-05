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

$args{genetics}    = "";
$args{threads}     = 1;
$args{simulations} = 100000;
$args{outfile}     = "";

my $options_okay = GetOptions(
    'genetics=s'    => \$args{genetics},
    'threads=i'     => \$args{threads},
    'simulations=i' => \$args{simulations},
    'outfile=s'     => \$args{outfile},
);

croak
"No genetics file! Please specify -g to define parents, poor quality individuals and marker types\n"
  if ( $args{genetics} eq "" );

my %genetics;

open my $geneticsfile, "<", $args{genetics}
  or croak "Can't open marker type file $args{genetics}: $OS_ERROR\n";

my $infoline;
while ( $infoline = <$geneticsfile> ) {
    chomp $infoline;
    last if ( $infoline =~ /Type/ );

    if ( $infoline =~ /^Ignore/ ) {
        my ( $ignore, $ind ) = split /\t/, $infoline;
        $genetics{ignore}{$ind} = 0;
    }
    elsif ( $infoline =~ /^Parents/ ) {
        my ( $header, $parents ) = split /\t/, $infoline;
        my @parents = split /,/, $parents;
        $genetics{parents} = \@parents;
    }
    elsif ( $infoline =~ /^Female/ or $infoline =~ /^Male/ ) {
        my ( $sex, $samples ) = split /\t/, $infoline;
        my @samples = split /,/, $samples;
        $genetics{sex}{$sex} = \@samples;
        map { $genetics{samplesex}{$_} = $sex } @samples;
    }
}

my $numdist = 0;
while ( my $marker_type = <$geneticsfile> ) {
    chomp $marker_type;
    my ( $parents, $males, $females, $type ) = split /\t/, $marker_type;
    my $f2s = "$males:$females";
    $genetics{distributions}{$f2s}++;
    $numdist++;
}

close $geneticsfile;

my %data;

my $part_pm = new Parallel::ForkManager( $args{threads} );
$part_pm->set_max_procs( $args{threads} );

$part_pm->run_on_finish(
    sub {
        my $pid       = shift;
        my $exit      = shift;
        my $part      = shift;
        my $childdata = pop;
        merge( $childdata, \%data );
        print STDERR "$part:Merged\n";
    }
);

sub merge {
    my ( $data, $all ) = @_;
    my $f2    = ( keys %{$data} )[0];
    my @vals  = sort { $b <=> $a } keys %{ $data->{$f2} };
    my $inc   = $args{simulations} / 100;
    my $v     = 0;
    my $psims = 0;
    my %pvals;
    while ( $psims < $args{simulations} ) {
        my $p = sprintf "%.2f", $psims / $args{simulations};
        $pvals{$p} = $vals[$v];
        $data->{$f2}{ $vals[$v] } -= $inc;
        while ( $data->{$f2}{ $vals[$v] } < 0 ) {
            $data->{$f2}{ $vals[ $v + 1 ] } += $data->{$f2}{ $vals[$v] };
            $v++;
        }
        $psims += $inc;
    }
    $all->{$f2} = \%pvals;
}

my @dists = sort keys %{ $genetics{distributions} };

foreach my $i ( 1 .. @dists ) {
    $part_pm->start($i) and next;
    my $f2 = $dists[ $i - 1 ];
    printf "%2d %s\n", $i, $f2;
    my ( $male_exp, $female_exp ) = split /:/, $f2;
    my $rms = get_rms_distribution(
        $male_exp, $female_exp,
        scalar @{ $genetics{sex}{Male} },
        scalar @{ $genetics{sex}{Female} }
    );
    $part_pm->finish( 0, { $f2 => $rms } );
}

$part_pm->wait_all_children;

my $outfilename;
if ($args{outfile} ne "") {
    $outfilename = $args{outfile};
}
else {
    $outfilename = $args{genetics};
    $outfilename = $1 if ( $outfilename =~ /(.+)\.txt/ );
    $outfilename .= ".rms_pval.txt";
}

open my $outfile, ">", $outfilename
  or croak "Can't open $outfilename: $OS_ERROR\n";

print $outfile "P value";
map { print $outfile "\t$_" } sort keys %data;
print $outfile "\n";

my $pval = 0.00;
while ( $pval < 1.00 ) {
    $pval = sprintf "%.2f", $pval;
    print $outfile "$pval";
    foreach my $f2 ( sort keys %data ) {
        print $outfile "\t$data{$f2}{$pval}";
    }
    print $outfile "\n";
    $pval += 0.01;
}
close $outfile;

sub get_rms_distribution {
    my ( $male_exp, $female_exp, $males, $females ) = @_;

    my %exp;
    my %male_p;
    my %female_p;
    get_expected_classes( \%exp, 'M', $males,   $male_exp,   \%male_p );
    get_expected_classes( \%exp, 'F', $females, $female_exp, \%female_p );

    # Add small probability of errors
    my $error_classes = 8 - keys %exp;
    for my $c ( 'MA', 'MB', 'MH', 'M.', 'FA', 'FB', 'FH', 'F.' ) {
        $exp{$c} = defined $exp{$c} ? $exp{$c} * 0.9 : 0.1 / $error_classes;
    }

    my $trials = $args{simulations};
    my %rms;
    for my $i ( 1 .. $trials ) {
        printf "%8d /%8d\t%s:%s\n", $i, $args{simulations}, $male_exp,
          $female_exp
          if ( $i % ( $args{simulations} / 10 ) == 0 );
        my %model;
        map {
            my $class = 'F' . get_random_class( \%female_p );
            $model{$class}++
        } ( 1 .. $females );
        map { my $class = 'M' . get_random_class( \%male_p ); $model{$class}++ }
          ( 1 .. $males );

        my $rms_model = sprintf "%.2f", get_rms( \%model, \%exp );
        $rms{$rms_model}++;
    }
    return \%rms;
}

sub get_rms {
    my ( $obs, $exp ) = @_;

    my $rms;
    for my $c ( 'MA', 'MB', 'MH', 'M.', 'FA', 'FB', 'FH', 'F.' ) {
        my $obsv = $obs->{$c} // 0;
        my $expv = $exp->{$c} // 0;
        $rms += ( $obsv - $expv )**2;
    }
    $rms = sqrt( $rms / 8 );

    return $rms;
}

sub get_random_class {
    my ($p) = @_;
    my $rand = rand();
    for my $gt ( keys %{$p} ) {
        return $gt if ( $p->{$gt}{min} < $rand && $rand <= $p->{$gt}{max} );
    }
}

sub get_expected_classes {
    my ( $exp_ref, $sex, $samples, $valid, $p ) = @_;
    my @exp = split /,/, $valid;
    my $shares = 0;
    for my $e (@exp) {
        if ( $e =~ /^(\d)?([ABH.])$/ ) {
            my $share = $1 // 1;
            $exp_ref->{"$sex$2"} = $share;
            $shares += $share;
        }
    }
    my $cum_prob = 0;
    for my $e (@exp) {
        $e = $2 if ( $e =~ /^(\d)([ABH.])$/ );
        my $prob = $exp_ref->{"$sex$e"} / $shares;
        $exp_ref->{"$sex$e"} = $prob * $samples;
        $p->{"$e"}{min} = $cum_prob;
        $cum_prob += $prob;
        $p->{"$e"}{max} = $cum_prob;
    }
}
