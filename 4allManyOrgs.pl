# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.01';
use lib '.';
use FileGz;
use DomStratStats;
use Hmmer3ScanTab;
use ParsePfam;
use strict;

# some constants
my $qCutSeq = 1e-4; # as in DomStratStats paper
my $comp = 'gzip';
# params, hardcoded for now for convenience
our ($fCut, $lCut) = qw(.5 40); # new permissive overlap params

# this script runs the whole method, from input sequences to GA, DSS and TSS predictions, on many organisms!
# the idea is to 
my ($hmmscan, $pfamA, $fiPfamADat, @orgs) = @ARGV;

unless (@orgs) {
    print "# $0 $VERSION - Get final domain predictions from multiple sequence files
# DomStratStats   ".(sprintf '%0.2f', $DomStratStats::VERSION)." - https://github.com/alexviiia/DomStratStats
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <hmmscan> <Pfam-A.hmm> <Pfam-A.hmm.dat> <orgs>...

The required inputs are
    <hmmscan>        the path to the HMMER3 hmmscan executable.
    <Pfam-A.hmm>     the path to your HMM library of choice (in HMMER3 format).
    <Pfam-A.hmm.dat> Pfam-specific annotations (GA thresholds, clans, nesting).
    <orgs>...        list of prefixes to identify inputs and outputs (see below).

This complex script runs the entire pipeline on several organisms, producing three kinds of 
    outputs that one may wish to compare, namely the Standard Pfam (which use the GA 
    thresholds and other processing), the Domain Stratified Statistics (without thresholds),
    and the Tiered Stratified q-values (with a hardcoded threshold of $qCutSeq per tier).

For each prefix ORG provided, the input is assumed to be ORG.fa or ORG.fa.gz (either will be
    read correctly), and the outputs will be (all compressed with gzip by default, GZ
    extension ommited below)

    ORG.txt       Raw HMMER3 output with p-values and using permissive thresholds
    ORG.noOvs.txt HMMER3 output with overlaps removed (in the permissive sense)
    ORG.dss.txt   Domain Stratified Statistics output, adds 3 columns
    ORG.tsq.txt   Tiered Stratified Q-values output, filtered and adds 2 columns
    ORG.ga.txt    Standard Pfam output, with 'GA' thresholds and intra-clan overlaps removed

Note that, while most steps are run independently per organism, all stratified statistics
    are estimated by pooling the p-value distributions across organisms (always stratified
    by domain family), which are usually better than single-organism estimates.  Moreover,
    the outputs remain separated by organism, which may be more convenient than simply
    merging the initial proteomes into a single sequence file (which produces a single
    output for all organisms).

See the online manual for more info.
"; # there should be at least one organism
    exit 0;
}

# for nesting net and GA thresholds
# in this script, the ile is mandatory
ParsePfam::dat($fiPfamADat); # populates $ParsePfam::nestingNet and $ParsePfam::acc2ds2ga, among other things

# some of the paths we'll need, exactly what is an input or an output is ill-defined (some are both, sorry)
my @fisSeq = map { $_.'.fa' } @orgs; # get clean genome FASTA file 
my @fis = map { $_.'.txt' } @orgs; # rawest data
my @fos = map { $_.'.noOvs.txt' } @orgs; # overlaps removed
my @fosQ = map { $_.'.dss.txt' } @orgs; # dom strat stats
my @fosQ2 = map { $_.'.tsq.txt' } @orgs; # tiered strat stats
my @fosG = map { $_.'.ga.txt' } @orgs; # Standard Pfam

print "Running hmmscan...\n";
for (my $i=0; $i<@fis; $i++) {
    my $fi = $fis[$i];
    my $org = $orgs[$i];
    if ( FileGz::realFile($fi) ) { print "  Skipping $org!\n"; }
    else {
	print "  Doing $org...\n";
	Hmmer3ScanTab::runHmmscan($hmmscan, $pfamA, $fisSeq[$i], $fi);
	# compress raw table when done
	my $ex = system 'gzip', $fi;
	die "Error: 'gzip $fi' returned $ex!\n" if $ex;
    }
}

# remove overlaps!!!
print "Removing overlaps...\n";
DomStratStats::removeOverlaps(\@fis, \@fos, $comp, $ParsePfam::nestingNet, $fCut, $lCut);
# this makes the q-value versions, without setting thresholds
DomStratStats::domStratStats(\@fisSeq, \@fos, \@fosQ, $comp);
# ditto tiered q-values
DomStratStats::tieredStratStats(\@fisSeq, \@fos, \@fosQ2, $comp, $qCutSeq);
# similar output but using standard pfam (for comparison)
DomStratStats::stdPfam(\@fos, \@fosG, $comp, $ParsePfam::acc2ds2ga, $ParsePfam::acc2clan);
