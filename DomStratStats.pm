# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

package DomStratStats;
our $VERSION = 1.03;

use lib '.';
use Hmmer3ScanTab;
use Qvalue;
use QvalueLocal;
use ParseFasta;
use Domains;
use strict;

# core functions that can be reused in different contexts

sub stdPfam {
    # sets the gathering threhsolds and other pfam-specific stuff (overlap removal within clans only)
    # produces outputs in the same style as DomStratStats, for comparison to our method (and practicality otherwise)
    # true to the original procedure, overlaps are treated strictly (not permissively like in the rest of this package)
    my ($fis, $fos, $comp, $acc2ds2ga, $acc2clan) = @_;
    
    # figure out inputs
    my @fis = getArr($fis);
    my @fos = getArr($fos);
    die "Error: input and output arrays are not the same length (".(scalar @fis)." vs ".(scalar @fos)."\n" unless @fis == scalar @fos;
    
    # here each file can be analyzed separately
    for (my $i=0; $i<@fis; $i++) {
	my $fi = $fis[$i];
	my $fo = $fos[$i];
	
	# now read the full files entry-by-entry
	my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
	# prepare output
	my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
	while ($protNext) {
	    # read next protein
	    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
	    # apply gathering thresholds filter
	    $hits = Domains::ga($hits, $acc2ds2ga);
	    # this removes remaining overlaps using clans.  Also sorts hits by start
	    $hits = Domains::filterPfamClans($hits, $acc2clan);
	    # print to output
	    Hmmer3ScanTab::printProt($fho, $hits);
	}
	
	# done with Hits input and output files
	close $fhi;
	close $fho;
    }
}

sub removeOverlaps {
    my ($fis, $fos, $comp, $nestingNet, $fCut, $lCut) = @_;
    
    # figure out inputs
    my @fis = getArr($fis);
    my @fos = getArr($fos);
    die "Error: input and output arrays are not the same length (".(scalar @fis)." vs ".(scalar @fos)."\n" unless @fis == scalar @fos;
    
    # here each file can be analyzed separately
    for (my $i=0; $i<@fis; $i++) {
	my $fi = $fis[$i];
	my $fo = $fos[$i];
	
	# now read the full files entry-by-entry
	my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
	# prepare output
	my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
	while ($protNext) {
	    # read next protein
	    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
	    # remove overlaps, overwritting hits.  Also sorts hits by start
	    $hits = Domains::removeOverlapsByScore($hits, $nestingNet, $fCut, $lCut);
	    # print to output
	    Hmmer3ScanTab::printProt($fho, $hits);
	}
	
	# done with Hits input and output files
	close $fhi;
	close $fho;
    }
}

sub domStratStats {
    # in some cases we want to consider multiple organisms (with separate files) in combination, let's handle that transparently here!
    # sequence files do not need to be consistent with HMMER files (in number of elements), but input and output HMMER3 files should have the same length!
    my ($fisSeq, $fis, $fos, $comp, $qCut, $lCut, $flCut) = @_;
    
    # figure out inputs
    my @fisSeq = getArr($fisSeq);
    my @fis = getArr($fis);
    my @fos = getArr($fos);
    die "Error: input and output arrays are not the same length (".(scalar @fis)." vs ".(scalar @fos)."\n" unless @fis == scalar @fos;
    
    # 1) Count all proteins tested, some of which may not be in HMMER3 output, hence need source FASTA file
    print "Counting proteins in FASTA file(s)...\n";
    my $m = ParseFasta::countSeqs(\@fisSeq);
    
    # 2) HMMER3 input FIRST PASS: only collect p-value distributions, separately per family
    # this code does it very efficiently
    print "Collecting p-value distributions...\n";
    my ($name2p2c, $name2m) = Hmmer3ScanTab::parsePValues(\@fis);
    
    # 3) compute q-values and local FDRs
    print "Computing q-values and local FDRs...\n";
    my %name2p2str; # all we really need is the formatted string to add, and here we set filters too!
    my %l2name2qc; # a robust way to compute FDR_l from local FDR threshold data
    while (my ($name, $p2c) = each %$name2p2c) {
	# process adjustment to the number of tests
	my $mName = $name2m->{$name}; # get adjustment for repeating domains, this may be undefined
	$mName = 0 unless defined $mName; # set to zero if undefined
	$mName += $m; # add the global count, now this is fine
	# actually compute q-values and local FDRs with my special packages
	my $p2q = Qvalue::p2q($p2c, $mName);
	my $p2l = QvalueLocal::p2q($p2c, $mName);
	# now we consolidate these numbers into a single string, AND we filter if requested
	my %p2str; # the map we want
	my @ps = sort {$a <=> $b} keys %$p2q; # should be the same for q-values and local FDR
	foreach my $p (@ps) {
	    # get both stats
	    my $q = $p2q->{$p};
	    my $l = $p2l->{$p};
	    my $c = $p2c->{$p}; # for FDR|lFDR
	    # we have to watch out with overwriting data from different p-values but the same lFDR!
	    my $lStr = sprintf '%9.2g', $l;
	    my $qc = $l2name2qc{$lStr}{$name}; # get previous data, but it will be undefined if there was no previous data
	    if ($qc) { # update data in this case
		$qc->[0] = $q; # I think in this case the larger of the q-values is best (if we used the q-value data directly from QvalueLocal, like Strimmer's method, then we'd have consistency I think, but meh).    Anyway, in this order of navigating p-values, we should automatically overwrite to have the largest q-values anyway, so this works!  Doing this is pesimistic (or conservative, unless it's just right), compared to the alternative.
		# however, the smaller q-value will give a more optimistic FDR_l, maybe this is the correct behavior (in this case do not update q in $qc)
		# I hope in most cases where different p-values give the same lFDR, that they give the same q-value too or at least very very similar ones
		# and note this decision doesn't affect the cumulative
		$qc->[1] += $c; # sum of counts, easy
	    }
	    else { $qc = [$q, $c]; } # save for first time
	    $l2name2qc{$lStr}{$name} = $qc; # overwrite updated data
	    # apply thresholds; if not met, then this p-value won't be in map (and will further save memory?)
	    if ($qCut) { next if $q > $qCut; }
	    if ($lCut) { next if $l > $lCut; }
	    # now format and store!
	    $p2str{$p} = sprintf '%9.2g %9.2g', $q, $l; # matches HMMER3's domain table format here
	}
	# add to larger structure
	$name2p2str{$name} = \%p2str;
    }
    
    # 3.5) finish off FDR_l (for local FDR thresholds).
    my $l2fdr = QvalueLocal::fdrOfLocalThreshFromQ(\%l2name2qc);
    
    # 4) HMMER3 input SECOND PASS: add stats to output, also filter if requested
    print "Generating output file with new columns (q, lFDR, and FDR|lFDR)...\n";
    # here each file can be analyzed separately
    for (my $i=0; $i<@fis; $i++) {
	my $fi = $fis[$i];
	my $fo = $fos[$i];
	# now read the full files entry-by-entry
	my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
#	# matches HMMER3's domain table format
#	my $headerExtra1 = sprintf '%19s', '--- strat stats ---';
#	my $headerExtra2 = sprintf '%9s %9s', 'q-value', 'local FDR';
#	my $headerExtra3 = sprintf '%9s %9s', '-' x 9, '-' x 9;
	# matches HMMER3's domain table format
	my $headerExtra1 = sprintf '%29s', '---- domain strat stats -----';
	my $headerExtra2 = sprintf '%9s %9s %9s', 'q-value', 'local FDR', 'FDR|lFDR';
	my $headerExtra3 = sprintf '%9s %9s %9s', '-' x 9, '-' x 9, '-' x 9;
	$header = Hmmer3ScanTab::addColsHeader($header, $headerExtra1, $headerExtra2, $headerExtra3);
	# prepare output
	my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
	while ($protNext) {
	    # read next protein
	    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
	    # add columns or omit from output!
	    my @hitsOut; # since we may set thresholds, let's keep only passing domains in this list
	    foreach my $h (@$hits) {
		my $str = $name2p2str{$h->{name}}{$h->{E}}; # get new substring (contains new columns only)
		next unless $str; # skip this domain if there's no substring (means it didn't pass thresholds)
		# add fdr_l here!
		# recover local FDR just to query it in $l2fdr
		my $l = substr $str, 10, 9; # I think this will work, to extract correct value regardless of whitespace situation (just as it was stored as key)
		my $fdr_l = $l2fdr->{$l}; # query, if something went wrong this'll be undefined
		die "Error: local FDR '$l' was not found in hash l2fdr.  Contact authors to correct bug!\n" unless defined $fdr_l; # consider this a serious error!!!
		# only now can we set the FDR_l threshold, let's set it if needed
		if ($flCut) { next if $fdr_l > $flCut; } 
		$str .= ' '.sprintf '%9.2g', $fdr_l; # add FDR_l to string with correct formatting
		$h->{origLine} = Hmmer3ScanTab::addCols($h->{origLine}, $str); # this adds the new columns correctly
		push @hitsOut, $h; # now addd to list of passing hits!
	    }
	    # print to output
	    Hmmer3ScanTab::printProt($fho, \@hitsOut);
	}
	
	# done with Hits input and output files
	close $fhi;
	close $fho;
    }
    
    print "Done!\n";
}

sub tieredStratStats {
    my ($fisSeq, $fis, $fos, $comp, $qCutSeq, $qCutDom) = @_;
    # if domain threshold is not specified, use same as for sequences
    $qCutDom = $qCutSeq unless defined $qCutDom;
    
    # constants
    my $pi0 = 1; # always the case for domains!
    my $doTiered = 1; # a named boolean
    
    # figure out inputs
    my @fisSeq = getArr($fisSeq);
    my @fis = getArr($fis);
    my @fos = getArr($fos);
    die "Error: input and output arrays are not the same length (".(scalar @fis)." vs ".(scalar @fos)."\n" unless @fis == scalar @fos;
    
    # 1) Count all proteins tested, some of which may not be in HMMER3 output, hence need source FASTA file
    print "Counting proteins in FASTA file...\n";
    my $m = ParseFasta::countSeqs(\@fisSeq);
    
    # 2) HMMER3 input FIRST PASS: only collect p-value distributions, separately per family
    # this code does it very efficiently
    print "Collecting tiered p-value distributions...\n";
    my ($name2pSeq2c, $name2pSeq2pDom2c) = Hmmer3ScanTab::parsePValues(\@fis, $doTiered);
    
    # 3) compute tiered q-values
    print "Computing tiered q-values...\n";
    my %name2pSeq2q; # these are sequence q-values
    my %name2pDom2q; # these are domain q-values conditional on the sequence q-value being set
    while (my ($name, $pSeq2c) = each %$name2pSeq2c) {
	# compute q-values for sequence p-values with my special package
	my $pSeq2q = Qvalue::p2q($pSeq2c, $m); # pi0=1 implicitly
	# now navigate tiered domain data and set first threshold
	my $pSeq2pDom2c = $name2pSeq2pDom2c->{$name}; # copy down structure
	my %pDom2c; # the conditional domain p-values
	while (my ($pSeq, $pDom2c) = each %$pSeq2pDom2c) {
	    if ($pSeq2q->{$pSeq} <= $qCutSeq) { # if first tier passes, let's analyze these domains!
		while (my ($pDom, $c) = each %$pDom2c) {
		    $pDom2c{$pDom} += $c;
		}
	    }
	}
	# now compute the domain q-values conditional on the sequence q-value threshold
	my $pDom2q = Qvalue::p2q(\%pDom2c, undef, $pi0); # this will compute number of tests automatically from data, but will keep a conservative pi0=1 (not implicitly set in this case!)
	# add to larger structures
	$name2pSeq2q{$name} = $pSeq2q;
	$name2pDom2q{$name} = $pDom2q;
    }
    
    # 4) HMMER3 input SECOND PASS: add stats to output, also filter
    print "Generating filtered output file with new columns (qSeq and qDom|Seq)...\n";
    # here each file can be analyzed separately
    for (my $i=0; $i<@fis; $i++) {
	my $fi = $fis[$i];
	my $fo = $fos[$i];
	# now read the full files entry-by-entry
	my ($fhi, $header, $protNext, $hitNext) = Hmmer3ScanTab::parseInit($fi);
	# matches HMMER3's domain table format
	my $headerExtra1 = sprintf '%19s', 'tiered strat stats';
	my $headerExtra2 = sprintf '%9s %9s', 'qSeq', 'qDom|Seq';
	my $headerExtra3 = sprintf '%9s %9s', '-' x 9, '-' x 9;
	$header = Hmmer3ScanTab::addColsHeader($header, $headerExtra1, $headerExtra2, $headerExtra3);
	# prepare output
	my $fho = Hmmer3ScanTab::printInit($fo, $header, $comp);
	while ($protNext) {
	    # read next protein
	    (my $prot, my $hits, $protNext, $hitNext) = Hmmer3ScanTab::parseProt($fhi, $protNext, $hitNext);
	    # filter and add columns!
	    my @hitsOut; # list of passing domains
	    foreach my $h (@$hits) {
		my $name = $h->{name}; # copy things that may be reused
		my $qSeq = $name2pSeq2q{$name}{$h->{ESeq}}; # sequence q-value
		print "Prot $prot, name $name, pSeq $h->{ESeq} has undefined qSeq!\n" unless defined $qSeq; ### DEBUG
		next if $qSeq > $qCutSeq; # set sequence threshold on this hit
		my $qDom = $name2pDom2q{$name}{$h->{E}}; # domain q-value conditional on sequence threshold
		next if $qDom > $qCutDom; # set domain threshold on this hit
		my $str = sprintf '%9.2g %9.2g', $qSeq, $qDom; # matches HMMER3's domain table format here
		$h->{origLine} = Hmmer3ScanTab::addCols($h->{origLine}, $str); # this adds the new columns correctly
		push @hitsOut, $h; # now addd to list of passing hits!
	    }
	    # print to output
	    Hmmer3ScanTab::printProt($fho, \@hitsOut);
	}
	
	# done with Hits input and output files
	close $fhi;
	close $fho;
    }
    
    print "Done!\n";
}

sub getArr {
    my ($arr) = @_;
    if ('ARRAY' eq ref $arr) { return @$arr; } # dereference array
    elsif (!ref $arr) { return ($arr); } # make array with a single scalar element
    else { die "Error: expected array but instead got ".(ref $arr)." !\n"; } # die if things were unexpected
}

1;
