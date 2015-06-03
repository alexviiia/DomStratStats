# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

package QvalueLocal;
our $VERSION = 1.01;

use lib '.';
use Qvalue;
use strict;

# my implementations of q-value procedures, most interesting of course is incomplete p-value version (that assumes pi0==1), but I reimplemented other things just for fun
# this second version adds Grenander density estimate (via PAV algorithm), which changes q-values slightly in addition to providing robust local FDRs!

### updated for treatment of densities (Modified Grenander density estimator)
# should work like original "fdrtool" when full data is provided (except here I only have crummy pi0 estimates), but handle incomplete domain data as well!!!

sub p2q {
    # general p2q mapping from p-value list or hash histogram, works with little change for complete and incomplete cases... 
    my ($p2c, $m, $pi0) = Qvalue::defaults(@_);
    # these two hacky steps ensure p=0 and p=1 will be part of Grenander density analysis.  Second case will be forced into correct value by minimum density requirement
    # in both cases don't overwrite existing counts, zero only forces presence
    $p2c->{0} = 0 unless $p2c->{0};
    $p2c->{1} = 0 unless $p2c->{1};
    # sort p-values from best to worst
    my @ps = sort {$a <=> $b} keys %$p2c;
    # let's construct the main cumulative (the original one), AND the "modified Grenander" version that keeps cumulative growing at least as fast as uniform random component!
    my %p2n; # cumulative
    my %p2n2; # ditto, but this one is corrected with minimum expected growth from uniform (random) component
    my $n = 0; # running cumulative
    # compute cumulative of observed p-values
    foreach my $p (@ps) {
	$n += $p2c->{$p}; # cumulative in counts
	$p2n{$p} = $n; # regular cumulative, in terms of p
	# enforce minimum cumulative requirement for second count
	my $n2 = $n;
	my $nMin = $p*$m*$pi0;
	$n2 = $nMin if $n2 < $nMin;
	if ($pi0 < 1) { # only apply top bound if pi0 isn't one, otherwise this destroys the data!  (and we don't want this for domains)
	    my $nMax = ($p*$pi0 + 1-$pi0)*$m; # this is value of top constraint
	    $n2 = $nMax if $n2 > $nMax;
	}
	$p2n2{$p} = $n2; # corrected cumulative
    }
    # let's also construct the "least concave majorant" LCM version using the PAV (pool adjacent violators) algorithm (actually consists of deleting points from second "corrected" cumulative)
    pav(\%p2n2, \@ps); # edits %p2n2 by reference
    # compute slopes from PAV output, simple case
    my $p2m2 = slopes(\%p2n2, $pi0, $m);
    my $p2m = slopes(\%p2n, $pi0, $m); # also compute original slopes, for plots (shouldn't really be used!)
    # now density should be fully corrected, q-values can be obtained directly and will be monotonic as desired, and local FDRs are also there!
    # only annoyance is we must interpolate, since PAV deleted points
    interpolate(\%p2n2, $p2m2, \@ps); # this interpolates very cleanly and easily
    my %p2dat; # let's put all things together in a single table, maybe that makes more sense...
    my %p2ql; # return only basic map in scalar context, it's useful for quick on-the-fly maps!
    for (my $i=0; $i<@ps; $i++) {
	# all values should be defined now and are ready to be merged into a single table!
	my $p = $ps[$i];
	my $n = $p2n{$p};
	my $n2 = $p2n2{$p};
	my $m1 = $p2m->{$p}; # don't overwrite $m (the number of tests)
	my $m2 = $p2m2->{$p};
	# actual new computations!
	my $q = $p ? $p*$m*$pi0/$n2 : 0; # q-value (really direct FDR, but cumulative should work appropriately now because of Grenander density estimator).  Only stupid special case is for p=0, which might not have been observed (always set q to zero regardless)
	my $l = $m*$pi0/$m2; # local FDR using this data
	$l = 1 if $l > 1; # ensure we don't output non-sensical lFDRs due to numerical oddities!
	$p2dat{$p} = [$q, $l, $n, $n2, $m1, $m2]; # store data we want!  Q and L are most important, list first.  Rest is for potential analysis later on
	$p2ql{$p} = $l; # save only local FDR separately!
    }
    # done, return data
    return wantarray ? (\%p2dat, $pi0, \@ps) : \%p2ql;
}

sub fdrOfLocalThresh {
    # for a given local FDR threshold, the resulting FDR is theoretically the average local FDR of everything declared significant
    # so if we have a sorted list of local FDRs, the FDR at a threshold is the average local FDR observed above that threshold
    # NOTE: this method performs poorly numerically (it's not a good idea to average lots of tiny numbers, each of which has roundoff error and might not be estimated very accurately in the first place).  Please use the other method below, which is based on q-values, it should be much more robust!!!
    my ($l2c) = @_; # input is a local FDR histogram (should be more compact than a list)
    my @ls = sort { $a <=> $b } keys %$l2c; # sort unique local FDRs ascending numerically
    my %l2fdr; # the relationship we want
    # we want the cumulative sum of local FDRs and their counts, to get the average we want
    my $lSum = 0;
    my $cSum = 0;
    foreach my $l (@ls) {
	my $c = $l2c->{$l};
	$lSum += $c*$l; # we want this local FDR summed with this multiplicity (because local FDRs are converted to text in keys, they become discretized so some values will be observed multiple times)
	$cSum += $c; # only add count here
	$l2fdr{$l} = $cSum ? $lSum/$cSum : 0; # this is the FDR at this local FDR threshold!  Note p=0 always gets a non-zero lFDR value, even if c==0 (and hence cSum==0), handle it appropriately.
    }
    return \%l2fdr; # the mapping we want!
}

sub fdrOfLocalThreshFromQ {
    # although averaging local FDR to get the FDR is technically correct (the previous function above), it is numerically fraught.
    # FDR should also be average q-value weighted by preds given the local FDR threshold, so we should be able to compute this FDR more robustly from here
    # unfortunately we need a more memory-hungry structure to compute these things
    my ($l2acc2qc) = @_;
    
    my %l2fdr; # the data we want
    # temp vars
    my %acc2qcCum; # cumulative counts per family and q-value thresholds
    my $cCum = 0; # the total cumulative sum of counts, most convenient to keep it separate instead of relooping through data above, even though it is contained there too
    # sort local FDRs
    my @ls = sort { $a <=> $b } keys %$l2acc2qc; # sort numerically ascending
    foreach my $l (@ls) {
	my $acc2qc = $l2acc2qc->{$l}; # get the data at this local FDR threshold
	# potentially this lFDR threshold affects many accs, let's traverse each of them
	while (my ($acc, $qc) = each %$acc2qc) {
	    my ($q, $c) = @$qc; # extract from array
	    # update each of the data in the cumulative structure
	    my $qcCum = $acc2qcCum{$acc}; # copy down structure, but it may be undefined!
	    $qcCum->[0] = $q; # the first element is the q-value threshold, just replace
	    $qcCum->[1] += $c; # second element is cumulative count for this family, update!
	    $acc2qcCum{$acc} = $qcCum; # overwrite structure, only needed if $qcCum was previously undefined!
	    $cCum += $c; # update global count too
	}
	# cool, now we compute the weighted q-value we want...
	if ($cCum == 0) { $l2fdr{$l} = 0; } # easy case for p=0 interpolations
	else {
	    my $fdr = 0;
	    foreach my $qcCum (values %acc2qcCum) { # accs don't matter in this loop
		my ($q, $c) = @$qcCum; # extract from array
		$fdr += $q*$c; # add this weighted term, unnormalized yet
	    }
	    $l2fdr{$l} = $fdr/$cCum; # normalize and store
	}
    }
    # done, return the mapping we want
    return \%l2fdr;
}


### internal subs, no need to call them outside usually!

sub pav {
    # applies pool-adjacent-violators (PAV) algorithm
    # violators simply get deleted from hash by reference!  Easiest way to deal with this!
    my ($p2n, $ps) = @_;
    # because we'll be deleting points, let's remember last two indexes not deleted.  These are the values before deleting has happened!
    my $i2 = 1;
    my $i3 = 0; 
    # compute cumulative of observed p-values
    for (my $i=2; $i < @$ps; $i++) {
	# let's compute last two slopes and merge points if slopes increased instead of decreasing, as is expected!  (or reverse if we're moving in reverse)
	my $p1 = $ps->[$i];
	my $p2 = $ps->[$i2];
	my $p3 = $ps->[$i3];
	my $n1 = $p2n->{$p1};
	my $n2 = $p2n->{$p2};
	my $n3 = $p2n->{$p3};
	my $m2 = ($n1-$n2)/($p1-$p2); # slope between current and last point
	my $m3 = ($n2-$n3)/($p2-$p3); # slope between last and second to last point
	if ($m3 < $m2) { # this is a violation of monotoniticy (slopes are supposed to decrease in this problem, but this is an increase; unless we're doing reverse)
	    delete $p2n->{$p2}; # the middle point gets removed, this will make things work out
#	    print join("\t", 'DEL+', $i2, $p2, $m3, $m2)."\n"; ###DEBUG
	    # in this case i3 stays the same
	    # unfortunately this triggers checking backwards for more monotonicity violations, all must be fixed before continuing forward!
	    my $isBad = 1; # boolean that tells me if there's violators still!
	    while ($isBad && $i3 > 0) { # things must still be bad AND we must still be in range (or there's nothing to test anymore)
		my $i4 = $i3-1; # search for next lowest point that hasn't been deleted, start with this one
		my $p4 = $ps->[$i4];
		my $n4 = $p2n->{$p4};
		while (!defined $n4) { # if this was deleted, keep looking back!
		    $i4--; # step back again.
		    $p4 = $ps->[$i4]; # update values again, retest!
		    $n4 = $p2n->{$p4};
		}
		# we're guaranteed i4>=0 because p[i=0] is always defined (never deleted), so this should always behave
		# cool, so we found last point that wasn't deleted!
		# compute these slopes and see if they're fine or not (overwrite old slopes, not needed anymore, but note the points are different because i2 was deleted!)
		$m2 = ($n1-$n3)/($p1-$p3); # slope between current and last point
		$m3 = ($n3-$n4)/($p3-$p4); # slope between last and second to last point
		if ($m3 < $m2) { # another violation was found, keep going backward!
		    delete $p2n->{$p3}; # now p3 is the point in the middle we need to delete!
#		    print join("\t", 'DEL-', $i3, $p3, $m3, $m2)."\n"; ###DEBUG
		    $i3 = $i4; # and i3 needs to be pushed back to last non-deleted point (which was i4 here)
		    $p3 = $ps->[$i3]; # recompute these for loop to continue!
		    $n3 = $p2n->{$p3};
		    # and isBad=1 still.  This will keep backward loop going!
		}
		else {
#		    print join("\t", 'STOP-', $i3, $p3, $m3, $m2)."\n"; ###DEBUG
		    $isBad = 0; # can stop loop, things are good!  Note i3 stays the same here!
		}
	    }
	}
	else { $i3 = $i2; } # in regular case i3 advances to old i2
	$i2 = $i; # in all cases i2 becomes current i.  $i will be incremented automatically
    }
}

sub slopes {
    # computes and stores slopes from PAV output (reduced cumulative) into parallel structure
    # easy to do here because no data is missing
    my ($p2n, $pi0, $m) = @_;
    my @ps = sort {$b <=> $a} keys %$p2n; # get p-values to analyze, smallest first (order doesn't really matter here)
    my %p2m; # the data we want
    # for worst p-value (at zero index) slope is pi0*m, doesn't come directly from data but is assumed by model
    my $i=0;
    my $p = $ps[$i];
    my $n = $p2n->{$p};
    my $m2 = $pi0*$m;
    $p2m{$p} = $m2;
    for ($i=1; $i<@ps; $i++) { # for every subsequent point (skip i=0 cause it was one already)
	my $pi = $ps[$i];
	my $ni = $p2n->{$pi};
	$m2 = ($n-$ni)/($p-$pi); # slope between this and previous point
	$p2m{$pi} = $m2; # store slope in structure.  Note slope of a pair goes to the left point (with smaller p-value)
	$p = $pi; $n = $ni; # replace prev point with this point, for next round
    }
    return \%p2m; # done, here's slopes structure!!!
}

sub interpolate {
    # get values for n (cumulative) and m (density, slope) for every p in @$ps (complete p-value list), given incomplete p2n,p2m structures
    # n's get linearly interpolated, m's are piecewise constant (copied to the right)
    # assume @$ps are sorted such that smallest p-values come first
    # input structures %$p2n, %$p2m are updated directly by reference, nothing is returned!
    my ($p2n, $p2m, $ps) = @_;
    # get data from very first point, will store previous points too.  It should be defined in structure too because the first point is never eliminated
    my $i = 0;
    my $p = $ps->[$i];
    my $n = $p2n->{$p};
    my $m = $p2m->{$p};
    for ($i=1; $i<@$ps; $i++) { # navigate rest of p-values
	my $pi = $ps->[$i]; # always defined
	my $ni = $p2n->{$pi}; # may be missing
	my $mi = $p2m->{$pi}; # ditto
	if (defined $ni) {
	    # no interpolation is needed, but copy values for next round
	    $p = $pi; $n = $ni; $m = $mi;
	}
	else {
	    # interpolate!  And old values stay the same (because they might continue to be needed, interpolation equations will work with same values)
	    $p2n->{$pi} = $n + ($pi-$p)*$m; # this successfully reverses slope and last point (and this point's pi) into this point's ni (linear interpolation)
	    $p2m->{$pi} = $m; # the slope stays constant (same as last point)
	}
    }
    # done!  Edited references so nothing to return!
}

1;
