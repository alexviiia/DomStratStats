# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

package Qvalue;
our $VERSION = 1.01;
use strict;

# my implementations of q-value procedures, most interesting of course is incomplete p-value version (that assumes pi0==1), but I reimplemented other things just for fun

sub p2q {
    # general p2q mapping from p-value list or hash histogram, works with little change for complete and incomplete cases... 
    my ($p2c, $m, $pi0) = defaults(@_);
    # sort p-values from best to worst
    my @ps = sort {$a <=> $b} keys %$p2c;
    # construct p2q map and other cumulatives of interest, no assumptions (in domains case we fill in some unobserved values automatically, here we don't bother)
    my %p2q; # the final map we want!
    my %p2n; # cumulative, also want but only for analysis.
    my %q2n; # a cumulative for actual, final q-values.
    # add either endpoint (p=0 and/or p=1) if it wasn't observed, useful for interpolations!
    unless ($p2c->{0}) {
	$p2n{0} = 0;
	$p2q{0} = 0;
	$q2n{0} = 0;
    }
    unless ($p2c->{1}) { # if it was observed, we have to assume the final count will be correct because the data HAS TO BE complete!  Our assumptions will break badly if p=1 was observed but the data was incomplete in some other sense.
	$p2n{1} = $m; # observe all data at that threshold!
	$p2q{1} = $pi0; # and maximum FDR is pi0!
	$q2n{$pi0} = $m; # and this is cumulative for max q-value
    }
    # compute cumulative of observed p-values, and estimate FDR assuming uniform (this is a preliminary q-value)
    my $n = 0;
    foreach my $p (@ps) {
	$n += $p2c->{$p}; # cumulative in counts
	$p2n{$p} = $n; # regular cumulative, in terms of p
	$p2q{$p} = $n ? $p*$m*$pi0/$n : 0; # this is the preliminary q-value (like a bonferroni correction but lessened by $n and $pi0)
    }
    # to actually get the correct q-values, we need to make them the minimum they can be for any given set (i.e. make monotonic)
    # we can easily get this minimum per threshold by navigating the list backwards
    my $qMin = $pi0; # a good worst case to start with
    foreach my $p (reverse @ps) {
	my $q = $p2q{$p}; # get non-monotonic q
	if ($qMin > $q) {
	    $qMin = $q; # replace current minimum with this q-value if we need to (true most of the time)
	    $q2n{$q} = $p2n{$p}; # since map didn't change, copy over cumulative from p to q (1-to-1 here)
	}
	elsif ($qMin < $q) {
	    $p2q{$p} = $qMin; # here q increased instead of decreased (rare but allowed non-monotonic case we're here to correct), overwrite with previous minimum
	    # note $qMin has already been mapped to an n before, and that previous time the cumulative has had to be larger than for the current p (because we're now going backwards).  Therefore doing nothing is appropriate.
	}
    }
    # done, return data... only return map if scalar context!
    return wantarray ? (\%p2q, \%p2n, \%q2n, $pi0) : \%p2q;
}

sub defaults {
    my ($p2c, $m, $pi0) = @_;
    # convert array to hash if needed
    $p2c = arr2hist($p2c) if ref($p2c) eq 'ARRAY';
    # if the number of tests $m is provided, assume the list is incomplete, so pi0 must be set pessimistically, unless it was set outside
    if ($m) { $pi0 = 1 unless defined $pi0; }
    else {
	# if number of tests is missing, must assume hist is complete and get from there!
	foreach my $c (values %$p2c) {
	    $m += $c;
	}
	# get pi0 from data, because it's complete (unless it was set outside)!
	# Right now this is stupid method, might introduce more options later.
	$pi0 = getPi0($p2c, $m) unless defined $pi0;
    }
    return ($p2c, $m, $pi0); # done, return processed data with no missing values anymore
}

sub getPi0 {
    my ($p2c, $m, $lambda) = @_;
    # this assumes p-value histogram is complete.  computes retarded estimate of pi_0 using only one lambda (default 0.5)
    $lambda = 0.5 unless defined $lambda;
    my $c = 0; # count p-values bigger than lambda
    while (my ($p, $cp) = each %$p2c) { # horrible notation, meh...
	$c += $cp if $p > $lambda; # here we might increase multiple counts at once depending on what $cp is!
    }
    my $pi0 = $c/$m/(1-$lambda); # this is the estimate
    $pi0 = 1 if $pi0 > 1; # for bizarre reasons (highly non-uniform histogram biased for large p-values, or just random chance) this may be larger than 1, obviously wrong, return biggest possible value instead
    return $pi0;
}

sub getPi0s {
    # this assumes p-value list is complete.
    # this plural version tries to be a bit more efficient at the plurality... doesn't take one lambda value, but analyzes all bins of a fixed resolution of 100x in [0,1] (which is what they do at Storey and Tibshirani 2003).
    # However, whole list is returned.  It'd be nicer to just return the one "best" estimate from a smoother near p=0, but I don't need that for my purposes.
    my ($ps) = @_;
    my $m = @$ps; # number of tests
    # build histogram with 100x resolution
    my @hist;
    foreach my $p (@$ps) {
	$hist[int $p*100]++; # floor of percent (like threshold), put in histogram
    }
    # now the histogram can be turned into a cumulative on the fly and into the sampled pi0 estimates
    my @pi0s; # list of estimates we produce, parallel to percent thresholds (same as index)
    # navigate from the top down, that's how cumulative works
    my $c = 0; # the cumulative count
    my $ci; # copy down hist counts, but things could be undefined sometimes
    for (my $i=@hist-1; $i >= 0; $i--) {
	$ci = $hist[$i]; # copy down value
	$c += $ci if $ci; # ignore undefineds
	if ($i < 100) { # can't do for i=100 (div by zero).
	    my $pi0 = $c*100/$m/(100-$i); # this is the estimate of pi0 at this threshold, but will be zero if there were no observations.  
	    $pi0 = 1 if $pi0 > 1; # in some cases this is actually larger than 1, obviously a bad estimate, return biggest possible value instead
	    $pi0s[$i] = $pi0; # now store!
	}
    }
    return \@pi0s; # return the estimates across all data!
}

sub arr2hist {
    # converts a list of p-values into a histogram hash, which is most convenient for very large datasets, and all my code is based on this setup
    my ($ps) = @_; # input array reference
    my %p2c; # the hash we want
    foreach my $p (@$ps) { $p2c{$p}++; }
    return \%p2c;
}

1;
