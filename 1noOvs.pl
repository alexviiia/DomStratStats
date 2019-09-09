# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.03';
use lib '.';
use Domains;
use DomStratStats;
use ParsePfam;
use strict;

# this script accepts as input a hmmer3 hmmscan output file, and outputs a table file with overlaps removed by p-value ranking
my ($fi, $fo, $fiPfamADat) = @ARGV;

unless ($fo) {
    print "# $0     $VERSION - Removes overlapping domains ranking by p-value
# DomStratStats ".(sprintf '%0.2f', $DomStratStats::VERSION)." - https://github.com/alexviiia/DomStratStats
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <input table> <output table> [<Pfam-A.hmm.dat>]

The required inputs are
    <input table>    the output from hmmscan (previous script).
    <output table>   input with some domains removed. Format is identical to input.

This is an additional optional input, which is specific to Pfam.
    <Pfam-A.hmm.dat> only used to identify domains that may nest.

Input files may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
";
    exit 0;
}

# constants
my $comp = 'gzip';
# params, hardcoded for now for convenience
our ($fCut, $lCut) = qw(.5 40); # new permissive overlap params

# if file was provided, get nesting net
if (FileGz::realFile($fiPfamADat)) {
    ParsePfam::dat($fiPfamADat); # populates $ParsePfam::nestingNet, among other things
}
else {
    print "Warning: file '$fiPfamADat' was not found, no domain nesting network will be used! (is this a mistake?)\n" if $fiPfamADat;
    # proceed quietly if $fiPfamADat was an empty quote or a zero (or some other obviously false value?)
}

# this does all the magic!
DomStratStats::removeOverlaps($fi, $fo, $comp, $ParsePfam::nestingNet, $fCut, $lCut);
