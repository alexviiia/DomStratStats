# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

my $VERSION = '1.00';
use lib '.';
use DomStratStats;
use strict;

# this script accepts as input a hmmer3 hmmscan output file (ideally with overlaps removed), and outputs a table file with filters corresponding to the tiered stratified q-value approach, and new columns too (qSeq and qDom|Seq)
my ($fiSeq, $fi, $fo, $qCutSeq, $qCutDom) = @ARGV;

unless ($qCutSeq) {
    print "# $0 $VERSION - Computes and adds q-values for sequences and domains.
# DomStratStats ".(sprintf '%0.2f', $DomStratStats::VERSION).", viiia.org/domStratStats
# Alejandro Ochoa, John Storey, Manuel Llinás, and Mona Singh.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $0 <FASTA input> <input table> <output table> \\ 
         <seq q-value> [<dom q-value>]

The required inputs are
    <FASTA input>  the FASTA sequence file, used to count all proteins.
    <input table>  the output from 1noOvs.pl (previous script).
    <output table> every input line that passes thresholds has two columns added.
    <seq q-value>  stratified sequence q-value threshold. Sets the desired FDR at the
                   sequence level, which can be interpreted as the FDR of calling proteins
                   as having at least one member of a given domain family. This threshold
                   is mandatory because they define conditional domain q-values. If no 
                   domain q-value threshold is set, then this same value is applied to 
                   conditional domain q-values, in which case the domain FDR 
                   (unconditionally) is approximately twice the q-value threshold passed.

You can optionally specify a domain q-value threshold that differs from the sequence
q-value threshold.
    <dom q-value>  threshold on the stratified domain q-value computed conditionally on the 
                   sequence threshold. The expected domain FDR (unconditionally) is 
                   approximately the sum of the sequence and domain q-value thresholds.

Input files may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

See the online manual for more info.
";
    exit 0;
}

# constants
my $comp = 'gzip';

# this sub does all the magic
DomStratStats::tieredStratStats($fiSeq, $fi, $fo, $comp, $qCutSeq, $qCutDom);
