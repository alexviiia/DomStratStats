# Copyright (C) 2014-2019 Alejandro Ochoa, Singh Research Group, Princeton University
# This file is part of DomStratStats.
# DomStratStats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DomStratStats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with DomStratStats.  If not, see <http://www.gnu.org/licenses/>.

# core Perl modules
use FindBin ();
# local modules
use lib '.';
use DomStratStats;
use strict;

# this script accepts as input a hmmer3 hmmscan output file (ideally with overlaps removed), and outputs a table file with columns added that correspond to the q-values and local FDRs of these predictions, and additionally filters files given a threshold
my ($fiSeq, $fi, $fo, $qCut, $lCut, $flCut) = @ARGV;

unless ($fo) {
    print "# $FindBin::Script: Compute and add domain q-values, local FDRs, and FDR|lFDR
# " . DomStratStats::version_string() . "
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Usage: perl -w $FindBin::Script <FASTA input> <input table> <output table> \\
         [<q-value> <local FDR> <FDR|lFDR>]

The required inputs are
    <FASTA input>  the FASTA sequence file, used to count all proteins.
    <input table>  the output from 1noOvs.pl (previous script).
    <output table> every input line has three columns added (unless thresolds are set).

The following optional thresholds may be set. Either or all can be set simultaneously, 
    which is supported although it is unusual and hard to interpret, so avoid using this 
    feature.  Use 1 for any threshold you don't wish to set but must list (for example, to 
    only set an FDR|lFDR threshold, set both the q-value and local FDR thresholds to 1).
    <q-value>      stratified domain q-value threshold. Sets FDR per family.
    <local FDR>    stratified domain local FDR threshold. Sets minimum posterior error
                   probability per family.
    <FDR|lFDR>     threshold on combined FDR across domain families implied by a stratified
                   domain local FDR threshold. Sets a stratified domain local FDR threshold 
                   indirectly, by specifying the resulting combined FDR instead.

Input files may be compressed with gzip, and may be specified with or without the .gz 
extension.  Output file will be automatically compressed with gzip.

This is the most basic approach to FDR control.  If you just want domain predictions at an 
FDR equal to Pfam's, set q <= 4e-4 or some other q-value threshold.  Regard lFDRs and 
FDR|lFDR as experimental features worthy of further research, at least for the moment.

Note that the FDR|lFDR is monotonic with the stratified domain lFDR, but not with stratified 
domain q-values.  When lFDRs are accurate (this is a big if!), they provide the optimal way 
of controlling the combined FDR through FDR|lFDR.
";
    exit 0;
}

# constants
my $comp = 'gzip';

# this does all the magic!
DomStratStats::domStratStats($fiSeq, $fi, $fo, $comp, $qCut, $lCut, $flCut);
