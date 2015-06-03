# Copyright (C) 2014 Alejandro Ochoa, Singh Research Group, Princeton University
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

package ParseFasta;
our $VERSION = 1.00;
use Compress::Zlib; # core module!
use lib '.';
use FileGz;
use strict;

# lightweight FASTA parser in a simple structure
# some variants emphasize low memory footprint
# there are no output writers because I've never needed them (FASTA is super easy to output if no max column limits are needed)

sub parseIntoHash {
    # reads regular fasta into a hash
    my ($fi, $doCompress) = @_;
    
    # a few variables
    my %id2seq;
    my $currId = '';
    my $currSeq = '';
    
    # read in file
    my $fhi = FileGz::getInFh($fi);
    while (<$fhi>) {
	chomp;
	if (s/^>//) {
	    if ($currSeq) { # so that nothing happens for the first protein name
		$currSeq = compress $currSeq if $doCompress; # compress if needed
		$id2seq{$currId} = $currSeq; # store last id/sequence pair
		$currSeq = ''; # now erase it to make way for next sequence
	    }
	    $currId = $_; # update protein name
	}
	else {
	    s/\s+//g; # remove all whitespace
	    $currSeq .= $_; # add to previous buffer
	}
    }
    close $fhi;
    # process/store very last id/sequence pair
    $currSeq = compress $currSeq if $doCompress; # compress if needed
    $id2seq{$currId} = $currSeq;
    
    # done, return
    return \%id2seq;
}

sub parseIds {
    # reads regular fasta, returns a "set" of IDs only
    # idea is to only get what we want because file is huge (low-mem setting)
    # works with hits files too!!!!!
    my ($fi) = @_;
    my %ids; # the set we want
    # read in file
    my $fhi = FileGz::getInFh($fi);
    while (<$fhi>) {
	chomp;
	$ids{$_} = 1 if s/^>//;
    }
    close $fhi;
    return \%ids; # done, return
}

sub parseInit {
    # this returns filehandle and the first ID (need to give that to next sub)
    # if no IDs were found, will return a null ID
    my $fi = shift;
    my $fhi = FileGz::getInFh($fi);
    my $firstId;
    # this will read until the first ID is found (ignore everything before that)
    while (<$fhi>) {
	chomp;
	if (s/^>//) {
	    $firstId = $_; # update protein name
	    last;
	}
    }
    # awesome, return
    return ($fhi, $firstId);
}

sub parseNext {
    # this takes the last ID read and returns the next ID and sequence from fasta file
    # if no IDs were found, will return a null ID
    my ($fhi, $id) = @_;
    my $seq = '';
    my $idNext;
    # read next entry
    while (<$fhi>) {
	chomp;
	if (s/^>//) { $idNext = $_; last; } # this marks the beginning of the next protein, stop loop
	else {
	    s/\s+//g; # remove all whitespace
	    $seq .= $_; # add to previous buffer
	}
    }
    return ($id, $seq, $idNext);
}

sub countSeqs {
    my ($fis) = @_;
    # counts number of IDs across multiple FASTA files potentially
    my $c = 0; # the number we want
    # go through every file in input, handling scalar case especially
    foreach my $fi (ref $fis ? @$fis : $fis) {
	my $fhi = FileGz::getInFh($fi);
	while (<$fhi>) { $c++ if /^\>/; } # found a new protein
	close $fhi;
    }
    return $c;
}

sub countLetters {
    my ($fis) = @_;
    # counts number of letters across multiple FASTA files potentially
    # works for amino acids and bases.  All non-ID and non-whitespace is counted (so don't feed it weird stuff cause no errors will be caught!).
    my $c = 0; # the number we want
    # go through every file in input, handling scalar case especially
    foreach my $fi (ref $fis ? @$fis : $fis) {
	my $fhi = FileGz::getInFh($fi);
	while (<$fhi>) {
	    chomp;
	    unless (/^>/) { # ignore ID lines
		s/\s+//g; # remove all whitespace
		$c += length; # add length of $_ (current line minus whitespace)
	    }
	}
	close $fhi;
    }
    return $c;
}

1;
