#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\n".
    "\t############################################################\n".
    "\t# $0 table_of_coords                         #\n".
    "\t#                                                          #\n".
    "\t# Runs 'trunc_seq.pl' with a tab-delimited input file for  #\n".
    "\t# several sequence files and coordinates. Each sequence    #\n".
    "\t# file and its coordinates have to be in a seperate line   #\n".
    "\t# in the input file. All indicated files have to be        #\n".
    "\t# present in the working directory.                        #\n".
    "\t#                                                          #\n".
    "\t# table_of_coords: seq-file\\tstart\\tstop\\t[outfile-format] #\n".
    "\t# The outfile-format argument is optional, see             #\n".
    "\t# trunc_seq.pl!                                            #\n".
    "\t#                                                          #\n".
    "\t# version 0.1                                   A Leimbach #\n".
    "\t# 08.02.2013                         aleimba[at]gmx[dot]de #\n".
    "\t############################################################\n\n";


### Shift arguments from @ARGV or give usage
my $tbl = shift or die $usage;
if($tbl =~ /-h/) {
    die $usage;
} elsif (!-e $tbl) {
    die "\tfatal error:\t\'$tbl\' doesn't exist!\n\n";
}


### Parse tab-delimited coord-file and run 'trunc_seq_features.pl' with all given files
open(TBL, "<$tbl") or die "Failed to open file \'$tbl\': $!\n";
while (<TBL>) {
    chomp;
    if (/^(.+\.\w+)\t(\d+)\t(\d+)\t(\w+)/) {
	system ("perl trunc_seq_v0.1.pl $1 $2 $3 $4"); # $1: seq-file, $2: start, $3: stop, $4: format
    } elsif (/^(.+\.\w+)\t(\d+)\t(\d+)/) { # without outfile-format
	system ("perl trunc_seq_v0.1.pl $1 $2 $3");
    }
}
close TBL;

exit;
