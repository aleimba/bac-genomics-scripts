#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output

my $usage = "\n".
   "\t####################################################################\n".
   "\t# $0 infile infileformat outfileformat        #\n". #$0 = program name
   "\t#                                                                  #\n".
   "\t# Converts a (multi-)sequence file of a specific format to a       #\n".
   "\t# differently formatted output file, with the help of bioperl.     #\n".
   "\t# www.bioperl.org                                                  #\n".
   "\t# formats e.g.: embl, exp, fasta, genbank                          #\n".
   "\t#                                                                  #\n".
   "\t# Adjust unix loop to run seq_format_converter.pl with several     #\n".
   "\t# files:                                                           #\n".
   "\t# for i in *.gbk; do perl seq_format_converter.pl \$i genbank embl; #\n".
   "\t# done                                                             #\n".
   "\t#                                                                  #\n".
   "\t# version 0.1                                           A Leimbach #\n".
   "\t# 10.11.2011                                 aleimba[at]gmx[dot]de #\n".
   "\t####################################################################\n\n";

# shift arguments from @ARGV
my $infile = shift or die $usage;
my $infileformat = shift or die $usage;
my $outfileformat = shift or die $usage;
# print usage if -h|--h|--help is given as argument
if ($infile =~ m/-h/) {
    die $usage;
}

# SeqIO objects for input and output
my $seq_in = Bio::SeqIO->new(-file => "<$infile", -format => $infileformat); # a Bio::SeqIO object
my $outfile = $infile;
$outfile =~ s/\..*$//; # delete extension from infile name
#print "infile: $infile\toutfile: $outfile\t\$&: $&\n";
my $seq_out = Bio::SeqIO->new(-file => ">$outfile.$outfileformat", -format => $outfileformat); # a Bio::SeqIO object

# write sequence to different format
while (my $seqobj = $seq_in->next_seq) { # a Bio::Seq object
    $seq_out->write_seq($seqobj);
}

print "\n\tCreated new file $outfile.$outfileformat!\n\n";
exit;
