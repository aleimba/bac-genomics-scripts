#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqUtils; # bioperl module with additional methods (including features) for Bio::Seq objects

my $usage = "\n".
   "\t###########################################################\n".
   "\t# $0 seq_file [outfile-format]                 #\n". #$0 = program name
   "\t#                                                         #\n".
   "\t# Reverse complement the sequence in a sequence file and  #\n".
   "\t# write the result to an outfile in the same sequence     #\n".
   "\t# format. For RichSeq files with features/annotation the  #\n".
   "\t# features will be adapted accordingly. Optionally, a     #\n".
   "\t# different output file format can be specified           #\n".
   "\t# (fasta/embl/genbank).                                   #\n".
   "\t# The script uses bioperl (www.bioperl.org).              #\n".
   "\t#                                                         #\n".
   "\t# Adjust unix loop to run the script with all files in    #\n".
   "\t# the current working directory:                          #\n".
   "\t# for i in *.embl; do revcom_seq.pl \$i genbank; done      #\n".
   "\t#                                                         #\n".
   "\t# version 0.1                                  A Leimbach #\n".
   "\t# 2013-02-08                        aleimba[at]gmx[dot]de #\n".
   "\t###########################################################\n\n";

### Shift argument from @ARGV or give usage
my $seq_file = shift or die $usage;
my $format = shift;
if ($seq_file =~ /-h/) {
    die $usage;
}


### Bio::SeqIO/Seq object to reverse complement
print "\nReverse complementing \"$seq_file\" ...\n";
my $seqin = Bio::SeqIO->new(-file => "<$seq_file"); # Bio::SeqIO object; no '-format' given, leave it to bioperl guessing
my $seq_obj = $seqin->next_seq; # Bio::Seq object
my $revcom = Bio::SeqUtils->revcom_with_features($seq_obj);


### Write reverse complemented sequence (and its features) to output Bio::SeqIO object
my $seqout; # Bio::SeqIO object
if ($format) { # true if defined
    $seq_file =~ s/^(.+)\.\w+$/$1_revcom\.$format/;
    $seqout = Bio::SeqIO->new(-file => ">$seq_file", -format => "$format");
} else {
    $seq_file =~ s/^(.+)(\.\w+)$/$1_revcom$2/;
    $seqout = Bio::SeqIO->new(-file => ">$seq_file");
}
$seqout->write_seq($revcom);
print "Created new file \"$seq_file\"!\n\n";

exit;
