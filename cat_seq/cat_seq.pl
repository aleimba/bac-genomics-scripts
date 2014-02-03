#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqUtils; # bioperl module with additional methods (including features) for Bio::Seq objects

my $usage = "\n".
   "\t#################################################################\n".
   "\t# $0 multi-seq_file [outfile-format]                    #\n". #$0 = program name
   "\t#                                                               #\n".
   "\t# The script merges RichSeq sequences (embl or genbank, but     #\n".
   "\t# also fasta) in a multi-sequence file to one artificial        #\n".
   "\t# sequence. The first sequence in the file is used as a         #\n".
   "\t# foundation to add the subsequent sequences (along with        #\n".
   "\t# features and annotations). Optionally, a different output     #\n".
   "\t# file format can be specified (fasta/embl/genbank).            #\n".
   "\t# The script uses bioperl (www.bioperl.org).                    #\n".
   "\t#                                                               #\n".
   "\t# Adjust unix loop to run the script with all multi-seq files   #\n".
   "\t# in the current working directory, e.g.:                       #\n".
   "\t# for i in *.embl; do cat_seq.pl \$i genbank; done               #\n".
   "\t#                                                               #\n".
   "\t# version 0.1                                        A Leimbach #\n".
   "\t# 08.02.2013                              aleimba[at]gmx[dot]de #\n".
   "\t#################################################################\n\n";

### Shift arguments from @ARGV or give usage
my $multi_seq = shift or die $usage;
my $format = shift;
if ($multi_seq =~/-h/) {
    die $usage;
}


### Bio::SeqIO/Seq objects to concat the seqs
print "\nConcatenating multi-sequence file \"$multi_seq\" to an artificial sequence file ...\n";
my $seqin = Bio::SeqIO->new(-file => "<$multi_seq"); # Bio::SeqIO object; no '-format' given, leave it to bioperl guessing
my @seqs; # store Bio::Seq objects for each seq in the multi-seq file
while (my $seq = $seqin->next_seq) { # Bio::Seq object
    push(@seqs, $seq);
}
Bio::SeqUtils->cat(@seqs);
my $cat_seq = shift @seqs; # the first sequence in the array ($seqs[0]) was modified!


### Write the artificial/concatenated sequence (with its features) to output Bio::SeqIO object
my $seqout; # Bio::SeqIO object
if ($format) { # true if defined
    $multi_seq =~ s/^(.+)\.\w+$/$1_artificial\.$format/;
    $seqout = Bio::SeqIO->new(-file => ">$multi_seq", -format => "$format");
} else {
    $multi_seq =~ s/^(.+)(\.\w+)$/$1_artificial$2/;
    $seqout = Bio::SeqIO->new(-file => ">$multi_seq");
}
$seqout->write_seq($cat_seq);
print "Created new file \"$multi_seq\"!\n\n";

exit;
