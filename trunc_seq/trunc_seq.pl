#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqUtils; # bioperl module with additional methods (including features) for Bio::Seq objects

my $usage = "\n".
   "\t###########################################################\n".
   "\t# $0 seq_file start stop [outfile-format]       #\n". #$0 = program name
   "\t#                                                         #\n".
   "\t# Truncates a sequence file (optionally with features/    #\n".
   "\t# annotations) according to the given coordinates. The    #\n".
   "\t# truncated sequence is written to an outfile in the same #\n".
   "\t# sequence format (embl/genbank/fasta). Alternatively, a  #\n".
   "\t# specific outfile-format can be given.                   #\n".
   "\t# The script uses bioperl (www.bioperl.org).              #\n".
   "\t#                                                         #\n".
   "\t# version 0.1                                  A Leimbach #\n".
   "\t# 08.02.2013                        aleimba[at]gmx[dot]de #\n".
   "\t###########################################################\n\n";


### Shift arguments from @ARGV or give usage
my $seq_file = shift or die $usage;
my $start = shift or die $usage;
my $stop = shift or die $usage;
my $format = shift;
if ($seq_file =~ m/-h/) {
    die $usage;
}


### Bio::SeqIO/Seq object to truncate
print "\nTruncating \"$seq_file\" to coordinates $start..$stop ...\n";
my $seqin = Bio::SeqIO->new(-file => "<$seq_file"); # Bio::SeqIO object; no '-format' given, leave it to bioperl guessing
my $seqobj = $seqin->next_seq; # Bio::Seq object
my $truncseq = Bio::SeqUtils->trunc_with_features($seqobj, $start, $stop);


### Write truncated sequence to output Bio::SeqIO object
my $seqout; # Bio::SeqIO object
$seq_file =~ s/.+\/(.+)$/$1/;
if ($format) { # true if defined
    $seq_file =~ s/^(.+)\.\w+$/$1_$start\_$stop\.$format/;
    $seqout = Bio::SeqIO->new(-file => ">$seq_file", -format => "$format");
} else {
    $seq_file =~ s/^(.+)(\.\w+)$/$1_$start\_$stop$2/;
    $seqout = Bio::SeqIO->new(-file => ">$seq_file");
}
$seqout->write_seq($truncseq);
print "Created new file \"$seq_file\"!\n\n";

exit;
