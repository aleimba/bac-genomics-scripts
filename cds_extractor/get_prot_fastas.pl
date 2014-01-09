#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SeqFeatureI;  # bioperl module to handle features in a sequence

my $usage = "\n".
   "\t###############################################################\n".
   "\t# $0 seq-file locus_tag-file                  #\n". #$0 = program name
   "\t#                                                             #\n".
   "\t# Extracts all protein sequences (CDSs) from a (multi)-embl or#\n".
   "\t# -genbankfile and writes them to a multi-fasta file. The ID  #\n".
   "\t# line includes either the locus tag, if that's not available #\n".
   "\t# protein id, gene or product (in this order). Pseudogenes    #\n".
   "\t# (tagged by /pseudo) are not included/skipped!               #\n".
   "\t# Optional, a file with locus tags can be given to extract    #\n".
   "\t# only these proteins (each locus tag in a new line).         #\n".
   "\t#                                                             #\n".
   "\t# version 0.2 (update 04.09.2012)    Andreas Leimbach         #\n".
   "\t# 24.05.2012        andreas.leimbach\@uni-wuerzburg.de         #\n".
   "\t###############################################################\n\n";

### ToDo
# - Get all features with translation if they don't have a locus_tag/protein_id/gene/product ... tag, maybe number through self (or just call ">unknown") and print out the translations! (Because otherwise features might be missed
# - Make a choice to print out other stuff from features? Maybe nucleotide sequences ... (plus flanking regions if needed)?
# - make a choice to use locus_tags or genes ... in the ID lines
# - Have a look at FeatureExtract from CBS for further ideas: http://www.cbs.dtu.dk/services/FeatureExtract/ (also available as download: http://www.cbs.dtu.dk/services/FeatureExtract/download.php or the virtual ribosome: http://www.cbs.dtu.dk/services/VirtualRibosome/


### Print usage if -h|--h|--help is given as argument or file-extension is not given
my $file = shift;
my $locus_tags = shift;
if (!defined $file) {
    die $usage;
} elsif ($file =~ m/-h/) {
    die $usage;
}


### Create a bioperl SeqIO object with the file
my $seqio_object = Bio::SeqIO->new(-file => "<$file");


### Print the input and output file which are processed and written, resp.
print "\nInput: $file\t";
$file =~ s/(.+)\.\w+$/$1/;
print "Output: $file.fasta\n";
open(FILE, ">$file.fasta");


### Get the protein fasta sequences from the translation tag in the embl/genbank file
while (my $seq_object = $seqio_object->next_seq) {
    foreach my $feat_object ($seq_object->get_SeqFeatures) {
	if ($feat_object->primary_tag eq 'CDS') {
	    if ($feat_object->has_tag('pseudo')) { # skip pseudogene, that don't have a translation!
		next;
	    }
	    if ($feat_object->has_tag('locus_tag')) { # Use a different order of id-tags here or an option which one to use (see ToDo above)?
		if (defined $locus_tags) {
		    open(LOCUS, "<$locus_tags");
		    while (my $locus = <LOCUS>) {
			chomp $locus;
			my ($feat) = $feat_object->get_tag_values('locus_tag'); # values always returned as ARRAYS!
			if ($locus eq $feat) {
			    print FILE ">", $feat_object->get_tag_values('locus_tag'), "\n";
			    print FILE $feat_object->get_tag_values('translation'), "\n";
			}
		    }
		    next;
		}
		print FILE ">", $feat_object->get_tag_values('locus_tag'), "\n";
	    } elsif ($feat_object->has_tag('protein_id') && !defined $locus_tags) {
		print FILE ">", $feat_object->get_tag_values('protein_id'), "\n";
	    } elsif ($feat_object->has_tag('gene') && !defined $locus_tags) {
		print FILE ">", $feat_object->get_tag_values('gene'), "\n";
	    } elsif  ($feat_object->has_tag('product') && !defined $locus_tags) {
		print FILE ">", $feat_object->get_tag_values('product'), "\n";
	    }
	    print FILE $feat_object->get_tag_values('translation'), "\n";
	}
    }
}
if (defined $locus_tags) {
    close LOCUS;
}
close FILE;

exit;
