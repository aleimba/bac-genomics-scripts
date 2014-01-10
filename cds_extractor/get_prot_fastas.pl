#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SeqFeatureI;  # bioperl module to handle features in a sequence

my $usage = "\n".
   "\t#################################################################\n".
   "\t# $0 seq-file (locus_tag-file)                  #\n". #$0 = program name
   "\t#                                                               #\n".
   "\t# Extracts all protein sequences (CDSs) from a (multi)-embl or  #\n".
   "\t# -genbankfile and writes them to a multi-fasta file. The fasta #\n".
   "\t# ID line includes either the locus tag (plus g=gene, p=product,#\n".
   "\t# o=organism; if existent), if that's not available protein id  #\n".
   "\t# (plus g=gene, p=product, o=organism; if existent), gene       #\n".
   "\t# (plus p=product, o=organism), product (plus o=organism), or   #\n".
   "\t# an internal CDS counter (in this order). Pseudogenes (tagged  #\n".
   "\t# by /pseudo) are not included (except in the CDS counter!).    #\n".
   "\t# CDS without locus tags are written to the error file          #\n".
   "\t# locus_tag_errors.txt (if the file already exists, errors will #\n".
   "\t# be appended to it)!                                           #\n".
   "\t# Optionally, a file with locus tags can be given to extract    #\n".
   "\t# only these proteins (each locus tag in a new line).           #\n".
   "\t#                                                               #\n".
   "\t# version 0.4 (update 06.02.2013)              Andreas Leimbach #\n".
   "\t# 24.05.2012                              aleimba[at]gmx[dot]de #\n".
   "\t#################################################################\n\n";

### ToDo
# - Make a choice to print out other stuff from features? Maybe nucleotide sequences ... (plus flanking regions if needed)? --> compare to extractfeat from the EMBOSS package!
# - Have a look at FeatureExtract from CBS for further ideas: http://www.cbs.dtu.dk/services/FeatureExtract/ (also available as download: http://www.cbs.dtu.dk/services/FeatureExtract/download.php or the virtual ribosome: http://www.cbs.dtu.dk/services/VirtualRibosome/


### Print usage if -h|--h|--help is given as argument or sequence file is not given
my $file = shift;
my $locus_list = shift;
if (!defined $file) {
    die $usage;
} elsif ($file =~ m/-h/) {
    die $usage;
}


### Create a bioperl SeqIO object with the file
my $seqio_object = Bio::SeqIO->new(-file => "<$file"); # no '-format' to leave to bioperl guessing


### Print the input and output file which are processed and written, resp.
print "\nInput: $file\t";
$file =~ s/(.+)\.\w+$/$1/;
print "Output: $file\_prot.fasta\n";
open(FILE, ">$file\_prot.fasta");


### Get the protein fasta sequences from the translation tag in the embl/genbank file
my $organism; # store the organism information of the genome file to include in each header, include plasmid names
my $no_locus_tag = 0; # print ONE error statement (i.e. only for the first) if no locus tag is found on a CDS
while (my $seq_object = $seqio_object->next_seq) {
    my $count = 0; # for CDSs without any of the features asked for below, counts also 'pseudo' CDSs
    foreach my $feat_object ($seq_object->get_SeqFeatures) {
	if ($feat_object->primary_tag eq 'source') {
	    ($organism) = $feat_object->get_tag_values('organism'); # values always returned as ARRAYS!
	    $organism =~ s/\s/\_/g;
	    if ($feat_object->has_tag('plasmid')) { # 'eval_sub' also possible, but 'has_tag' more elegantly
		my ($plasmid) = $feat_object->get_tag_values('plasmid');
		$plasmid =~ s/\s/\_/g;
		$organism = $organism.'-plasmid_'.$plasmid;
	    }
	}
	if ($feat_object->primary_tag eq 'CDS') {
	    $count++;
	    if ($feat_object->has_tag('pseudo')) { # skip pseudogenes, they don't include /translations!
		next;
	    }
	    if ($feat_object->has_tag('locus_tag')) {
		if (defined $locus_list) { # if locus_list is given only get those CDSs
		    open(LOCUS, "<$locus_list");
		    while (my $locus = <LOCUS>) {
			chomp $locus;
			my ($feat_locus) = $feat_object->get_tag_values('locus_tag'); # values always returned as ARRAYS!
			if ($locus eq $feat_locus) {
			    my $product = eval_sub($feat_object, 'product'); # subroutine to evaluate the existence of the tag and catch the error if not existent (for 'product' also replaces whitespaces with '_'); method 'has_tag' not possible, because variable $product would not be defined in print below (with 'eval_sub' the variable is just empty, but defined)
			    my $gene = eval_sub($feat_object, 'gene');
			    print FILE ">", $feat_object->get_tag_values('locus_tag'), " g=$gene p=$product o=$organism\n";
			    print FILE $feat_object->get_tag_values('translation'), "\n";
			}
		    }
		    next;
		}
		my $product = eval_sub($feat_object, 'product');
		my $gene = eval_sub($feat_object, 'gene');
		print FILE ">", $feat_object->get_tag_values('locus_tag'), " g=$gene p=$product o=$organism\n";
	    } elsif (defined $locus_list) { # in case a list of locus tags is given, no need to look at CDSs that doesn't have a locus tag
		next;
	    } elsif ($feat_object->has_tag('protein_id')) {
		my ($protein_id) = $feat_object->get_tag_values('protein_id');
		locus_tag_info($no_locus_tag, $organism, $protein_id); # subroutine to inform no locus tag is found
		my $product = eval_sub($feat_object, 'product');
		my $gene = eval_sub($feat_object, 'gene');
		print FILE ">$protein_id g=$gene p=$product o=$organism\n";
		$no_locus_tag++;
	    } elsif ($feat_object->has_tag('gene')) {
		my ($gene) = $feat_object->get_tag_values('gene');
		locus_tag_info($no_locus_tag, $organism, $gene);
		my $product = eval_sub($feat_object, 'product');
		print FILE ">$gene p=$product o=$organism\n";
		$no_locus_tag++;
	    } elsif  ($feat_object->has_tag('product')) {
		my ($product) = $feat_object->get_tag_values('product');
		locus_tag_info($no_locus_tag, $organism, $product);
		print FILE ">$product o=$organism\n";
		$no_locus_tag++;
	    } else { # if none of the above tags are existent use the internal counter
		$count = 'CDS'.$count;
		locus_tag_info($no_locus_tag, $organism, $count);
		print FILE ">$count o=$organism\n";
		$no_locus_tag++;
	    }
	    print FILE $feat_object->get_tag_values('translation'), "\n";
	}
    }
}
if (defined $locus_list) {
    close LOCUS;
}
if (-e 'locus_tag_errors.txt') {
    close ERR;
}
close FILE;

exit;

#############
#Subroutines#
#############

### Subroutine to evaluate the existence of tags and catch the error
sub eval_sub {
    my ($feat_object, $tag) = @_;
    my $value = '';
    eval {($value) = $feat_object->get_tag_values($tag);}; # catch error if tag doesn't exist
    if ($tag eq 'product') {
	$value =~ s/\s/_/g;
    }
    return $value;
}

### Subroutine to inform if locus_tags are missing
sub locus_tag_info {
    my ($no_locus_tag, $organism, $tag) = @_;
    my $err = 'locus_tag_errors.txt';
    if ($no_locus_tag == 0) { # Give only one warning per sequence to STDOUT
	print "###\'$organism\' has at least one CDS without a locus tag, the CDSs are written to the error file \'$err\'!\n";
	if (-e $err) {
	    print "### The error file \'$err\' already exists, further errors will be appended to the existing file!\n";
	}
    }
    open (ERR, ">>$err") or die "Failed to create file \'$err\': $!\n";
    print ERR "$organism\t$tag\n";
    close ERR;
    return 1;
}
