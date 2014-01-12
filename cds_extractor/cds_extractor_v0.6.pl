#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

cds_extractor.pl                                           24-05-2012

=head1 SYNOPSIS

C<perl cds_extractor.pl -i seq_file.embl -p>

=head1 DESCRIPTION

Extracts protein or DNA sequences of CDS features from a
(multi)-embl or -genbankfile and writes them to a multi-fasta file.
The fasta ID line includes either the locus tag (plus g=gene,
p=product, o=organism; if existent), if that's not available protein
id (plus g=gene, p=product, o=organism; if existent), gene (plus
p=product, o=organism), product (plus o=organism), or an internal
CDS counter (in this order). The organism info includes also
possible plasmid names. Pseudogenes (tagged by /pseudo) are not
included (except in the CDS counter).

CDSs without locus tags are written to an error file
I<locus_tag_errors.txt> (if the file already exists, errors will
be appended to it)!

Optionally, a file with locus tags can be given to extract only
these proteins (each locus tag in a new line).

The Perl script runs on BioPerl (L<http://www.bioperl.org>).

=head1 OPTIONS

=head2 Mandatory options

=over 23

=item B<-i>=I<str>, B<-input>=I<str>

Input RichSeq sequence file including annotation (embl or genbank)

=item B<-p>, B<-protein>

Extract protein sequence for each CDS feature, excludes
option B<-n>

=item B<-n>, B<-nucleotide> 

Extract nucleotide sequence for each CDS feature, excludes
option B<-p>

=back

=head2 Optional options

=over 28

=item B<-h>, B<-help>

Run Perldoc on POD

=item B<-l>=I<str>, B<-locustag_list>=I<str>

List of locus tags to extract only those

=item B<-u>=I<int>, B<-upstream>=I<int>

Include given number of flanking nucleotides upstream of each
CDS feature, forces option B<-n>

=item B<-d>=I<int>, B<-downstream>=I<int>

Include given number of flanking nucleotides downstream of each
CDS feature, forces option B<-n>

=item B<-f>, B<-full_header>

Include full ID header for downstream I<blast_prot_finder.pl>
analysis (all IDs include 'g=', 'p=', and 'o=' for parsing)

=back

=head1 OUTPUT

=over 23

=item F<*_cds_aa.fasta>

Multi-fasta file of CDS protein sequences

=item F<*_cds_nucl.fasta>

Multi-fasta file of CDS DNA sequences

=item (F<locus_tag_errors.txt>)

Lists CDS features without locus tags

=back

=head1 EXAMPLES

=over

=item C<perl cds_extractor.pl -i seq_file.gbk -p -l
locus_tags.txt>

=item C<perl cds_extractor.pl -i seq_file.embl -n -l
locus_tags.txt -u 100 -d 20>

=item C<perl cds_extractor.pl -i seq_file.gbk -p -f>

=back

=head1 VERSION

0.6                                                update: 06-06-2013

=head1 AUTHOR

Andreas Leimbach                                aleimba[at]gmx[dot]de

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 (GPLv3) of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see L<http://www.gnu.org/licenses/>.

=cut


########
# MAIN #
########

use warnings;
use strict;
use Getopt::Long; # module to get options from the command line
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::Seq; # bioperl module to play with the sequence and its features
use Bio::SeqFeatureI; # bioperl module to handle features in a sequence



### Get the options with Getopt::Long, works also abbreviated and with two "--": -i, --i, -input ...
my $seq_file = ''; # RichSeq sequence file including feature annotation
my $protein = ''; # extract protein sequences for each CDS feature; excludes option '-n'
my $nucleotide = ''; # extract nucleotide sequences for each CDS feature; excludes option '-p'
my $upstream = 0; # include given number of flanking nucleotides upstream of each CDS feature; forces option '-n'
my $downstream = 0; # include given number of flanking nucleotides downstream of each CDS feature; forces option '-n'
my $locustag_list; # list of locus tags to extract only those
my $full_header; # include a full ID header for downstream 'blast_prot_finder.pl' analysis (needed for regex in subroutine 'split_ID')
my $help = ''; # run perldoc on POD
GetOptions ('input=s' => \$seq_file, 'protein' => \$protein, 'nucleotide' => \$nucleotide, 'upstream:i' => \$upstream, 'downstream:i' => \$downstream, 'locustag_list:s' => \$locustag_list, 'full_header' => \$full_header, 'help' => \$help);



### Enforce mandatory option '-n' if '-u' or '-d' are given
if ($upstream || $downstream) {
    $nucleotide = 1;
    undef $protein;
}



### Run perldoc on POD for option help or if essential options are not given
my $usage = "perldoc $0";
if (!$seq_file) {
    die system($usage);
} elsif (!$protein && !$nucleotide) {
    die system($usage);
} elsif ($protein && $nucleotide) {
    die system($usage);
} elsif ($help) {
    die system($usage);
}



### Create a bioperl SeqIO object with the file
my $seqio_object = Bio::SeqIO->new(-file => "<$seq_file"); # no '-format' to leave to bioperl guessing



### Print the input and output file which are processed and written, respectively
print "\nInput: $seq_file\t";
my $ori_filename = $seq_file; # store original filename to give error statements
$seq_file =~ s/(.+)\.\w+$/$1/;
if ($protein) {
    $seq_file = $seq_file.'_cds_aa.fasta';
    print "Output: $seq_file\n";
    open (OUT, ">$seq_file") or die "Failed to create file \'$seq_file\': $!\n";
} elsif ($nucleotide) {
    $seq_file = $seq_file.'_cds_nucl.fasta';
    print "Output: $seq_file\n";
    open (OUT, ">$seq_file") or die "Failed to create file \'$seq_file\': $!\n";
}



### Get the list of locus tags if $locustag_list defined
my %locus_tags; # store locus tags and indicate if present in seq_file
if (defined $locustag_list) {
    open (LOCUS, "<$locustag_list") or die "Failed to open file \'$locustag_list\': $!\n";
    while (my $locus = <LOCUS>) {
	chomp $locus;
	$locus_tags{$locus} = 0; # changes to 1 if the locus tag was found in seq_file
    }
    close LOCUS;
}



### Get the protein fasta sequences from the translation tag or the nucleic acid subseq
my $organism; # store the organism information of the genome file to include in each header, include plasmid names
my $no_locus_tag = 0; # print ONE error statement (i.e. only for the first) if no locus tag is found on a CDS
my %double_id; # control if a locus tag or protein id is double in the file, then exit (they should be unambiguous); see sub 'control_double'
while (my $seq_object = $seqio_object->next_seq) {
    my $count = sprintf("%04d", 0); # for CDSs without any of the features asked for below, counts also 'pseudo' CDSs
    my @feat_objects = $seq_object->get_SeqFeatures; # slurp all features from the seq to check CDS tag annotation
    if (grep ($_->primary_tag eq 'CDS', @feat_objects) == 0) {
	close OUT;
	unlink $seq_file or die "Couldn't remove empty result file \'$seq_file\': $!\n";
	die "\n###Fatal error\nSequence with no CDS annotation in \'$ori_filename\'!\nExiting program!\n\n";
    }
    foreach my $feat_object (@feat_objects) { # a Bio::SeqFeatureI object
	if ($feat_object->primary_tag eq 'source') {
	    $organism = eval_sub($feat_object, 'organism');
	    $organism =~ s/\s/\_/g;
	    if ($feat_object->has_tag('plasmid')) { # 'eval_sub' also possible, but 'has_tag' more elegantly
		my ($plasmid) = $feat_object->get_tag_values('plasmid'); # values always returned as ARRAYS; might use 'join(' ', $feat_object->get_tag_values($tag))' instead
		$plasmid =~ s/\s/\_/g;
		$organism = $organism.'-plasmid_'.$plasmid;
	    }
	}
	if ($feat_object->primary_tag eq 'CDS') {
	    $count++;
	    if ($feat_object->has_tag('pseudo')) { # skip pseudogenes, they don't include '/translation'
		next;
	    }
	    if ($feat_object->has_tag('locus_tag')) {
		control_double($feat_object->get_tag_values('locus_tag'), 'Locus tag'); # subroutine to control if locus tag already exists
		if (defined $locustag_list) { # if locus_list is given only get those CDSs
		    if (grep ($locus_tags{$_} == 1, keys %locus_tags) == keys %locus_tags) { # scalar context
			last; # skip residual feature objects if all locus tags are found
		    }
		    my ($feat_locus) = $feat_object->get_tag_values('locus_tag'); # values always returned as ARRAYS
		    my ($locus) = grep (/$feat_locus/, keys %locus_tags);
		    if ($locus) {
			$locus_tags{$locus} = 1;
			my $product = eval_sub($feat_object, 'product'); # subroutine to evaluate the existence of the tag and catch the error if not existent (for 'product' also replaces whitespaces with '_'); method '->has_tag' not possible, because variable $product would not be defined in print below (with 'eval_sub' the variable is empty, but defined)
			my $gene = eval_sub($feat_object, 'gene');
			print OUT ">$feat_locus g=$gene p=$product o=$organism\n";
			print_seq($feat_object, $seq_object); # subroutine to print the protein or nucleic sequence
		    }
		    next; # jump to the next feature object
		}
                ### '$locustag_list' is not defined:
		my $product = eval_sub($feat_object, 'product');
		my $gene = eval_sub($feat_object, 'gene');
		print OUT ">", $feat_object->get_tag_values('locus_tag'), " g=$gene p=$product o=$organism\n";
	    } elsif (defined $locustag_list) { # in case a list of locus tags is given, no need to look at CDSs that don't have a locus tag
		next; # jump to the next feature object
	    } elsif ($feat_object->has_tag('protein_id')) {
		my ($protein_id) = $feat_object->get_tag_values('protein_id');
		control_double($protein_id, 'Protein id');
		locus_tag_info($no_locus_tag, $organism, $protein_id); # subroutine to inform no locus tag is found
		my $product = eval_sub($feat_object, 'product');
		my $gene = eval_sub($feat_object, 'gene');
		print OUT ">$protein_id g=$gene p=$product o=$organism\n";
		$no_locus_tag++;
	    } elsif ($feat_object->has_tag('gene')) {
		my ($gene) = $feat_object->get_tag_values('gene');
		locus_tag_info($no_locus_tag, $organism, $gene);
		my $product = eval_sub($feat_object, 'product');
		if ($full_header) { # option '-f' induces full ID print out, needed for 'blast_prot_finder.pl'
		    print OUT ">$gene g=$gene p=$product o=$organism\n";
		} else {
		    print OUT ">$gene p=$product o=$organism\n";
		}
		$no_locus_tag++;
	    } elsif  ($feat_object->has_tag('product')) {
		my ($product) = $feat_object->get_tag_values('product');
		locus_tag_info($no_locus_tag, $organism, $product);
		if ($full_header) {
		    print OUT ">$product g= p=$product o=$organism\n";
		} else {
		    print OUT ">$product o=$organism\n";
		}
		$no_locus_tag++;
	    } else { # if none of the above tags are existent use the internal counter
		my $cds_count = 'CDS'.$count;
		locus_tag_info($no_locus_tag, $organism, $cds_count);
		if ($full_header) {
		    print OUT ">$cds_count g= p= o=$organism\n";
		} else {
		    print OUT ">$cds_count o=$organism\n";
		}
		$no_locus_tag++;
	    }
	    print_seq($feat_object, $seq_object);
	}
    }
    if (defined $locustag_list) {
	if (grep ($locus_tags{$_} == 1, keys %locus_tags) == keys %locus_tags) {
	    last; # skip residual sequence objects if all locus tags are found
	}
    }
}
if (-e 'locus_tag_errors.txt') {
    close ERR;
}
close OUT;



### Print locus tags that were not found in seq_file
if (defined $locustag_list) {
    my @missed = grep ($locus_tags{$_} == 0, keys %locus_tags);
    if (@missed) {
	print "### The following locus tags were not found in \'$ori_filename\':\n";
	foreach (sort @missed) {
	    print "$_\t";
	}
	print "\n";
    }
}


exit;


###############
# Subroutines #
###############

### Control if a locus tag or protein id is double in the seq_file and exit (they should be unambiguous) 
sub control_double {
    my ($id, $type) = @_;
    if ($double_id{$id}) {
	die "\n###Fatal error!\n$type \'$id\' exists at least two times in \'$organism\' of \'$ori_filename\', but should be unambiguous! Please modify all repetitive occurences.\nExiting program!\n\n";
    } else {
	$double_id{$id} = 1;
    }
    return 1;
}



### Evaluate the existence of tags and catch the error, also replace spaces with '_' for proteins
sub eval_sub {
    my ($feat_object, $tag) = @_;
    my $value = '';
    eval {($value) = $feat_object->get_tag_values($tag);}; # catch error if tag doesn't exist; values always returned as ARRAYS
    if ($tag eq 'product') {
	$value =~ s/\s/_/g;
    }
    return $value;
}



### Inform if locus_tags are missing
sub locus_tag_info {
    my ($no_locus_tag, $organism, $tag) = @_;
    my $err = 'locus_tag_errors.txt';
    if ($no_locus_tag == 0) { # Give only one warning per sequence to STDOUT
	print "###\'$organism\' has at least one CDS without a locus tag! The respective CDSs are written to the error file \'$err\'!\n";
	if (-e $err) {
	    print "### The error file \'$err\' already exists, the current errors will be appended to the existing file!\n";
	}
	open (ERR, ">>$err") or die "Failed to create file \'$err\': $!\n";
    }
    print ERR "$organism\t$tag\n";
    return 1;
}



### Print the protein or nucleic sequence to the result file
sub print_seq {
    my ($feat_object, $seq_object) = @_;
    if ($protein) {
	print OUT $feat_object->get_tag_values('translation'), "\n";
    } elsif ($nucleotide) {
	if ($feat_object->strand == 1) { # feature on leading strand
	    my $start = $feat_object->start - $upstream;
	    my $stop =  $feat_object->end + $downstream;
	    if ($start < 0) { # 'seq_obj->subseq' doesn't work with negative numbers, but with numbers > lenghts of seq_obj
		$start = $seq_object->length - abs($start);
		$stop = $seq_object->length + $stop;
	    }
	    my $subseq = $seq_object->subseq($start, $stop);
	    print OUT substr($subseq, 0, $upstream); # upstream bases, not possible to include in '$feat_obj->spliced_seq'
	    print OUT $feat_object->spliced_seq->seq; # to include also features with "join" locations; '$feat_object->spliced_seq' returns a Bio::Seq object thus call funtion '->seq' for the seq string
	    print OUT substr($subseq, length($subseq) - $downstream, $downstream), "\n"; # downstream bases
	} elsif ($feat_object->strand == -1) { # feature on lagging strand
	    my $start = $feat_object->start - $downstream;
	    my $stop = $feat_object->end + $upstream;
	    if ($start < 0) {
		$start = $seq_object->length - abs($start);
		$stop = $seq_object->length + $stop;
	    }
	    my $trunc_obj = $seq_object->trunc($start, $stop); # a Bio::Seq object, needed for revcom
	    my $rev_obj = $trunc_obj->revcom; # a Bio::Seq object
	    print OUT substr($rev_obj->seq, 0, $upstream);
	    print OUT $feat_object->spliced_seq->seq; # to include also features with "join" locations
	    print OUT substr($rev_obj->seq, $rev_obj->length - $downstream, $downstream), "\n";
	}
    }
    return 1;
}
