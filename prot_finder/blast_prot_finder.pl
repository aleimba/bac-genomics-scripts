#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

blast_prot_finder.pl                                       03-09-2012

=head1 SYNOPSIS

This script is intended to search for homologous proteins in annotated
bacterial genomes. Therefore, a previous BLASTP
(L<http://blast.ncbi.nlm.nih.gov/Blast.cgi>) needs to be run with
known query protein sequences against a BLAST database of subject
proteins (e.g. all proteins from several E. coli genomes).
For this purpose the script I<get_prot_fastas_for_blast_prot_finder.pl>
can be used to create multi-fasta files of all non-pseudo CDS from
genome sequence files to create the needed BLAST database.
The BLAST report file, the query protein fasta file, and the subject
multi-fasta file are then given to I<blast_prot_finder.pl>. Significant
BLAST hits are filtered and the result is written to a tab-separated
file. The subject hits are also concatenated in a multi-fasta file
for each query sequence (including the query sequence). Optionally,
Clustal-Omega (L<http://www.clustal.org/omega/>) can be called
(has to be in the C<$PATH> or change C<$clustal_call>) to create an
alignment (fasta format) for each of the concatenate multi-fasta files.
The alignments can then be used to calculate phylogenies. Use I<SeaView>
(L<http://pbil.univ-lyon1.fr/software/seaview.html>) to convert the
alignment format to Phylip for I<RAxML>
(L<http://sco.h-its.org/exelixis/software.html>).

Run the script I<get_prot_fastas_for_blast_prot_finder.pl> and the
BLASTPs manually or use the shell scripts I<blast_prot_finder*.sh> (see
below examples) to do the whole pipeline (including I<blast_prot_finder.pl>
and Clustal-Omega alignments) consecutively in one folder.

Additionally, the result file 'blast_hits.txt' can be given to the
script I<blastp_iTOL_binary.pl> to create a presence/abscence binary matrix
of the results. This comma-separated file can be loaded into iTOL
(L<http://itol.embl.de/>) to associate the data with a phylogenetic tree.

The Perl script runs on BioPerl (L<http://www.bioperl.org>).

=head1 OPTIONS

=head2 Mandatory options

=over

=item B<-r>, B<-report>          BLASTP report/output

=item B<-q>, B<-query>           Query sequence file [fasta format]

=item B<-s>, B<-subject>         Subject sequence file [fasta format]

=back

=head2 Optional options

=over

=item B<-h>, B<-help>          Help (perldoc POD)

=item B<-i>=I<int>, B<-ident>=I<int>   Identity cutoff [default 70]

=item B<-c>=I<int>, B<-cov>=I<int>     Coverage cutoff [default 70]

=item B<-a>, B<-align>           Call Clustal-Omega for alignment

=back

=head1 OUTPUT

=over 17

=item F<blast_hits.txt>

Contains the filtered BLAST hits, tab-separated

=item F<acc#_hits.fasta>

Multi-fasta protein files of subject hits for each query protein (with acc#),
includes the query protein sequence

=item F< *.idx>

Index file of the subject protein file for fast sequence retrieval (can be
deleted if no further BLASTs are needed with this subject sequences)

=item (F<acc#_aln.fasta>)

Optional, Clustal-Omega alignment of subject hits for a query

=item (F<acc#_tree.nwk>)

Optional, Clustal-Omega NJ-guide tree in Newick format

=back

=head1 EXAMPLES

=over

=item C<perl blast_prot_finder.pl -q query_prot.fasta -s subject_prot.fasta
-r blast-report.out -i 50 -c 50>

=item C<perl blast_prot_finder.pl -q query_prot.fasta -s subject_prot.fasta
-r blast-report.out -a>

=item C<./blast_prot_finder_legacy.sh subject_file-extension query_prot.fasta
ident_cutoff cov_cutoff>

=item C<./blast_prot_finder_plus.sh subject_file-extension query_prot.fasta
ident_cutoff cov_cutoff>

=back

=head1 VERSION

0.4                                                update: 24-01-2013

=head1 AUTHOR

Andreas Leimbach                                aleimba[at]gmx[dot]de

=head1 LICENSE

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 (GPLv3) of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut


########
# MAIN #
########

use strict;
use warnings;
use Getopt::Long; # module to get options from the command line
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SearchIO; # bioperl module to handle blast reports
use Bio::Index::Fasta; # bioperl module to create an index for a multi-fasta file for faster retrieval of sequences



### ToDo
# - Implement 'tblastn' to search proteins in unannotated/bad annotated genomes (subject of the BLAST DB would be (multi-)nucleotide fasta)? --> Get the coords for the hits: look for overlapping CDSs (in annotated genomes), translate the nucleotide sequence and find the nearest Met if hit not 100% (see bioperl HowTo Beginners: '$seq_obj->translate(-complete => 1);' and set the codon table to bacteria!)
# - Use 'needle' from the EMBOSS package to do global Needleman-Wunsch alignments to calculate the coverage? So far the corresponding shell-scripts use the 'calculate an optimal Smith-Waterman alignment' option of BLAST (legacy: -s T. blast+: -use_sw_tback) as suggested in Moreno-Hagelsieb et al 2008


### Get the options with Getopt::Long, works also abbreviated and with two "--": -r, --r, -report ...
my $usage = "perldoc $0";
my $report = ''; # name of the blast report/output file
my $query_file = ''; # needed to include the query proteins in the result/hit multi-fasta protein file for subsequent alignment
my $subject_file = ''; # multi-fasta protein file from 'get_prot_fastas_for_blast_prot_fnder.pl' which was used to create the BLAST DB (the subjects)
my $ident_cutoff = 70; # stores the identity cutoff, or use standard value 70%
my $cov_cutoff = 70; # stores the coverage cutoff, or use standard value 70%
my $align = ''; # Optionally, align the sequences with Clustal-Omega
my $help = ''; # run perldoc on the POD
GetOptions ('report=s' => \$report, 'query=s' => \$query_file, 'subject=s' => \$subject_file,'ident_cutoff:i' => \$ident_cutoff, 'cov_cutoff:i' => \$cov_cutoff, 'align' => \$align, 'help' => \$help);


### Run perldoc on POD
if (!$report || !$query_file || !$subject_file) {
    die system($usage);
} elsif ($help) {
    die system($usage);
}


### Parse the blast report/output file
my $hit_outfile = "blast_hits.txt"; # tab-separated file to store significant blast hits
file_exist($hit_outfile); # subroutine to test for file existence, as opened in append mode
open (HIT, ">>$hit_outfile") or die "Failed to create file \'$hit_outfile\': $!\n";
print HIT "Subject_organism\tSubject_gene\tSubject locus_tag\tSubject protein_function\tQuery_accession\tQuery_description\tQuery_coverage [%]\tIdentities [%]\tE-value\tBit-Score\n";
my %blast_hits; # store significant hits in here with the query_acc (key) and an array reference (several queries can have the same subject locus tag as hit) on hit subject locus tags (values) for retrieval afterwards
my $parser = new Bio::SearchIO(-file => "<$report", -format => 'blast'); # Bio::SearchIO object
while (my $result = $parser->next_result) { # several query sequences possible (result = entire analysis for a single query seq)
    my $no_blasthit = 0; # Report if no significant hits were found for a query
    my @locus_tags; # array to store the subject locus tags hits for each query, fed into %blast_hits via reference
    my $query_acc = $result->query_accession;
    $query_acc =~ s/\.$//; # rm a '.' if present at the end of the string (for non-NCBI fasta headers)
    my $query_desc = $result->query_description;
    my $query_length = $result->query_length;
    while (my $hit = $result->next_hit) { # several subject sequences in the database might have hits
        my $hit_locus_tag = $hit->name;
        my ($gene, $product, $organism) = split_ID($hit->description); # subroutine to split the subject fasta ID lines (see get_prot_fastas_for_blast_prot_finder.pl)
        while (my $hsp = $hit->next_hsp) { # each hit might have one or more hsps (the alignments shown in a blast report)
            my $hsp_length = $hsp->length('hit'); # reference/subject hsp length
            my $coverage = ($hsp_length/$query_length)*100;
            $coverage = sprintf("%.2f", $coverage); # round coverage to 2 digits after decimal point
            my $perc_identity = $hsp->percent_identity;
            $perc_identity = sprintf("%.2f", $perc_identity);
            if ($perc_identity >= $ident_cutoff && $coverage >= $cov_cutoff) { # report only significant BLAST hits
                $no_blasthit++;
                push(@locus_tags, $hit_locus_tag); # store significant subject hit in @locus_tags
                my $evalue = $hsp->evalue;
                $evalue =~ s/\,$//; # rm ',' from the end of the evalue
                print HIT "$organism\t";
                if ($gene =~ /.+/) { # print 'empty' tab if gene tag doesn't exist
                    print HIT "$gene\t";
                } else {
                    print HIT "\t";
                }
                print HIT "$hit_locus_tag\t$product\t$query_acc\t$query_desc\t$coverage\t$perc_identity\t$evalue\t", $hsp->bits, "\n";
            }
        }
    }
    $blast_hits{$query_acc} = \@locus_tags; # the same locus tag can be a hit for different queries, thus locus tags are not unique and an array reference data structure is needed
    if ($no_blasthit == 0) {
        print "\nNo significant BLASTP hit was found for query: $query_acc\n";
    }
}
close HIT;


### Create index for multi-fasta subject protein file for faster retrieval of protein sequences, indeces have to be unique (which works fine for locus tags)
my $inx = Bio::Index::Fasta->new(-filename => $subject_file . ".idx", -write_flag => 1);
$inx->make_index($subject_file); # by default the fasta indexing code will use the string following the > character as a key, in this case the locus tags


### Get the significant BLAST hits from the indexed multi-fasta protein subject file and the query protein file (w/o index)
my @fasta_files; # store all created fasta files
foreach my $query_acc (sort keys %blast_hits) {
    if (!scalar @{$blast_hits{$query_acc}}) { # skip a query if it has no subject hits
        next;
    }
    my $fasta_outfile = "$query_acc\_hits.fasta";
    file_exist($fasta_outfile);
    push (@fasta_files, $fasta_outfile);
    my $seqio_outobj = Bio::SeqIO->new(-file => ">>$fasta_outfile"); # write a multi-fasta file of hits for each query, format not needed, as everything is and should be fasta anyway
    my $query_seqioobj = Bio::SeqIO->new(-file => "<$query_file", -format => 'fasta'); # Bio::SeqIO object to retrieve the current query sequence and write it as the first seq into the hit multi-fasta file
    while (my $seq_inobj = $query_seqioobj->next_seq) { # Bio::Seq object; get the query protein seq
        if ($seq_inobj->display_id =~ /$query_acc/) { # get the query seq w/o index
            $seqio_outobj->write_seq($seq_inobj);
        }
    }
    foreach my $locus_tag (@{$blast_hits{$query_acc}}) { # locus tags for each query stored as array reference from @locus_tags above
        my $seq_obj = $inx->fetch($locus_tag); # a Bio::Seq object; fetch subject seq from index
        my ($gene, $product, $organism) = split_ID($seq_obj->desc);
        if ($gene =~ /.+/) {
            $seq_obj->desc("$organism $gene"); # set the description of the fasta ID line to the new one, only if a gene name exists
        } else {
            $seq_obj->desc("$organism");
        }
        $seqio_outobj->write_seq($seq_obj);
    }
}


### OPTIONAL method to extract the BLAST hit protein sequences without bioperl and without an index
# open (SUBJECT, "<$subject_file") or die "Failed to open file \'$subject_file\': $!\n";
# open (QUERY, "<$query_file") or die "Failed to open file \'$subject_file\': $!\n"; # include the query protein seqs in the result files
# my @fasta_files; # store all created fasta files
# foreach my $query_acc (sort keys %blast_hits) {
    # if (!scalar @{$blast_hits{$query_acc}}) { # skip a query if it has no subject hits
        # next;
    # }
    # my $fasta_outfile = "$query_acc\_hits.fasta";
    # file_exist($fasta_outfile);
    # push (@fasta_files, $fasta_outfile);
    # open (OUT, ">>$fasta_outfile") or die "Failed to create file \'$fasta_outfile\': $!\n";
    # while (my $line = <QUERY>) { # write the correct query protein sequence to the subject hit multi-fasta
        # if ($line =~ /$query_acc/) {
            # chomp $line;
            # print OUT "$line\n";
            # $line = <QUERY>;
            # while ($line !~ /^>/ && $line !~ /^$/) { # read in sequence until the next fasta header, or EOF
                # chomp $line;
                # print OUT "$line\n";
                # $line = <QUERY>;
                # if (eof) { # get out of the loop if the next $line is the end of the file
                    # print OUT "$line\n"; # still print the last line, as it might still have sequence
                    # last;
                # }
            # }
            # last; # jump out of the loop if locus tag is found (the rest of the file doesn't need to be parsed)
        # }
    # }
    # foreach my $locus_tag (@{$blast_hits{$query_acc}}) { # locus tags for each query stored as array reference from @locus_tags above
        # while (my $line = <SUBJECT>) {
            # if ($line =~ /$locus_tag/) {
                # chomp $line;
                # $line =~ s/>.+\s(g=.+)$/$1/; # get rid of the locus tag for the subroutine split_ID (actually here should 'my @desc = split (' ', $line);' work [omitting the split_ID sub] also as all the spaces are replaced by '_' in 'get_prot_fastas_for_blast_prot_fnder.pl'
                # my ($gene, $product, $organism) = split_ID($line);
                # print OUT ">$locus_tag $organism ";
                # if ($gene =~ /.+/) {
                    # print OUT "$gene\n";
                # } else {
                    # print OUT "\n";
                # }
                # $line = <SUBJECT>;
                # while ($line !~ /^>/ && $line !~ /^$/) { # read in sequence until the next fasta header, or EOF
                    # chomp $line;
                    # print OUT "$line\n";
                    # $line = <SUBJECT>;
                    # if (eof) { # get out of the loop if the next $line is the end of the file
                        # print OUT "$line\n"; # still print the last line, as it might still have sequence
                        # last;
                    # }
                # }
                # seek SUBJECT, 0, 0; # set filepointer for the filehandle SUBJECT back to zero for the next locus_tag
                # last;
            # }
        # }
    # }
    # close OUT;
# }
# close SUBJECT;
# close QUERY;


### State which files were created or warn if no BLAST hits were found at all
if (-s $hit_outfile < 150) { # smaller than just the header, which should be 133 bytes
    print "No significant BLAST hits could be found!\n";
    unlink $hit_outfile;
    exit;
} else {
    print "\n###########################################################################\n";
    print "The following files were created:\n";
    print "\tA summary of the BLAST results were written to \'$hit_outfile\'!\n";
}
print "\tThe protein sequences of the BLAST hits were written to:\n";
foreach my $fasta (@fasta_files) {
    print "\t\t\t\t$fasta\n";
}
print "\tIf no further BLASTs with these subject sequences is needed\n\tthe index file \'$subject_file.idx\' can be deleted!\n";
print "###########################################################################\n";


### Align with Clustal-Omega if option '--align' (-a|--a) is given
if ($align) {
    print "\nStarting Clustal-Omega alignment with file ...\n";
    foreach my $fasta (@fasta_files) {
    print "\t$fasta\n";
    my $out = $fasta;
    $out =~ s/\_hits.fasta//;
    my $clustal_call = "clustalo -i $fasta -o $out\_aln.fasta --verbose --guidetree-out=$out\_tree.nwk";
    system ($clustal_call);
    }
}


exit;


###############
# Subroutines #
###############

### Subroutine to split the header of the protein multi-fasta files from 'get_prot_fastas_for_blast_prot_finder.pl'
sub split_ID {
    my $hit_desc = shift;
    my ($gene, $product, $organism) = '';
    if ($hit_desc =~ /^(g=.*)(p=.*)(o=.*)$/) { # If a product tag is too long, BLAST will introduce a space in the report, thus cannot split "@desc = split (' ', $hit_desc);" and have to use regex instead
        $gene = $1;
        $product = $2;
        $organism = $3;
        $gene =~ s/^g=//;
        $gene =~ s/\s//g; # get rid of all optionally introduced spaces
        $product =~ s/^p=//;
        $product =~ s/\s//g;
        $product =~ tr/_/ /; # replace the '_' back to spaces, as this was changed in the script 'get_prot_fastas_for_blast_prot_fnder.pl'
        $organism =~ s/^o=//; # don't replace the '_' back, no spaces might be better for phylogenetic programs
        $organism =~ s/\s//g;
    }
    return ($gene, $product, $organism);
}


### Subroutine to test for file existence and delete old file, as the blast hits are appended to files
sub file_exist {
    my $file = shift;
    if (-e $file) {
        print "\nThe result file \'$file\' exists already and will be overwritten!\n";
        unlink $file;
    }
}
