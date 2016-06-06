#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<prot_finder.pl> - search for query protein homologs in annotated
bacterial genomes with BLASTP

=head1 SYNOPSIS

C<perl prot_finder.pl -r report.blastp -s subject.faa E<gt> blast_hits.tsv>

=head1 DESCRIPTION

This script is intended to search for homologous proteins in
annotated bacterial genomes. For this purpose, a previous
L<B<BLASTP>|http://blast.ncbi.nlm.nih.gov/Blast.cgi>), either
L<legacy or plus|https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>,
needs to be run with query protein sequences against
a B<BLASTP> database of subject proteins (e.g. all proteins from
several I<Escherichia coli> genomes).

The script L<C<cds_extractor.pl>|/cds_extractor> (with options B<-p
-f>) can be used to create multi-FASTA protein files of all
non-pseudo CDS from RichSeq genome files to create the needed
subject B<BLASTP> database. Present locus tags will be used as FASTA
IDs, but see L<C<cds_extractor.pl>|/cds_extractor> for a description
of the format. Query protein sequences for the B<BLASTP> need a
unique FASTA ID.

The B<BLASTP> report file (option B<-r>), the subject protein
multi-FASTA file (option B<-s>), and optionally the query protein
(multi-)FASTA file (option B<-q>) are then given to
C<prot_finder.pl>. Significant B<BLASTP> subject hits are filtered
according to the given cutoffs (options B<-i>, B<-cov_q>, and
B<-cov_s>) and the result is printed as an informative tab-separated
result table to C<STDOUT>. To apply global identity/coverage cutoffs
to subject hits high-scoring pairs (HSPs) are tiled (see
L<http://www.bioperl.org/wiki/HOWTO:Tiling> and
L<http://search.cpan.org/dist/BioPerl/Bio/Search/Hit/GenericHit.pm>).
Additionally, the subject protein sequences with significant query
hits are written to result multi-FASTA files, named according to the
respective query FASTA IDs (optionally including the query sequence
with option B<-q>).

Optionally, L<B<Clustal Omega>|http://www.clustal.org/omega/> can be
called (option B<-a> with optional B<-p>) to create multiple
alignments (FASTA format) for each of the resulting multi-FASTA
files. These alignments can be used to calculate phylogenies e.g.
with L<B<RAxML>|http://sco.h-its.org/exelixis/software.html> or
L<B<MEGA>|http://www.megasoftware.net/>.

Run the script L<C<cds_extractor.pl>|/cds_extractor> (with options
B<-p -f>) and the B<BLASTP> manually or use the bash shell wrapper
script C<prot_finder_pipe.sh> (see below L</"EXAMPLES">) to execute
the whole pipeline including C<prot_finder.pl> (with optional option
B<-q>). For a description of the pipeline and additional options see
option B<-h> of the shell script. Be aware that some options in
C<prot_finder_pipe.sh> corresponding to options in C<prot_finder.pl>
have different names. If L<C<cds_extractor.pl>|/cds_extractor> is
used in the pipeline (option B<-f> of the shell script) the working
folder has to contain the annotated bacterial genome subject files
(in RichSeq format, e.g. EMBL or GENBANK format).

At last, the resulting tab-separated table can be given to the
script C<prot_binary_matrix.pl> to create a presence/absence matrix
of the query proteins for each genome. Again see option B<-h> of
C<prot_binary_matrix.pl> for additional info. The presence/absence
matrix can also be transposed with script C<transpose_matrix.pl>
(see its help with B<-h>). These presence/absence matri(x|ces) can
e.g. be loaded into L<B<iTOL>|http://itol.embl.de/> to associate the
data with a phylogenetic tree. Also, you can use
C<binary_group_stats.pl> to calculate presence/absence statistics
for groups of columns and not simply single columns of the matrix.
C<binary_group_stats.pl> also has a comprehensive manual with its
option B<-h>.

=head1 OPTIONS

=head2 Mandatory options

=over 22

=item B<-r>=I<str>, B<-report>=I<str>

Path to B<BLASTP> report/output

=item B<-s>=I<str>, B<-subject>=I<str>

Path to subject multi-FASTA protein sequence file (*.faa) created
with L<C<cds_extractor.pl>|/cds_extractor> (and its options B<-p
-f>), which was used to create the B<BLASTP> database

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-d>=I<str>, B<-dir_result>=I<str>

Path to result folder [default = query identity and coverage
cutoffs, './results_i#_cq#']

=item B<-f>, B<-force_dir>

Force output to an existing result folder, otherwise ask user to
remove content of existing folder. Careful, files from a previous
analysis might not be overwritten if different to current analysis.

=item B<-q>=I<str>, B<-query>=I<str>

Path to query (multi-)FASTA protein sequence file (*.faa) with
B<unique> FASTA IDs, which was used as query in the B<BLASTP>. Will
include each query protein sequence in the respective multi-FASTA
F<query-ID_hits.faa> result file.

=item B<-b>, B<-best_hit>

Give only the best hit (i.e. highest identity) for each subject
sequence if a subject has several hits with different queries

=item B<-i>=I<int>, B<-ident_cutoff>=I<int>

Query identity cutoff for significant hits (not including gaps), has
to be an integer number >= 0 and <= 100 [default = 70]

=item B<-cov_q>=I<int>, B<-cov_query_cutoff>=I<int>

Query coverage cutoff, has to be an integer >= 0 and <= 100 [default
= 70]

=item B<-cov_s>=I<int>, B<-cov_subject_cutoff>=I<int>

Subject/hit coverage cutoff, has to be an integer >= 0 and <= 100
[default = 0]

=item B<-a>, B<-align_clustalo>

Call L<B<Clustal Omega>|http://www.clustal.org/omega/> for multiple
alignment of each F<query-ID_hits.faa> result file

=item B<-p>=I<str>, B<-path_clustalo>=I<str>

Path to executable B<Clustal Omega> binary if not present in global
C<PATH> variable; requires option B<-a>

=item B<-t>=I<int>, B<-threads_clustalo>=I<int>

Number of threads for B<Clustal Omega> to use; requires option B<-a>
[default = all processors on system]

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 17

=item C<STDOUT>

The resulting tab-delimited output table with the significant
subject B<BLASTP> hits is printed to C<STDOUT>. Redirect (e.g. to a
file in the result directory, options B<-d -f>) or pipe into another
tool as needed (e.g. C<prot_binary_matrix.pl>).

=item F<./results_i#_cq#>

All result files are stored in a result folder

=item F<./results_i#_cq#/query-ID_hits.faa>

Multi-FASTA protein files of significant subject hits for each query
protein (named after the respective query FASTA ID), optionally
include the respective query protein sequence (with option B<-q>)

=item F<subject.faa.idx>

Index file of the subject protein file for fast sequence retrieval
(can be deleted if no further B<BLASTPs> are needed with these
subject sequences)

=item (F<./results_i#_cq#/queries_no_blastp-hits.txt>)

Lists all query sequence IDs without significant subject hits; with
option B<-b> includes also queries with significant hits but
I<without> a best blast hit for a subject

=item (F<./results_i#_cq#/clustal_omega.log>)

Optional log file of verbose B<Clustal Omega> C<STDOUT/STDERR> messages

=item (F<./results_i#_cq#/query-ID_aln.fasta>)

Optional B<Clustal Omega> multiple alignment of each
F<query-ID_hits.faa> result file in FASTA alignment format

=item (F<./results_i#_cq#/query-ID_tree.nwk>)

Optional B<Clustal Omega> NJ-guide tree in Newick format

=back

=head1 EXAMPLES

=head2 L<C<cds_extractor.pl>|/cds_extractor>

=over

=item C<for file in *.(gbk|embl); do perl cds_extractor.pl -i "$file" -p -f; done>

=item C<cat *.faa E<gt> subject.faa>

=item C<rm !(subject).faa>

=back

=head2 Legacy B<BLASTP>

=over

=item C<formatdb -p T -i subject.faa -n prot_finder_db>

=item C<blastall -p blastp -d prot_finder_db -i query.faa -o prot_finder.blastp -e 1e-10 -F F -s T -b 500>

=back

B<or>

=head2 B<BLASTP+>

=over

=item C<makeblastdb -dbtype prot -in subject.faa -out prot_finder_db>

=item C<blastp -db prot_finder_db -query query.faa -out prot_finder.blastp -evalue 1e-10 -seg no -use_sw_tback -num_alignments 500>

=back

=head2 C<prot_finder.pl>

=over

=item C<perl prot_finder.pl -r prot_finder.blastp -s subject.faa -cov_s 80 E<gt> blast_hits.tsv>

=back

B<or>

=over

=item C<perl prot_finder.pl -r prot_finder.blastp -s subject.faa -d result_dir -f -q query.faa -i 50 -cov_q 50 -b -a -p ~/bin/clustalo -t 6 E<gt> result_dir/blast_hits.tsv>

=back

=head2 All-in-one with bash script pipeline

=over

=item C<./prot_finder_pipe.sh -q query.faa -s subject.faa E<gt> blast_hits.tsv>

=back

B<or>

=over

=item C<./prot_finder_pipe.sh -q query.faa -f (embl|gbk) -d result_dir -p legacy -e 0 -t 12 -i 50 -c 50 -k 30 -b -a -o ~/bin/clustalo -m E<gt> result_dir/blast_hits.tsv>

=back

=head1 DEPENDENCIES

=over

=item L<B<BioPerl>|http://www.bioperl.org>>

Tested with B<BioPerl> version 1.006923

=item L<B<Clustal Omega>|http://www.clustal.org/omega/>

Tested with B<Clustal Omega> version 1.2.1

=back

=head1 VERSION

 0.7.1                                             update: 05-04-2016
 0.1                                                       03-09-2012

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 (GPLv3) of the
License, or (at your option) any later version.

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

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO; # BioPerl module to handle sequence input/output
use Bio::SearchIO; # BioPerl module to handle BLAST reports
use Bio::Index::Fasta; # BioPerl module to create an index for a multi-FASTA file for faster sequence retrieval

my $Cmdline = "$0 @ARGV"; # used call command

### Get the options with Getopt::Long
my $Blastp_Report_File; # path to BLASTP report/output file
my $Subject_File; # multi-FASTA protein file from 'cds_extractor.pl' (with options '-p -f') which was used to create the BLASTP DB (the subjects)
my $Result_Dir; # path to result folder; default is set below to '"./results_i".$Ident_Cutoff."_cq".$Cov_Query_Cutoff'
my $Opt_Force_Result_Dir; # force output to existing results folder
my $Query_File; # optionally, needed to include the query proteins in the result/hit multi-FASTA protein file for subsequent alignment
my $Opt_Best_Hit; # optionally, give only the best hit (i.e. highest identity) for each subject sequence, even if a subject sequence has several hits with different queries; if option not given report all subjects hits for each query
my $Ident_Cutoff = 70; # query identity cutoff (without gaps)
my $Cov_Query_Cutoff = 70; # query coverage cutoff
my $Cov_Subject_Cutoff = 0; # subject/hit coverage cutoff
my $Opt_Align_Clustal; # optionally, align the result sequences with Clustal Omega
my $Clustal_Path; # optionally, give path to the Clustal Omega binary
my $Clustal_Threads; # optionally, give number of threads Clustal Omega will use
my $VERSION = '0.7.1';
my ($Opt_Version, $Opt_Help);
GetOptions ('report=s' => \$Blastp_Report_File,
            'subject=s' => \$Subject_File,
            'dir_result=s' => \$Result_Dir,
            'force_dir' => \$Opt_Force_Result_Dir,
            'query=s' => \$Query_File,
            'best_hit' => \$Opt_Best_Hit,
            'ident_cutoff=i' => \$Ident_Cutoff,
            'cov_query_cutoff=i' => \$Cov_Query_Cutoff,
            'cov_subject_cutoff=i' => \$Cov_Subject_Cutoff,
            'align_clustalo' => \$Opt_Align_Clustal,
            'path_clustalo=s' => \$Clustal_Path,
            'threads_clustalo=i' => \$Clustal_Threads,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD, enforce mandatory options, and check options
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);

if (!$Blastp_Report_File || !$Subject_File) {
    my $warning = "\n### Fatal error: Mandatory options '-r' and '-s' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}
die "\n### Fatal error:\nBLASTP report file '$Blastp_Report_File' does not exist: $!\n" if (!-e $Blastp_Report_File);
die "\n### Fatal error:\nSubject multi-FASTA protein sequence file '$Subject_File' does not exist: $!\n" if (!-e $Subject_File);
die "\n### Fatal error:\nQuery multi-FASTA protein sequence file '$Query_File' does not exist: $!\n" if ($Query_File && !-e $Query_File);

if (($Ident_Cutoff < 0 || $Ident_Cutoff > 100) || ($Cov_Query_Cutoff < 0 || $Cov_Query_Cutoff > 100) || ($Cov_Subject_Cutoff < 0 || $Cov_Subject_Cutoff > 100)) {
    my $warning = "\n### Fatal error:\nAll cutoff options ('-i', '-cov_q', and '-cov_s') require an integer number >= 0 and <= 100 as value!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}

die "\n### Fatal error:\nOption '-p' requires the path to an executable Clustal Omega binary as value, not '$Clustal_Path'!\n" if ($Clustal_Path && !-x $Clustal_Path); # without '-p': presence of 'clustalo' in global PATH checked by system call at end of script
if (!$Opt_Align_Clustal && ($Clustal_Path || $Clustal_Threads)) {
    warn "\n### Warning: Options '-p' and/or '-t' set without their required option '-a', forcing option '-a'!\n";
    $Opt_Align_Clustal = 1;
}

if ($Clustal_Threads) {
    my $max_cpus = qx{nproc --all 2> /dev/null} || die "\n### Fatal error:\nCouldn't run Unix GNU command 'nproc' to determine the overall number of processors on the local system!\n"; # get max number of processors on the system; '||' because success exit code is zero (in '$?')
    chomp $max_cpus;
    if ($Clustal_Threads < 0 || $Clustal_Threads > $max_cpus) {
        my $warning = "\n### Fatal error: Option '-t' requires an integer number > 0 and <= $max_cpus as value, not '$Clustal_Threads'!\n";
        pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
    }
}

print STDERR "\nScript call command: $Cmdline\n"; # print call command after '-h|-v'



### Create result folder
if (!$Result_Dir) { # can't give default before 'GetOptions' in case cutoffs are set by the user
    $Result_Dir = "./results_i".$Ident_Cutoff."_cq".$Cov_Query_Cutoff;
} else {
    $Result_Dir =~ s/\/$//; # get rid of a potential '/' at the end of $Result_Dir path
}
if (-e $Result_Dir && !$Opt_Force_Result_Dir) {
    empty_dir($Result_Dir); # subroutine to empty a directory with user interaction
} elsif (!$Opt_Force_Result_Dir) {
    mkdir $Result_Dir;
}



### Parse the BLASTP report/output file; print queries, subject ID (e.g. locus tag) hits, and stats; or store in hash for option '-b'
print STDERR "Parsing BLASTP report '$Blastp_Report_File' and looking for significant hits according to cutoffs ...\n"; # run status of script
print "# subject_organism\tsubject_ID\tsubject_gene\tsubject_protein_desc\tquery_ID\tquery_desc\tquery_coverage [%]\tquery_identities [%]\tsubject/hit_coverage [%]\te-value of best HSP\tbit-score of best HSP\n"; # header for output table with significant BLASTP hits

my %Query_Acc; # hash to check if query accessions/IDs are unique (which they have to be because of %Blast_Hits below) and to store all query_acc for '@Queries_No_Blasthit' below
my %Blast_Hits; # hash of array to store significant BLASTP hits for sequence retrieval afterwards; query_acc as key and subject ID hits as anonymous array
my %Best_Subject_Hit; # for option '-b'; hash of hash to store subject/query info/stats for only the best hit for each subject sequence ID (key)

my $Blast_Report = new Bio::SearchIO(-file => "<$Blastp_Report_File", -format => 'blast'); # Bio::SearchIO object
while (my $result = $Blast_Report->next_result) { # Bio::Search::Result::GenericResult object; several query sequences possible ($result = entire analysis for a single query seq)
    my $hit_present = 0; # indicates if significant BLASTP hit found for a query
    my @subject_ids; # array to store ALL significant subject ID hits for each query, for %Blast_Hits anonymous array

    my $query_desc = $result->query_description;
    my $query_acc = $result->query_accession =~ s/\.$//r; # = query ID; rm a '.' if present at the end of the string (for non-NCBI FASTA headers); non-destructive modifier causes the result of the substitution to be returned instead of modifying $_ (see http://perldoc.perl.org/perlrequick.html#Search-and-replace)
    die "\n### Fatal error:\nQuery accession/ID '$query_acc' is not unique in the query file, but has to be. Please edit all repetitive occurences and rerun BLASTP or the bash script pipeline!\n" if ($Query_Acc{$query_acc});
    $Query_Acc{$query_acc} = 1;

    while (my $hit = $result->next_hit) { # Bio::Search::Hit::GenericHit object; several subject sequences in the database might have hits for a query
        my $perc_identity = $hit->frac_identical('query'); # ignores gaps; method will call (requires) BioPerls HSP tiling 'tile_hsps()' to get value (see http://search.cpan.org/dist/BioPerl/Bio/Search/Hit/GenericHit.pm and for 'tile_hsps()' http://search.cpan.org/~cjfields/BioPerl-1.6.924/Bio/Search/SearchUtils.pm)
        $perc_identity *= 100;
        my $query_cov = $hit->frac_aligned_query; # method requires hsp tiling
        $query_cov *= 100;
        my $hit_cov = $hit->frac_aligned_hit; # = subject coverage; method requires hsp tiling
        $hit_cov *= 100;

        # "significant" hit according to cutoffs
        if ($perc_identity >= $Ident_Cutoff && $query_cov >= $Cov_Query_Cutoff && $hit_cov >= $Cov_Subject_Cutoff) {
            $hit_present = 1;
            my $hit_id = $hit->name; # = subject_id
            my ($gene, $product, $organism) = split_fasta_header($hit_id, $hit->description); # subroutine to split the subject FASTA ID lines (see cds_extractor.pl with option '-f')
            $gene = '' if (!$gene); # empty string if gene tag doesn't exist for print; $product and $organism should always exist
            my $evalue = $hit->significance;
            $evalue =~ s/\,$//; # rm ',' from the end of the evalue

            if (!$Opt_Best_Hit) { # print all hits to STDOUT directly without option '-b'
                print "$organism\t$hit_id\t$gene\t$product\t$query_acc\t$query_desc\t$query_cov\t$perc_identity\t$hit_cov\t$evalue\t", $hit->bits, "\n";
                push(@subject_ids, $hit_id); # store all hits for current query

            } elsif ($Opt_Best_Hit) { # store only the best hit for each subject ID (need to be stored, not printed, to check all queries = $result); print to STDOUT is below
                if (!$Best_Subject_Hit{$hit_id} || $Best_Subject_Hit{$hit_id}->{'perc_identity'} < $perc_identity) { # either hit/subject ID doesn't exist yet or replace with hit of higher identity
                    $Best_Subject_Hit{$hit_id} = {'organism' => $organism,
                                                  'gene' => $gene,
                                                  'product' => $product,
                                                  'query_acc' => $query_acc,
                                                  'query_desc' => $query_desc,
                                                  'query_cov' => $query_cov,
                                                  'perc_identity' => $perc_identity,
                                                  'hit_cov' => $hit_cov,
                                                  'evalue' => $evalue,
                                                  'bit_score' => $hit->bits}; # hash of hash; subject/hit ID keys are unique with '-b'
                }
            }
        }
    }

    # only if significant hit for current query; needed after 'next_hit'-while to collect all (subject) hits
    $Blast_Hits{$query_acc} = \@subject_ids if ($hit_present && !$Opt_Best_Hit); # hash of array; the same hit/subject ID can be a hit for different queries (without option '-b'), thus subject IDs are not unique and hash of (anonymous) array data structure is suitable; done in the same way for option '-b' during print out below (not here to check all queries)
}



### Option '-b' given; print out only the best hit for each hit/subject sequence and store respective IDs in %Blast_Hits (as was already done without '-b' above)
if ($Opt_Best_Hit) {
    print STDERR "Printing only the best hit for each subject protein in the whole BLASTP analysis (option '-b') ...\n"; # run status of script
    my $skip = ''; # skip queries that have already been processed

    foreach my $hit_id (sort{lc $Best_Subject_Hit{$a}->{'query_acc'} cmp lc $Best_Subject_Hit{$b}->{'query_acc'}} keys %Best_Subject_Hit) { # sort hit/subject IDs (keys of %Best_Subject_Hit) by 'query_acc' to look at each query_acc only once by $skip-ing the others
        next if ($Best_Subject_Hit{$hit_id}->{'query_acc'} eq $skip);
        $skip = $Best_Subject_Hit{$hit_id}->{'query_acc'};

        my @subject_ids = sort{lc $Best_Subject_Hit{$a}->{'organism'} cmp lc $Best_Subject_Hit{$b}->{'organism'}} grep($Best_Subject_Hit{$_}->{'query_acc'} eq $Best_Subject_Hit{$hit_id}->{'query_acc'}, keys %Best_Subject_Hit); # get all hit/subject IDs for the current 'query_acc', sorted by 'organism'
        $Blast_Hits{$skip} = \@subject_ids; # store all hit/subject IDs in hash of array with query_acc as key ($skip here), compatible to above without '-b'

        foreach my $subject_id (sort @subject_ids) { # print best hits to STDOUT, corresponding to the print of all hits without '-b' above
            print "$Best_Subject_Hit{$subject_id}->{'organism'}\t".
                  "$subject_id\t".
                  "$Best_Subject_Hit{$subject_id}->{'gene'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'product'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'query_acc'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'query_desc'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'query_cov'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'perc_identity'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'hit_cov'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'evalue'}\t".
                  "$Best_Subject_Hit{$subject_id}->{'bit_score'}\n";
        }
    }
}

my @Queries_No_Blasthit = grep (!$Blast_Hits{$_}, sort{lc $a cmp lc $b} keys %Query_Acc); # get all 'query_acc' without a significant and (for '-b') best blast hit (can use '$hit_present' above only without '-b', but to make consistent use '%Blast_Hits' with and without '-b')
# With '-b' have to use the approach here in case a query protein has a significant hit but not the best hit in any subject, because then wouldn't be in the results (a significant query hit might also be overwritten by a higher identity subsequent hit to another query, thus $hit_present doesn't work with '-b')
if (keys %Query_Acc == @Queries_No_Blasthit) {
    rmdir $Result_Dir;
    die "\nNo significant BLASTP hits could be found, exiting!\n";
}


### Create index for multi-FASTA subject protein file for faster sequence retrieval; indeces have to be unique (which works fine for locus tags, the most probable subject IDs with cds_extractor.pl)
print STDERR "Indexing subject multi-FASTA file, '$Subject_File\.idx', for sequence retrieval ...\n"; # run status of script
my $Inx = Bio::Index::Fasta->new(-filename => "$Subject_File\.idx", -write_flag => 1); # see http://www.bioperl.org/wiki/HOWTO:Beginners#Indexing_for_Fast_Retrieval or http://www.bioperl.org/wiki/HOWTO:Local_Databases
$Inx->make_index($Subject_File); # by default the FASTA indexing code will use the string following the > character as a key, in this case the subject IDs



### Get the significant BLASTP hit sequences from the indexed multi-FASTA protein subject file and the query protein file (w/o index) and write them to the result dir
print STDERR "Using the index to retrieve subject protein sequences from significant BLASTP hits for each query "; # run status of script
my $Query_Seqioobj;
if ($Query_File) { # option '-q' with path to query multi-FASTA file; to retrieve respective query seq and write it as first seq into the hit multi-FASTA result files (*query*_hits.faa)
    print STDERR "including each query sequence from file '$Query_File' (option '-q') "; # run status of script
    $Query_Seqioobj = Bio::SeqIO->new(-file => "<$Query_File", -format => 'fasta'); # Bio::SeqIO object
}
print STDERR "...\n"; # run status of script

my @Fasta_Files; # store all hit result FASTA filenames for Clustal Omega
foreach my $query_acc (sort keys %Blast_Hits) {
    my $fasta_outfile = "$Result_Dir/$query_acc\_hits.faa";
    push (@Fasta_Files, $fasta_outfile);
    my $seqio_outobj = Bio::SeqIO->new(-file => ">$fasta_outfile"); # write a multi-FASTA file of hits for each query; format not needed, as everything is and should be FASTA anyway

    # get the corresponding query_acc sequence if option '-q'
    if ($Query_File) {
        my $query_found = 0;
        while (my $query_seqinobj = $Query_Seqioobj->next_seq) { # Bio::Seq object; index not needed should be small file
            if ($query_seqinobj->display_id eq $query_acc) {
                $query_found = 1;
                $seqio_outobj->write_seq($query_seqinobj); # print seq entry
                seek($Query_Seqioobj->_fh, 0, 0); # set filepointer back to zero for the next query_acc/ID
                last; # query accn/ID is found
            }
        }
        die "\n### Fatal error:\nQuery ID '$query_acc' not found in the query file '$Query_File' given with option '-q'. Sure this is the correct file which was used for the BLASTP?\n" if (!$query_found);
    }

    # get all hit-subject sequences for the query
    foreach my $subject_id (sort @{ $Blast_Hits{$query_acc} }) {
        my $subject_seqobj = $Inx->fetch($subject_id); # a Bio::Seq object; fetch subject seq from index
        die "\n### Fatal error:\nSubject ID '$subject_id' not found in the subject file '$Subject_File' given with option '-s'. Sure this is the correct file which was used for the BLASTP?\n" if (!$subject_seqobj);

        ## used to have the following out-commented code, but why not use original desc from cds_extractor?
        #my ($gene, $product, $organism) = split_fasta_header($subject_id, $subject_seqobj->desc);
        #if ($gene) {
            #$subject_seqobj->desc("$organism $gene"); # set the description of the FASTA ID line to a new one, if a gene name exists
        #} else {
            #$subject_seqobj->desc("$organism"); # w/o gene name
        #}

        $seqio_outobj->write_seq($subject_seqobj);
    }
}



### OPTIONAL method to extract the BLASTP hit protein sequences without BioPerl and without an index; works with out-commented sub 'read_fasta_entry'
#print STDERR "Retrieving subject protein sequences from significant BLASTP hits for each query "; # run status of script
#open (my $Subject_Fh, "<", $Subject_File);

#my $Query_Fh;
#if ($Query_File) { # option '-q' with path to query multi-FASTA file; to retrieve respective query seq and write it as first seq into the hit multi-FASTA result files (*query*_hits.faa)
    #print STDERR "including each query sequence from file '$Query_File' (option '-q') "; # run status of script
    #open ($Query_Fh, "<", $Query_File);
#}
#print STDERR "...\n"; # run status of script

#my @Fasta_Files; # store all hit result FASTA filenames for Clustal Omega
#foreach my $query_acc (sort keys %Blast_Hits) {
    #my $fasta_outfile = "$Result_Dir/$query_acc\_hits.faa";
    #push (@Fasta_Files, $fasta_outfile);
    #open (my $fasta_out_fh, ">", $fasta_outfile); # write a multi-FASTA file of hits for each query

    ## get the corresponding query_acc sequence if option '-q'
    #if ($Query_File) {
        #my $query_found = 0;
        #my $next_fasta_header; # for multi-line FASTA input files to store next entry header/ID line while parsing in subroutine 'read_fasta_entry'
        #while (<$Query_Fh>) {
            #chomp;
            #(my $seq_entry, $next_fasta_header) = read_fasta_entry($_, $Query_Fh, $next_fasta_header); # subroutine to read each FASTA seq entry of a multi-seq file separately; return header (->[0]) and seq (->[1]) as anonymous array in $seq_entry
            #if ($seq_entry->[0] =~ /^>$query_acc/) { # use '^' to anchor query acc/ID match (can't use ' ' as below for subjects, in case FASTA header has only one "word"/term after the '>')
                #$query_found = 1;
                #print $fasta_out_fh "$seq_entry->[0]\n$seq_entry->[1]\n\n"; # print seq entry
                #seek $Query_Fh, 0, 0; # set filepointer back to zero for the next query acc/ID
                #$. = 0; # set line number of seq file to 0 (seek doesn't do it automatically)
                #last; # query acc/ID found
            #}
        #}
        #die "\n### Fatal error:\nQuery sequence '$query_acc' not found in the query file '$Query_File' given with option '-q'. Sure this is the correct file which was used for the BLASTP?\n" if (!$query_found);
    #}

    ## get all hit-subject sequences for the current query
    #foreach my $subject_id (sort @{ $Blast_Hits{$query_acc} }) {
        #my $subject_found = 0;
        #my $next_fasta_header;
        #while (<$Subject_Fh>) {
            #chomp;
            #(my $seq_entry, $next_fasta_header) = read_fasta_entry($_, $Subject_Fh, $next_fasta_header); # subroutine
            #if ($seq_entry->[0] =~ /^>$subject_id /) { # use '^' and ' ' to force complete subject ID match (e.g. problem with ABU83972 'ECABU_c27750' and CFT073 'c2775')
                #$subject_found = 1;

                ### used to have the following out-commented code, but why not use original desc from cds_extractor?
                ## print FASTA header/ID for subject sequence
                ##$seq_entry->[0] =~ s/>.+\s(g=.+)$/$1/; # get rid of the subject ID for subroutine 'split_fasta_header' below
                ##my ($gene, $product, $organism) = split_fasta_header($subject_id, $seq_entry->[0]);
                ##print $fasta_out_fh ">$subject_id $organism ";
                ##if ($gene) {
                    ##print $fasta_out_fh "$gene\n";
                ##} else {
                    ##print $fasta_out_fh "\n";
                ##}
                ##print $fasta_out_fh "$seq_entry->[1]\n\n";

                #print $fasta_out_fh "$seq_entry->[0]\n$seq_entry->[1]\n\n";
                #seek $Subject_Fh, 0, 0; # for the next $subject_id
                #$. = 0;
                #last;
            #}
        #}
        #die "\n### Fatal error:\nSubject sequence '$subject_id' not found in the subject sequence file '$Subject_File' given with option '-s'. Sure this is the correct file which was used for the BLASTP?\n" if (!$subject_found);
    #}
    #close $fasta_out_fh;
#}
#close $Subject_Fh;
#close $Query_Fh if ($Query_File);



### Align with Clustal Omega if option '-a' is set
if ($Opt_Align_Clustal) {
    print STDERR "Starting Clustal Omega alignment (option '-a') with file\n"; # run status of script
    foreach my $fasta_file (@Fasta_Files) {
        print STDERR "  $fasta_file";
        my $out = $fasta_file =~ s/\_hits.faa$//r;
        my $clustal_call = " -i $fasta_file -o $out\_aln.fasta --verbose --guidetree-out=$out\_tree.nwk >> $Result_Dir/clustal_omega.log"; # redirect verbose STDOUT output to file (can't use clustalo option '-l' because will overwrite for each call)
        $clustal_call = " --threads=$Clustal_Threads" . $clustal_call if ($Clustal_Threads);
        if ($Clustal_Path) { # append path to Clustal Omega binary if option '-p'
            $clustal_call = $Clustal_Path . $clustal_call;
        } else { # otherwise 'clustalo' hopefully present in global path
            $clustal_call = 'clustalo' . $clustal_call;
        }
        system ($clustal_call) == 0 or die "\n### Fatal error:\nClustal Omega alignment with option '-a' for multi-FASTA output file '$fasta_file' could not be run. If Clustal Omega's binary 'clustalo' is not installed in your global \$PATH use option '-p' to give the path to the binary! $!\n";
    }
    print STDERR "\n";
}



### Final run status of script, state if queries without BLASTP hits
print STDERR "Result files were created in '$Result_Dir'"; # run status of script
if (@Queries_No_Blasthit) {
    my $no_blasthit_file = "$Result_Dir/queries_no_blastp-hits.txt";
    print STDERR ", including '$no_blasthit_file' listing all queries without a significant BLASTP hit"; # run status of script
    open (my $no_blasthit_fh, ">", "$no_blasthit_file");
    print $no_blasthit_fh join("\n", @Queries_No_Blasthit), "\n"; # print query_accs to file
    close $no_blasthit_fh;
}
print STDERR ".\n"; # run status of script
print STDERR "If no further BLASTPs with the subject sequences in file '$Subject_File' are needed the index file '$Subject_File.idx' can be deleted!\n";


exit;



###############
# Subroutines #
###############

### Subroutine to empty a directory with user interaction
sub empty_dir {
    my $dir = shift;
    print STDERR "\nDirectory '$dir' already exists! You can use either option '-d' to set a different result directory name, or do you want to replace the directory and all its contents [y|n]? ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
        unlink glob "$dir/*"; # remove all files in results directory
    } else {
        die "\nScript abborted!\n";
    }
    return 1;
}



### Read sequence entries from FASTA file
#sub read_fasta_entry {
    #my ($line, $fh, $next_fasta_header) = @_;

    ## possible multi-line seq in FASTA
    #my ($seq, $header);
    #if ($. == 1) { # first line of file
        #die "\n### Fatal error:\nNot a FASTA input file, first line of file should be a FASTA ID/header line and start with a '>':\n$line\n" if ($line !~ /^>/);
        #$header = $line;
    #} elsif ($next_fasta_header) {
        #$header = $next_fasta_header;
        #$seq = $line;
    #}
    #while (<$fh>) {
        #chomp;
        #if (/^>/) {
            #$next_fasta_header = $_; # store ID/header for next seq entry
            #return ([$header, $seq], $next_fasta_header); # return anonymous array with current header and seq
        #}
        #$seq .= $_; # concatenate multi-line seq
    #}
    #return ([$header, $seq], $next_fasta_header) if (eof);
    #return; # return undef
#}



### Subroutine to split the headers/IDs of the protein multi-FASTA files from 'cds_extractor.pl' and its '-f' option
sub split_fasta_header {
    my ($id, $desc) = @_;

    my ($gene, $product, $length, $organism, $ec); # $length and $ec not used
    $desc =~ s/^\s*//; # remove possible space at begin of desc string (might be introduced by BioPerl as explained below); needed in case no 'g=' but introduced space in front of 'p=', then split will not work below
    $desc =~ s/\s(?!p=|l=|o=|ec=)//g; # if a FASTA ID line is too long the BLAST report hit desc is over several lines and BioPerl will introduce space characters for a newline, thus get rid of these extra spaces with a lookahead in the regex but leave the spaces in the cds_extractor file format intact
    die "\n### Fatal error:\nFASTA 'annotation' from 'cds_extractor.pl' is not recognized on the following line. Please use 'cds_extractor.pl' with options '-p -f' to create your multi-FASTA protein subject files.\n>$id $desc\n" if ($desc !~ /(g=\w+)?\s?     # optional gene tag
                                                                         (p=\S+)?\s?     # optional product tag; \S+ instead of \w needed for non-alphanumeric characters (e.g. commas)
                                                                         l=\d+\.\.\d+\s? # location always included
                                                                         (o=\S+)?\s?     # optional organism tag
                                                                         (ec=\S+)?       # optional EC number tag
                                                                         /x); # for the format see 'cds_extractor.pl'

    foreach (split(/\s/, $desc)) {
        if (/g=(.*)$/) {
            $gene = $1;
        } elsif (/p=(.*)$/) {
            warn "\n### Warning:\nNo product annotation (p=) on line:\n>$id $desc\nProceeding ...\n" if (!$1);
            $product = $1;
            $product =~ tr/_/ /; # replace the '_' back to spaces, as this was changed in 'cds_extractor.pl'
        } elsif (/l=(\d+)\.\.(\d+)$/) { # $length not used
            $length = abs($2 - $1) + 1; # abs(stop - start) + 1
        } elsif (/o=(.*)$/) {
            $organism = $1; # don't replace the '_' back, no spaces might be better for phylogenetic programs
        } elsif (/ec=(.*)$/) { # $ec not used
            $ec = $1;
        } else {
            die "\n### Fatal error:\nFASTA 'annotation' is not recognized on the following line. Please use 'cds_extractor.pl' with options '-p -f' to create your multi-FASTA protein subject files.\n>$id $desc\n";
        }
    }
    return ($gene, $product, $organism);
}
