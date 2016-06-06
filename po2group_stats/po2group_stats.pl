#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<po2group_stats.pl> - categorize orthologs from Proteinortho5
output according to genome groups

=head1 SYNOPSIS

C<perl po2group_stats.pl -i matrix.proteinortho -d genome_fasta_dir/
-g group_file.tsv -p E<gt> overall_stats.tsv>

=head1 DESCRIPTION

Categorize the genomes in an ortholog/paralog output matrix (option
B<-i>) from a
L<B<Proteinortho5>|http://www.bioinf.uni-leipzig.de/Software/proteinortho/>
calculation according to group affiliations. The group
affiliations of the genomes are intended to get overall
presence/absence statistics for groups of genomes and not simply
single genomes (e.g. comparing 'marine', 'earth', 'commensal',
'pathogenic' etc. genome groups). Percentage inclusion (option
B<-cut_i>) and exclusion (option B<-cut_e>) cutoffs can be set to
define how strict the presence/absence of genome groups within an
orthologous group (OG) are defined. Of course groups can also hold
only single genomes to get single genome statistics. Group
affiliations are defined in a mandatory B<tab-delimited> group input
file (option B<-g>) with B<minimal two> and B<maximal four> groups.

Only alphanumeric (a-z, A-Z, 0-9), underscore (_), dash (-), and
period (.) characters are allowed for the B<group names> in the
group file to avoid downstream problems with the operating/file
system. As a consequence, also no whitespaces are allowed in these!
Additionally, B<group names>, B<genome filenames> (should be
enforced by the file system), and B<FASTA IDs> considering B<all>
genome files (mostly locus tags; should be enforced by
Proteinortho5) need to be B<unique>.

B<Proteinortho5> (PO) has to be run with option B<-singles> to
include also genes without orthologs, so-called singletons/ORFans,
for each genome in the PO matrix (see the
L<PO manual|http://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html>).
Additionally, option B<-selfblast> is recommended to enhance
paralog detection by PO.

To explain the logic behind the categorization, the following
annotation for example groups will be used. A '1' exemplifies a
group genome count in a respective OG E<gt>= the rounded inclusion
cutoff, a '0' a group genome count E<lt>= the rounded exclusion
cutoff. The presence and absence of OGs for the group affiliations
are structured in different categories depending on the number of
groups. For B<two groups> (e.g. A and B) there are five categories:
'A specific' (A:B = 1:0), 'B specific' (0:1), 'cutoff core' (1:1),
'underrepresented' (0:0), and 'unspecific'. Unspecific OGs have a
genome count for at least B<one> group outside the cutoffs
(exclusion cutoff E<lt> genome count E<lt> inclusion cutoff) and
thus cannot be categorized. These 'unspecific' OGs will only be
printed to a final annotation result file with option B<-u>. Overall
stats for all categories are printed to C<STDOUT> in a final
tab-delimited output matrix.

B<Three groups> (A, B, and C) have the following nine categories: 'A
specific' (A:B:C = 1:0:0), 'B specific' (0:1:0), 'C specific'
(0:0:1), 'A absent' (0:1:1), 'B absent' (1:0:1), 'C absent' (1:1:0),
'cutoff core' (1:1:1), 'underrepresented' (0:0:0), and 'unspecific'.

B<Four groups> (A, B, C, and D) are classified in 17 categories: 'A
specific' (A:B:C:D = 1:0:0:0), 'B specific' (0:1:0:0), 'C specific'
(0:0:1:0), 'D specific' (0:0:0:1), 'A-B specific' (1:1:0:0), 'A-C
specific' (1:0:1:0), 'A-D specific' (1:0:0:1), 'B-C specific'
(0:1:1:0), 'B-D specific' (0:1:0:1), 'C-D specific' (0:0:1:1), 'A
absent' (0:1:1:1), 'B absent' (1:0:1:1), 'C absent' (1:1:0:1), 'D
absent' (1:1:1:0), 'cutoff core' (1:1:1:1), 'underrepresented'
(0:0:0:0), and 'unspecific'.

The resulting group presence/absence (according to the cutoffs) can
also be printed to a binary matrix (option B<-b>) in the result
directory (option B<-r>), excluding the 'unspecific' category. Since
the categories are the logics underlying venn diagrams, you can also
plot the results in a venn diagram using the binary matrix (option
B<-p>). The 'underrepresented' category is exempt from the venn
diagram, because it is outside of venn diagram logics.

For an illustration of the logics have a look at the example venn
diagrams F<./pics/venn_diagram_logics.[svg|png]>.

There are two optional categories (which are only considered for the
final print outs and in the final stats matrix, not for the binary
matrix and the venn diagram): 'strict core' (option B<-co>) for OGs
where B<all> genomes have an ortholog. Of course all the 'strict
core' OGs are also included in the 'cutoff core' category ('strict
core' is identical to 'cutoff core' with B<-cut_i> 1 and B<-cut_e>
0). Option B<-s> activates the detection of 'singleton/ORFan' OGs
present in only B<one> genome. Depending on the cutoffs and number
of genomes in the groups, category 'underrepresented' includes most
of these singletons.

Additionally, annotation is retrieved from multi-FASTA files created
with L<C<cds_extractor.pl>|/cds_extractor>. See
L<C<cds_extractor.pl>|/cds_extractor> for a description of the
format. These files are used as input for the PO analysis and with
option B<-d> for C<po2group_stats.pl>. The annotations are printed
in category output files in the result directory.

Annotations are only pulled from one representative genome for each
category present in the current OG. With option B<-co> you can set a
specific genome for the representative annotation for category
'strict core'. For all other categories the representative genome is
chosen according to the order of the genomes in the group files,
depending on the presence in each OG. Thus, the best annotated
genome should be in the first group at the topmost position
(especially for 'cutoff core'), as well as the best annotated ones
at the top in all other groups.

In the result files, each orthologous group (OG) is listed in a row
of the resulting category files, the first column holds the OG
numbers from the PO input matrix (i.e. line number minus one). The
following columns specify the ID for each CDS, gene, EC number(s),
product, and organism are shown depending on their presence in the
CDS's annotation. The ID is in most cases the locus tag (see
L<C<cds_extractor.pl>|/cds_extractor>). If several EC numbers exist
for a single CDS they are separated by a ';'. If the representative
genome within an OG includes paralogs (co-orthologs) these will be
printed in the following row(s) B<without> a new OG number in the
first column.

The number of OGs in the category annotation result files are the
same as listed in the venn diagram and the final stats matrix.
However, since only annotation from one representative annotation is
used the CDS number will be different to the final stats. The final
stats include B<all> the CDS in this category in B<all> genomes
present in the OG in groups E<gt>= the inclusion cutoff (i.e. for
'strict core' the CDS for all genomes in this OG are counted). Two
categories are different, for 'unspecific' all unspecific groups are
included, for 'underrepresented' all groups E<lt>= the exclusion
cutoffs. This is also the reason, the 'pangenome' CDS count is
greater than the 'included in categories' CDS count in the final
stats matrix, as genomes in excluded groups are exempt from the CDS
counts for most categories. 'Included in categories' is the OG/CDS
sum of all non-optional categories ('*specific', '*absent', 'cutoff
core', 'underrepresented', and 'unspecific'), since the optional
categories are included in non-optionals. An exception to the
difference in CDS counts are the 'singletons' category where OG and
CDS counts are identical in the result files and in the overall
final output matrix (as there is only one genome), as well as in
group-'specific' categories for groups including only one genome.

At last, if you want the respective representative sequences for a
category you can first filter the locus tags from the result file
with Unix command-line tools:

C<grep -v "^#" result_file.tsv | cut -f 2 E<gt> locus_tags.txt>

And then feed the locus tag list to
L<C<cds_extractor.pl>|/cds_extractor> with option B<-l>.

As a final note, in the L<F<prot_finder>|/prot_finder> workflow is a
script, C<binary_group_stats.pl>, based upon C<po2group_stats.pl>,
which can calculate overall presence/absence statistics for column
groups in a delimited TEXT binary matrix (as with genomes here).

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Proteinortho (PO) result matrix (*.proteinortho or *.poff)

=item B<-d>=I<str>, B<-dir_genome>=I<str>

Path to the directory including the genome multi-FASTA PO input
files (*.faa or *.ffn), created with
L<C<cds_extractor.pl>|/cds_extractor>

=item B<-g>=I<str>, B<-groups_file>=I<str>

Tab-delimited file with group affiliation for the genomes with
B<minimal two> and B<maximal four> groups (easiest to create in a
spreadsheet software and save in tab-separated format). B<All>
genomes from the PO matrix need to be included. Group names can only
include alphanumeric (a-z, A-Z, 0-9), underscore (_), dash (-), and
period (.) characters (no whitespaces allowed either). Example
format with two genomes in group A, three in group B and D, and one
in group C:

group_A\tgroup_B\tgroup_C\tgroup_D
genome1.faa\tgenome2.faa\tgenome3.faa\tgenome4.faa
genome5.faa\tgenome6.faa\t\tgenome7.faa
\tgenome8.faa\t\tgenome9.faa

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-r>=I<str>, B<-result_dir>=I<str>

Path to result folder [default = inclusion and exclusion percentage
cutoff, './results_i#_e#']

=item B<-cut_i>=I<float>, B<-cut_inclusion>=I<float>

Percentage inclusion cutoff for genomes in a group per OG, has to be
E<gt> 0 and E<lt>= 1. Cutoff will be rounded according to the genome
number in each group and has to be E<gt> the rounded exclusion
cutoff in this group. [default = 0.9]

=item B<-cut_e>=I<float>, B<-cut_exclusion>=I<float>

Percentage exclusion cutoff, has to be E<gt>= 0 and E<lt> 1. Rounded
cutoff has to be E<lt> rounded inclusion cutoff. [default = 0.1]

=item B<-b>, B<-binary_matrix>

Print a binary matrix with the presence/absence genome group results
according to the cutoffs (excluding 'unspecific' category OGs)

=item B<-p>, B<-plot_venn>

Plot venn diagram from the binary matrix (except 'unspecific' and
'underrepresented' categories) with function C<venn> from B<R>
package B<gplots>, requires option B<-b>

=item B<-co>=(I<str>), B<-core_strict>=(I<str>)

Include 'strict core' category in output. Optionally, give a genome
name from the PO matrix to use for the representative output
annotation. [default = topmost genome in first group]

=item B<-s>, B<-singletons>

Include singletons/ORFans for each genome in the output, activates
also overall genome OG/CDS stats in final stats matrix for genomes
with singletons

=item B<-u>, B<-unspecific>

Include 'unspecific' category representative annotation file in
result directory

=item B<-a>, B<-all_genomes_overall>

Report overall stats for all genomes (appended to the final stats
matrix), also those without singletons; will include all overall
genome stats without option B<-s>

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The tab-delimited final stats matrix is printed to C<STDOUT>.
Redirect or pipe into another tool as needed.

=item F<./results_i#_e#>

All output files are stored in a results folder

=item F<./results_i#_e#/[*_specific|*_absent|cutoff_core|underrepresented]_OGs.tsv>

Tab-delimited files with OG annotation from a representative genome
for non-optional categories

=item (F<./results_i#_e#/[*_singletons|strict_core|unspecific]_OGs.tsv>)

Optional category tab-delimited output files with representative
annotation

=item (F<./results_i#_e#/binary_matrix.tsv>)

Tab-delimited binary matrix of group presence/absence results
according to cutoffs (excluding 'unspecific' category)

=item (F<./results_i#_e#/venn_diagram.pdf>)

Venn diagram for non-optional categories (except 'unspecific' and
'underrepresented' categories)

=back

=head1 EXAMPLES

=head2 L<C<cds_extractor.pl>|/cds_extractor>

=over

=item C<for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i [-p|-n]; done>

=back

=head2 L<B<Proteinortho5>|http://www.bioinf.uni-leipzig.de/Software/proteinortho/>

=over

=item C<proteinortho5.pl -graph [-synteny] -cpus=# -selfblast -singles -identity=50 -cov=50 -blastParameters='-use_sw_tback [-seg no|-dust no]' *.[faa|ffn]>

=back

=head2 C<po2group_stats.pl>

=over

=item C<perl po2group_stats.pl -i matrix.[proteinortho|poff] -d genome_fasta_dir/ -g group_file.tsv -r result_dir -cut_i 0.7 -cut_e 0.2 -b -p -co genome4.[faa|ffn] -s -u -a E<gt> overall_stats.tsv>

=back

=head1 DEPENDENCIES

=over

=item B<Statistical computing language L<R|http://www.r-project.org/>>

C<Rscript> is needed to plot the venn diagram with option B<-p>,
tested with version 3.2.2

=item B<L<gplots|https://cran.r-project.org/web/packages/gplots/index.html>>

Package needed for B<R> to plot the venn diagram with function
C<venn>. Tested with B<gplots> version 2.17.0.

=back

=head1 VERSION

 0.1.3                                             update: 06-06-2016
 0.1                                                       23-10-2015

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

my $Cmdline = "$0 @ARGV"; # used call command

### Get the options with Getopt::Long
my $PO_Matrix_File; # PO result matrix (*.proteinortho or *.poff)
my $Genome_Dir; # directory with genome multi-FASTA files (PO input)
my $Groups_File; # tab-separated file with group classification for the genomes
my $Result_Dir; # path to results folder; default is set below to '"./results_i".$Inclusion_Cutoff."c".$Exclusion_Cutoff'
my $Inclusion_Cutoff = 0.9; # inclusion percentage cutoff for genome count per OG (floating point number >0 && <=1)
my $Exclusion_Cutoff = 0.1; # exclusion percentage cutoff (floating point number >=0 && <1)
my $Opt_Binary; # print resulting binary matrix according to the group inclusion/exclusion results
my $Opt_Venn; # create venn diagram with R package 'gplots' function 'venn'; requires $Opt_Binary
my $Strict_Core; # report also strict core genome, i.e. OGs present in ALL genomes of PO matrix; either option flag (and use first genome in first group) or give genome file for result annotation print
my $Opt_Singletons; # report also singletons for each genome (outside of the group classifications)
my $Opt_Unspecific; # report also group-unspecific OGs ('exclusion < OG_group_genome_count < inclusion' for any group)
my $Opt_Genomes_Overall; # report overall stats also for genomes without singletons
my $VERSION = '0.1.3';
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$PO_Matrix_File,
            'dir_genome=s' => \$Genome_Dir,
            'groups_file=s' => \$Groups_File,
            'result_dir=s' => \$Result_Dir,
            'cut_inclusion=f' => \$Inclusion_Cutoff,
            'cut_exclusion=f' => \$Exclusion_Cutoff,
            'binary_matrix' => \$Opt_Binary,
            'plot_venn' => \$Opt_Venn,
            'core_strict:s' => \$Strict_Core,
            'singletons' => \$Opt_Singletons,
            'unspecific' => \$Opt_Unspecific,
            'all_genomes_overall' => \$Opt_Genomes_Overall,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD and enforce mandatory options
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);

if (!$PO_Matrix_File || !$Genome_Dir || !$Groups_File) {
    my $warning = "\n### Fatal error: Mandatory options '-i', '-d', or '-g' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}
$Genome_Dir =~ s/\/$//; # get rid of a potential '/' at the end of $Genome_Dir path
die "\n### Fatal error: Genome directory '$Genome_Dir' does not exist: $!\n" if (!-d $Genome_Dir);

if (($Inclusion_Cutoff <= 0 || $Inclusion_Cutoff > 1) || ($Exclusion_Cutoff < 0 || $Exclusion_Cutoff >= 1)) {
    my $warning = "\n### Fatal error:\nInclusion (0 < '-cut_i' <= 1) or exclusion (0 <= '-cut_e' < 1) cutoff(s) not chosen correctly!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}

print STDERR "Script call command: $Cmdline\n"; # print call command after '-h|-v'

if ($Opt_Venn && !$Opt_Binary) {
    warn "Option '-p' to draw the venn diagram set, but not it's required option '-b' for the binary matrix. Forcing option '-b'!\n";
    $Opt_Binary = 1;
}



### Parse genomes in groups file
print STDERR "Parsing group file '$Groups_File' "; # run status of script
open (my $Groups_File_Fh, "<", "$Groups_File");

my @Group_Names; # store group name column order
my %Genome_Groups; # hash of hash with group name as key, internal hash with resp. genomes anonymous array (key 'genomes'), or calculated integer inclusion/exclusion cutoffs (key 'inclusion' etc.)
while (<$Groups_File_Fh>) {
    chomp;
    next if (/^\s*$/); # skip empty lines

    # check groups file header and get group names
    if ($. == 1) { # header of groups file (first line)
        die "\n### Fatal error:\nGroups file '$Groups_File' does not have a correctly formatted header line with at least TWO and maximal FOUR tab-separated group name columns (without other whitespaces)!\n" if (!/^\S+\t\S+\t?\S*\t?\S*$/);
        foreach my $group_name (split(/\t/)) {
            check_file_system_conformity($group_name, 'group names'); # subroutine to check string only contains alphanumeric '0-9a-zA-Z', underscore '_', dash '-', and period '.' characters
            die "\n### Fatal error:\nGroup name '$group_name' is not unique in the group file '$Groups_File'! Please choose unique characters/names for the groups in the group file!\n" if (grep {/$group_name/} @Group_Names);
            push (@Group_Names, $group_name); # store groups array
        }
        next; # skip header line for subsequent code
    }

    # parse genomes in groups file
    die "\n### Fatal error:\nLine '$.' in the groups file '$Groups_File' does not have the correct format with maximal four tab-separated genome file names according to the groups in the header (without other whitespaces):\n$_\n" if (!/^\S*\t?\S*\t?\S*\t?\S*$/);
    my $column = -1; # column counter to associate genome filename (column) to correct group in header (@Group_Names array zero-based, thus start with '-1')
    foreach my $genome (split(/\t/)) {
        $column++;
        next if ($genome =~ /^$/); # skip empty columns in group file, in case the groups have uneven genome numbers
        die "\n### Fatal error:\nGenome '$genome' is not unique in the group file '$Groups_File'! However, the group file needs to contain all of the genomes from the header in the Proteinortho matrix '$PO_Matrix_File' and vice versa. Make sure the strings correspond exactly (also lower/uppercase) to the header of the matrix and the input files in directory '$Genome_Dir'!\n" if (check_genome_in_groups($genome)); # subroutine to check if a genome is present in %Genome_Groups
        push (@{ $Genome_Groups{$Group_Names[$column]}->{'genomes'} }, $genome); # push into anonymous array in hash of hash
    }
}
close $Groups_File_Fh;
print STDERR "with ", scalar map(@{ $Genome_Groups{$_}->{'genomes'} }, keys %Genome_Groups), " genomes ...\n"; # run status of script

foreach (@Group_Names) {
    die "\n### Fatal error:\nGroup '$_' does not contain a genome, please add a suitable genome filename or remove the group from the group file '$Groups_File'!\n" if (!$Genome_Groups{$_});
}



### Check $Strict_Core and set genome if not given
if ($Strict_Core) { # if a genome filename is given with option '-co'
    die "\n### Fatal error:\nThe given genome '$Strict_Core' with option '-co' is not present in the genomes of the group file '$Groups_File'. The string has to be exactly the same, please check!\n" unless (check_genome_in_groups($Strict_Core)); # subroutine to check if a genome is present in %Genome_Groups
} elsif (defined $Strict_Core && !$Strict_Core) { # if option '-co' is set without argument
    $Strict_Core = $Genome_Groups{$Group_Names[0]}->{'genomes'}[0];
    print STDERR "Genome for option '-co' was set automatically to the first genome in the first group, '$Strict_Core'!\n";
}



### Round inclusion and exclusion count cutoffs for each group
foreach (sort keys %Genome_Groups) {
    $Genome_Groups{$_}->{'inclusion'} = int((scalar @{ $Genome_Groups{$_}->{'genomes'} } * $Inclusion_Cutoff) + 0.5); # round positive number half up to integer ((s)printf rounds half to even), see e.g. https://stackoverflow.com/questions/178539/how-do-you-round-a-floating-point-number-in-perl
    print STDERR "Group '$_' with ", scalar @{ $Genome_Groups{$_}->{'genomes'} }, " genomes has a rounded inclusion cutoff of $Genome_Groups{$_}->{'inclusion'} genome(s) and ";

    $Genome_Groups{$_}->{'exclusion'} = int((scalar @{ $Genome_Groups{$_}->{'genomes'} } * $Exclusion_Cutoff) + 0.5);
    print STDERR "an exclusion cutoff of $Genome_Groups{$_}->{'exclusion'} genome(s)!\n";
    die "\n### Fatal error:\nThe rounded inclusion cutoff has to be greater than the exclusion cutoff for all groups, please choose appropriate cutoffs!\n" if ($Genome_Groups{$_}->{'inclusion'} <= $Genome_Groups{$_}->{'exclusion'});
}



### Create result folder
if (!$Result_Dir) { # can't give default before 'GetOptions' in case cutoffs are set by the user
    $Result_Dir = "./results_i".$Inclusion_Cutoff."_e".$Exclusion_Cutoff;
} else {
    $Result_Dir =~ s/\/$//; # get rid of a potential '/' at the end of $Result_Dir path
}
if (-e $Result_Dir) {
    empty_dir($Result_Dir); # subroutine to empty a directory with user interaction
} else {
    mkdir $Result_Dir;
}



### Parse OGs in input PO matrix
print STDERR "Parsing Proteinortho input matrix '$PO_Matrix_File' ...\n"; # run status of script
open (my $PO_Matrix_Fh, "<", "$PO_Matrix_File");

my @Genome_Files; # store genome filename column order
my %Ortho_Groups; # hash of hash with OG# as key, internal hash with genome filename as key and IDs (e.g. locus tags) in anonymous array
while (<$PO_Matrix_Fh>) {
    chomp;

    # check PO input file header and get genome file names
    if ($. == 1) { # header of PO matrix file (first line)
        die "\n### Fatal error:\nProteinortho input matrix '$PO_Matrix_File' does not have the mandatory header line, which starts with the first three tab-separated mandatory columns:\n# Species\tGenes\tAlg.-Conn.\n" if (!/^# Species\tGenes\tAlg\.-Conn\./);

        # get input genomes from PO matrix and check if fit to %Genome_Groups from $Groups_File
        foreach my $genome (split(/\t/)) {
            next if ($genome =~ /# Species|Genes|Alg\.-Conn\./); # non-genome header columns
            die "\n### Fatal error:\nGenome '$genome' present in the header of the Proteinortho input matrix '$PO_Matrix_File' but not in the groups file '$Groups_File'. However, the group file needs to contain all of the genomes from the matrix header and vice versa. Make sure the strings correspond exactly (also lower/uppercase) to the matrix header and the input files in directory '$Genome_Dir'!\n" unless (check_genome_in_groups($genome)); # subroutine to check if a genome is present in %Genome_Groups
            push (@Genome_Files, $genome);
        }
        die "\n### Fatal error:\nNot the same amount of genomes in the group file '$Groups_File' and the Proteinortho input matrix header '$PO_Matrix_File', but the group file needs to contain all of the genomes from the header in the Proteinortho matrix and vice versa!\n" if (@Genome_Files != map(@{ $Genome_Groups{$_}->{'genomes'} }, keys %Genome_Groups));
        next; # skip header line of PO matrix for subsequent code
    }

    # parse PO ortholog matrix
    my ($organisms, $genes, $alg_conn, @og) = split(/\t/); # $alg_conn not used
    my $organism_count = 0; # count number of organisms in OG to compare to '$organisms' (should be ok from PO)
    my $gene_count = 0; # count number of genes in OG to compare to '$genes' (should be ok from PO)
    my $column = -1; # column counter to associate column to correct genome file in OG (@Genome_Files array zero-based, thus start with '-1')
    foreach (@og) {
        $column++;
        die "\n###Fatal error:\nThere's an empty field in the Proteinortho input matrix '$PO_Matrix_File' in row '$.', column '$column+4'. This shouldn't happen, as all the fields of the matrix should either have a FASTA ID from the input genome files or a '*' in case the respective genome doesn't have an ortholog in this orthologous group!\n" if (/^$/);
        next if (/^\*$/); # skip columns without a FASTA ID (ortholog) in PO matrix (indicated by '*'), have to skip after '$column++' otherwise column won't fit
        $organism_count++;
        my @paralogs = split(',');
        foreach (sort {$a cmp $b} @paralogs) {
            $gene_count++;
            push(@{ $Ortho_Groups{$.-1}->{$Genome_Files[$column]} }, $_); # push into anonymous array in hash of hash ($.-1 = OG number, because of header)
        }
    }
    die "\n### Fatal error:\nThe indicated number of species ($organisms) does not fit to the counted number of organisms ($organism_count) in the orthologous group of line '$.' of the Proteinortho input matrix '$PO_Matrix_File'!\n" if ($organism_count != $organisms); # to check PO matrix (should be ok from PO)
    die "\n### Fatal error:\nThe indicated number of genes ($genes) does not fit to the counted number of orthologs/paralogs ($gene_count) in the orthologous group of line '$.' of the Proteinortho input matrix '$PO_Matrix_File'!\n" if ($gene_count != $genes); # to check PO matrix (should be ok from PO)
}
close $PO_Matrix_Fh;



### Parse annotations in genome multi-FASTAs
print STDERR "Parsing annotation from multi-FASTA CDS genome files in directory '$Genome_Dir' ...\n";
my %Annotation; # hash of hash to store the annotation of the genome files
foreach my $genome (@Genome_Files) { # @Genome_Files filled from PO matrix
    my $genome_file_path = "$Genome_Dir/$genome";
    check_file_exist($genome_file_path); # subroutine
    open (my $genome_fh, "<", "$genome_file_path");

    # parse anno in current genome file
    while (my $anno_line = <$genome_fh>) {
        chomp $anno_line;
        die "\n### Fatal error:\n'$genome_file_path' is not a FASTA input file. First line of the file should be a FASTA ID/header line and start with a '>':\n$anno_line\n" if ($anno_line !~ /^>/ && $. == 1);
        next if ($anno_line !~ /^>/); # skip non-ID FASTA lines

        my ($id, $gene, $product, $length, $organism, $ec); # $length not used
        foreach (split(/\s/, $anno_line)) {
            if (/^>(.+)$/) {
                $id = $1;
            } elsif (/g=(.*)$/) {
                $gene = $1;
            } elsif (/p=(.*)$/) {
                warn "\n### Warning:\nNo product annotation (p=) in file '$genome_file_path' on line:\n$anno_line\nProceeding ...\n" if (!$1);
                $product = $1;
                $product =~ s/\_/ /g;
            } elsif (/l=(.*)$/) { # $length not used
                my ($start, $stop) = split(/\.\./, $1);
                $length = abs($stop - $start) + 1;
            } elsif (/o=(.*)$/) {
                $organism = $1;
                $organism =~ s/\_/ /g;
            } elsif (/ec=(.*)$/) {
                $ec = $1;
            } else {
                die "\n### Fatal error:\nAnnotation in file '$genome_file_path' is not recognized on the following line. Please use 'cds_extractor.pl' to create your multi-FASTA protein genome files.\n$anno_line\n";
            }
        }
        die "\n### Fatal error:\nThe following FASTA ID of file '$genome_file_path' is not unique but has to be considering ALL input genome files (actually Proteinortho should have complained already). Please modify all repetitive occurences.\n$id\n" if ($Annotation{$id});
        #die "\n### Fatal error:\nThe FASTA ID '$id' is present in the multi-FASTA genome file '$genome_file_path' but not for this genome in the Proteinortho matrix '$PO_Matrix_File'. However, all IDs from the genome files in directory '$Genome_Dir' need to be present in the matrix and vice versa. Please run Proteinortho5 with option '-singles' to include also genes without orthologs, so-called singletons/ORFans (recommended is also option '-selfblast' to enhance paralog detection).\n" unless (check_po_id($id, $genome)); # subroutine to check if ID is present in %Ortho_Groups (from PO matrix); the reversed check (an ID that's present in PO matrix but not in the corresponding genome file) is in sub 'print_og_results' # this code was 1.) way too slow and 2.) somehow broke '$Query_Specific_OGs' and '$Query_Singletons' during implementation in 'po2anno.pl', no idea why ... (thus handle with care here)

        # fill annotation hash of hash
        $Annotation{$id} = {'genome' => $genome,
                            'gene' => $gene,
                            'product' => $product,
                            'organism' => $organism,
                            'ec' => $ec};
    }
    close $genome_fh;
}



### Check if CDS counts in multi-FASTA genome files and PO matrix are equal
print STDERR "Comparing annotation CDS counts to orthologous group CDS counts (test if everything is present) ...\n";
foreach my $genome (@Genome_Files) {
    my $ortho_cds_count = map ($Ortho_Groups{$_}->{$genome} ? @{ $Ortho_Groups{$_}->{$genome} } : (), keys %Ortho_Groups); # de-reference anonymous array in hash of hash
    # map evaluates BLOCK or EXPR in list context and returns the LIST value composed of the results of each such evaluation (each element of LIST may produce zero, one, or more elements in the returned value <=> grep returns only one value). In scalar context, returns the total number of elements so generated.
    # 'grep($Ortho_Groups{$_}->{$category} , keys %Ortho_Groups)' should also work (see result matrix)?

    # the map line above is a fancy way of writing
    #foreach (keys %Ortho_Groups) {
        #$ortho_cds_count += @{ $Ortho_Groups{$_}->{$genome} } if ($Ortho_Groups{$_}->{$genome});
    #}

    my $anno_cds_count = grep ($Annotation{$_}->{'genome'} eq $genome, keys %Annotation);
    die "\n### Fatal error:\nThere are more CDSs in file '$genome' than for this genome in the Proteinortho matrix '$PO_Matrix_File', but counts have to be equal. Please run Proteinortho5 with option '-singles' to include also genes without orthologs, so-called singletons/ORFans (recommended is also option '-selfblast' to enhance paralog detection).\n" if ($ortho_cds_count < $anno_cds_count); # overlaps with discarded (commented) sub 'check_po_id', but much more efficient (albeit not ID-specific)
    die "\n### Fatal error:\nThere are less CDSs in file '$genome' than for this genome in the Proteinortho matrix '$PO_Matrix_File', but counts have to be equal. Please check if the files are correct.\n" if ($ortho_cds_count > $anno_cds_count); # overlaps with ID-specific check in sub 'print_og_results', but here checks all IDs, the sub checks only those used for result file annotations
}



### Calculate stats from PO matrix according to the group cutoffs, categorize the results and print to result files
print STDERR "Calculating stats from the Proteinortho matrix according to the group cutoffs ...\n"; # run status of script

# open result file for binary matrix and print header
my $Binary_Matrix_File; # filename needed also for venn diagram R script below
my $Binary_Matrix_Fh;
if ($Opt_Binary) {
    $Binary_Matrix_File = "$Result_Dir/binary_matrix.tsv";
    open ($Binary_Matrix_Fh, ">", "$Binary_Matrix_File");
    print $Binary_Matrix_Fh join("\t", @Group_Names), "\n"; # print header for binary matrix; print array-elements tab-separated and after the last element a newline
}

my %Count_Group_OG; # hash of hash to store OG/CDS count and output filehandle for each category ('core', 'group-specific' etc.)
foreach my $og (sort { $a <=> $b } keys %Ortho_Groups) {
    my %present_og_genomes; # hash of array to store genomes of a group PRESENT in the current OG for cutoff tests and print outs

    # 'strict_core' category OG if present in ALL genomes (with option '-co')
    if ($Strict_Core && keys %{$Ortho_Groups{$og}} == @Genome_Files) {
        # always count OGs and CDSs for each category for final stats and plausability checks
        $Count_Group_OG{'strict_core'}->{'OG_count'}++;
        count_category_CDS(\@Genome_Files, 'strict_core', $og); # subroutine to count all CDS (also paralogs) in the current OG for ALL genomes in the indicated category (strict core category includes all genomes, thus @Genome_Files)
        open_result_file("$Result_Dir/strict_core_OGs.tsv", 'strict_core') if (!$Count_Group_OG{'strict_core'}->{'fh'}); # subroutine to open result file and store fh (if not opened yet)
        print_og_results($og, $Strict_Core, $Count_Group_OG{'strict_core'}->{'fh'}); # subroutine to print representative annotation for specified genome to result file (here genome is '$Strict_Core')
    }

    # count and save number of genomes in the current $og with their corresponding group affiliation
    foreach my $group (@Group_Names) { # keep initial group order in groups file (for representative annotation)
        foreach my $genome (@{ $Genome_Groups{$group}->{'genomes'} }) { # keep initial genome order for each group in groups file (for representative annotation)
            if ($Ortho_Groups{$og}->{$genome}) { # $genome present in the current $og
                push(@{ $present_og_genomes{$group} }, $genome); # push into anonymous array in hash; keeps initial genome order in groups

                # 'singleton/ORFan' category OG if only present in ONE genome (with option '-s')
                if ($Opt_Singletons && keys %{$Ortho_Groups{$og}} == 1) {
                    $Count_Group_OG{$genome}->{'OG_count'}++; # here: category = $genome
                    $Count_Group_OG{$genome}->{'CDS_count'} += @{ $Ortho_Groups{$og}->{$genome} }; # don't need sub 'count_category_CDS' because only one genome
                    if (!$Count_Group_OG{$genome}->{'fh'}) {
                        my $file_no_extension = $genome =~ s/\.\w+$//r; # remove file extensions from @Genome_Files to have nicer output files; non-destructive modifier causes the result of the substitution to be returned instead of modifying $_ (see http://perldoc.perl.org/perlrequick.html#Search-and-replace)
                        open_result_file("$Result_Dir/$file_no_extension\_singletons.tsv", $genome); # subroutine
                    }
                    print_og_results($og, $genome, $Count_Group_OG{$genome}->{'fh'}); # subroutine
                }
            }
        }
    }

    # filter genome counts of the groups in this OG according to inclusion-/exclusion-cutoffs
    my @OG_included_groups; # array for groups with number of genomes in this OG >= the inclusion cutoff
    my @OG_excluded_groups; # array for groups with genomes <= the exclusion cutoff
    my @OG_unspec_groups; # array for groups with 'exclusion < genome_count < inclusion'
    foreach my $group (@Group_Names) {
        if ($present_og_genomes{$group}) { # only run if group has a genome in this OG (otherwise undefined ARRAY reference)
            if (@{ $present_og_genomes{$group} } >= $Genome_Groups{$group}->{'inclusion'}) {
                push(@OG_included_groups, $group);
            } elsif (@{ $present_og_genomes{$group} } <= $Genome_Groups{$group}->{'exclusion'}) {
                push(@OG_excluded_groups, $group);
            } elsif (@{ $present_og_genomes{$group} } > $Genome_Groups{$group}->{'exclusion'} && @{ $present_og_genomes{$group} } < $Genome_Groups{$group}->{'inclusion'}) {
                push(@OG_unspec_groups, $group);
            }
        } else { # $present_og_genomes{$group} == 0 for group (undef anonymous array), i.e. always <= exclusion
            push(@OG_excluded_groups, $group);
        }
    }

    # 'unspecific' category OG if at least ONE group present in '@OG_unspec_groups'
    # no further categorization, i.e. skip rest of logics
    if (@OG_unspec_groups) {
        $Count_Group_OG{'unspecific'}->{'OG_count'}++;
        foreach my $group (@OG_unspec_groups) {
            count_category_CDS(\@{ $present_og_genomes{$group} }, 'unspecific', $og); # subroutine to count all CDS in the current OG for ALL genomes in the included groups (here unspec groups); %present_og_genomes contains only genomes present in this OG (see above 'push')
        }
        if ($Opt_Unspecific) { # print only to output file with option '-u'
            open_result_file("$Result_Dir/unspecific_OGs.tsv", 'unspecific') if (!$Count_Group_OG{'unspecific'}->{'fh'}); # subroutine
            print_og_results($og, @{ $present_og_genomes{$OG_unspec_groups[0]} }[0], $Count_Group_OG{'unspecific'}->{'fh'}); # subroutine; use first genome in first group which is present in current OG
        }
        next; # skip to next $og, because now unspec and not suitable for any further categories
    }

    # print binary matrix (in contrast to OGs in category 'unspecific' above, 'underrepresented' OGs below not skipped because package 'venn' can handle a row with only zeros)
    if ($Opt_Binary) {
        my $i = 0;
        foreach my $group (@Group_Names) {
            $i++;
            print $Binary_Matrix_Fh "1" if (grep {/$group/} @OG_included_groups);
            print $Binary_Matrix_Fh "0" if (grep {/$group/} @OG_excluded_groups);
            print $Binary_Matrix_Fh "\t" if ($i < @Group_Names);
        }
        print $Binary_Matrix_Fh "\n";
    }

    # 'cutoff_core' category OG if ALL groups >= inclusion (includes 'strict_core' OGs where ALL genomes are present in the OG)
    # 2 groups A/B [2 inclus :0 exclus], 3 groups A/B/C [3:0], 4 groups A/B/C/D [4:0]
    if (@OG_included_groups == @Group_Names && @OG_excluded_groups == 0) {
        $Count_Group_OG{'cutoff_core'}->{'OG_count'}++;
        foreach my $group (@OG_included_groups) {
            count_category_CDS(\@{ $present_og_genomes{$group} }, 'cutoff_core', $og); # subroutine
        }
        open_result_file("$Result_Dir/cutoff_core_OGs.tsv", 'cutoff_core') if (!$Count_Group_OG{'cutoff_core'}->{'fh'}); # subroutine
        print_og_results($og, @{ $present_og_genomes{$OG_included_groups[0]} }[0], $Count_Group_OG{'cutoff_core'}->{'fh'}); # subroutine; again use first genome in first group which is present in current OG

    # 'underrepresented' category OG if ALL groups <= exclusion_cutoff (depending on the cutoffs and number of genomes in the groups will include most singletons)
    # 2 groups [0:2], 3 [0:3], 4 [0:4]
    } elsif (@OG_included_groups == 0 && @OG_excluded_groups == @Group_Names) {
        $Count_Group_OG{'underrepresented'}->{'OG_count'}++;
        open_result_file("$Result_Dir/underrepresented_OGs.tsv", 'underrepresented') if (!$Count_Group_OG{'underrepresented'}->{'fh'}); # subroutine

        my $genome_found; # indicates if a genome with annotation found
        foreach my $group (@OG_excluded_groups) { # have to go through all groups in '@OG_excluded_groups' in case groups have no genomes ('$present_og_genomes{$group} == 0' above) in this OG; at least one group will have a genome (otherwise wouldn't be an OG in the PO matrix)
            if ($present_og_genomes{$group}) { # group with genome in this OG
                count_category_CDS(\@{ $present_og_genomes{$group} }, 'underrepresented', $og); # subroutine
                next if ($genome_found); # print annotation to output only for first genome found
                $genome_found = print_og_results($og, @{ $present_og_genomes{$group} }[0], $Count_Group_OG{'underrepresented'}->{'fh'}); # subroutine returns 1
            }
        }

    # 'group-specific' category OG if only ONE group >= inclusion_cutoff
    # 2 groups [1:1], 3 [1:2], 4 [1:3]
    } elsif (@OG_included_groups == 1 && @OG_excluded_groups == @Group_Names - 1) {
        $Count_Group_OG{"$OG_included_groups[0]_specific"}->{'OG_count'}++; # only one element in array thus [0]
        foreach my $group (@OG_included_groups) {
            count_category_CDS(\@{ $present_og_genomes{$group} }, "$OG_included_groups[0]_specific", $og); # subroutine
        }
        open_result_file("$Result_Dir/$OG_included_groups[0]_specific_OGs.tsv", "$OG_included_groups[0]_specific") if (!$Count_Group_OG{"$OG_included_groups[0]_specific"}->{'fh'}); # subroutine
        print_og_results($og, @{ $present_og_genomes{$OG_included_groups[0]} }[0], $Count_Group_OG{"$OG_included_groups[0]_specific"}->{'fh'}); # subroutine

    # 'group-absent' category OG if ONE group <= exclusion_cutoff (only for 3 and 4 groups)
    # 3 groups (A+B not C [2:1], A+C not B, and B+C not A) and 4 groups (A+B+C not D [3:1], A+B+D not C, A+C+D not B, B+C+D not A)
    } elsif (@Group_Names >= 3 && (@OG_included_groups == @Group_Names - 1 && @OG_excluded_groups == 1)) {
        $Count_Group_OG{"$OG_excluded_groups[0]_absent"}->{'OG_count'}++;
        foreach my $group (@OG_included_groups) {
            count_category_CDS(\@{ $present_og_genomes{$group} }, "$OG_excluded_groups[0]_absent", $og); # subroutine
        }
        open_result_file("$Result_Dir/$OG_excluded_groups[0]_absent_OGs.tsv", "$OG_excluded_groups[0]_absent") if (!$Count_Group_OG{"$OG_excluded_groups[0]_absent"}->{'fh'}); # subroutine
        print_og_results($og, @{ $present_og_genomes{$OG_included_groups[0]} }[0], $Count_Group_OG{"$OG_excluded_groups[0]_absent"}->{'fh'}); # subroutine

    # 'two group-specific' category OG if TWO groups >= inclusion_cutoff (only for 4 groups)
    # 4 groups (B+D not A or C [2:2], B+C not A or D, A+D not B or C, A+C not B or D, A+B not C or D, C+D not A or B)
    } elsif (@Group_Names == 4 && (@OG_included_groups == @OG_excluded_groups)) {
        $Count_Group_OG{"$OG_included_groups[0]-$OG_included_groups[1]_specific"}->{'OG_count'}++;
        foreach my $group (@OG_included_groups) {
            count_category_CDS(\@{ $present_og_genomes{$group} }, "$OG_included_groups[0]-$OG_included_groups[1]_specific", $og); # subroutine
        }
        open_result_file("$Result_Dir/$OG_included_groups[0]-$OG_included_groups[1]_specific_OGs.tsv", "$OG_included_groups[0]-$OG_included_groups[1]_specific") if (!$Count_Group_OG{"$OG_included_groups[0]-$OG_included_groups[1]_specific"}->{'fh'}); # subroutine
        print_og_results($og, @{ $present_og_genomes{$OG_included_groups[0]} }[0], $Count_Group_OG{"$OG_included_groups[0]-$OG_included_groups[1]_specific"}->{'fh'}); # subroutine

    } else {
        die "\n### Fatal error:\nNo category logic test fits for OG '$og', which shouldn't happen. Have you smuggled only one or five groups into the groups file '$Groups_File' (which also shouldn't be possible)? The script can handle only two, three, or four groups please edit the input!\n"; # shouldn't happen, since regex in group file parse checks for minimally two and maximally four tab-separated group names and all category logics should be covered
    }
}



### Close all output files
close $Binary_Matrix_Fh if ($Opt_Binary);
map { close $Count_Group_OG{$_}->{'fh'} } keys %Count_Group_OG;



### Plot the venn diagram with R package 'gplots' function 'venn'
if ($Opt_Venn) {
    my $tmp_r_script = "$Result_Dir/tmp_venn.R"; # filename of the R script
    my $venn_diagram_name = "$Result_Dir/venn_diagram.pdf"; # filename of the output PDF histogram
    print STDERR "Drawing venn diagram (option '-p') '$venn_diagram_name' ...\n"; # run status of script
    open (my $r_fh, ">", "$tmp_r_script");

    select $r_fh; # select fh for standard print/f output
    print "#!/usr/bin/Rscript --vanilla --slave\n"; # header of R script
    print "suppressMessages(library(gplots))\n"; # load library 'gplots' for function 'venn' and suppress STDOUT message from loading
    print "binary_matrix <- read.table(\"$Binary_Matrix_File\", header=TRUE, sep=\"\\t\")\n";
    print "pdf(\"$venn_diagram_name\")\n";
    print "venn(binary_matrix)\n"; # create the venn diagram
    print "out <- dev.off()\n"; # write histogram to PDF and suppress STDOUT output by diverting it
    select STDOUT; # change standard print/f output back to STDOUT
    close $r_fh;

    # execute tmp R script with 'Rscript'
    system("Rscript $tmp_r_script") == 0 or die "### Fatal error:\nPlotting the venn diagram didn't work. Either statistical programming language 'R' or the required R package 'gplots' is not installed, not in global \$PATH, or something wrong with the binary matrix '$Binary_Matrix_File' or the temporary R script '$tmp_r_script', which runs the drawing of the venn diagram. Function 'venn' of the 'gplots' package is required and needs to be installed!\n";
    unlink $tmp_r_script;
}



### Count overall category OGs and CDS, and plausibility check for above logics
my ($Total_Category_CDS, $Total_Category_OGs) = (0, 0);
foreach my $category (keys %Count_Group_OG) { # can't use a map in case of strict_core, singletons (categories which overlap non-optional categories)
    next if ($category =~ /strict_core/ || grep {$category eq $_} @Genome_Files); # skip optional 'strict_core' (option '-co') or singleton ('-s') categories

    $Total_Category_OGs += $Count_Group_OG{$category}->{'OG_count'};
    $Total_Category_CDS += $Count_Group_OG{$category}->{'CDS_count'};
}
die "\n### Fatal error:\nResulting categories don't include all present OGs, this shouldn't happen!\n" if (keys %Ortho_Groups != $Total_Category_OGs); # is this really always the case?



### Print final output stats matrix
print STDERR "Finished printing output files to result directory '$Result_Dir', printing final stats matrix to STDOUT ...\n"; # run status of script
print "# category\tOG_count\tCDS_count (all included genomes)\n"; # header for matrix
print "pangenome\t", scalar keys %Ortho_Groups, "\t", scalar keys %Annotation, "\n";
print "included_in_categories\t$Total_Category_OGs\t$Total_Category_CDS\n";
foreach my $category (sort {lc($a) cmp lc($b)} keys %Count_Group_OG) {
    print "$category";
    if (grep {$category eq $_} @Genome_Files) { # print overall stats for genomes with singletons (option '-s')
        print " genome\t", scalar grep($Ortho_Groups{$_}->{$category} , keys %Ortho_Groups), "\t", scalar grep($Annotation{$_}->{'genome'} eq $category, keys %Annotation), "\n"; # why does it work without the fancy map during CDS count check?
        print "$category singletons"; # include category "key word" singleton for category stats below
    }
    print "\t$Count_Group_OG{$category}->{'OG_count'}\t$Count_Group_OG{$category}->{'CDS_count'}\n";
}

# print out overall stats for all genomes even without singletons (option '-a')
if ($Opt_Genomes_Overall) {
    foreach my $genome (sort {lc($a) cmp lc($b)} @Genome_Files) {
        next if (grep {$genome eq $_} keys %Count_Group_OG); # skip genomes with singletons (option '-s'), already printed above
        print "$genome genome\t", scalar grep($Ortho_Groups{$_}->{$genome} , keys %Ortho_Groups), "\t", scalar grep($Annotation{$_}->{'genome'} eq $genome, keys %Annotation), "\n";
    }
}

exit;


#############
#Subroutines#
#############

### Subroutine to test for multi-FASTA genome file existence
sub check_file_exist {
    my $file = shift;
    die "\n### Fatal error:\nGenome file '$file' is listed in the Proteinortho matrix, but does not exist in directory '$Genome_Dir': $!\n" if (!-e $file);
    return 1;
}



### Subroutine to check if a string adheres to file system conformity, esp. for filenames
sub check_file_system_conformity {
    my ($string, $type) = @_;
    die "\n\n### Fatal error:\nOnly alphanumeric '0-9a-zA-Z', underscore '_', dash '-', and period '.' characters allowed for $type to avoid operating/file system problems. Thus, please adapt '$string' in the group file '$Groups_File' accordingly!\n" if ($string !~ /^[\w_\-\.]+$/);
    return 1;
}



### Subroutine to check if a genome filename is present in hash %Genome_Groups (hash filled from the $Groups_File)
sub check_genome_in_groups {
    my $genome = shift;
    foreach my $group (keys %Genome_Groups) {
        return 1 if (grep {$_ eq $genome} @{ $Genome_Groups{$group}->{'genomes'} });
    }
    return; # return 'undef'
}



### Subroutine to check if an ID is present in hash %Ortho_Group (hash filled from the PO matrix)
# this code was 1.) way too slow and 2.) somehow broke '$Query_Specific_OGs' and '$Query_Singletons' during implementation in 'po2anno.pl', no idea why ... (thus handle with care here)
#sub check_po_id {
    #my ($id, $genome) = @_;
    #foreach my $og (keys %Ortho_Groups) {
        #return 1 if (grep {$_ eq $id} @{ $Ortho_Groups{$og}->{$genome} });
    #}
    #return;
#}



### Subroutine to count and sum all CDS for ALL genomes in a certain category for a specific OG
sub count_category_CDS {
    my ($genomes_array_ref, $category, $og) = @_;
    foreach my $genome (@{ $genomes_array_ref }) { # de-reference array
        $Count_Group_OG{$category}->{'CDS_count'} += @{ $Ortho_Groups{$og}->{$genome} };
    }
    return 1;
}



### Subroutine to empty a directory with user interaction
sub empty_dir {
    my $dir = shift;
    print STDERR "\nDirectory '$dir' already exists! You can use either option '-r' to set a different output result directory name, or do you want to replace the directory and all its contents [y|n]? ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
        unlink glob "$dir/*"; # remove all files in results directory
    } else {
        die "\nScript abborted!\n";
    }
    return 1;
}



### Open a result file and print the header
sub open_result_file {
    my ($filename, $category) = @_;
    open ($Count_Group_OG{$category}->{'fh'}, ">", "$filename"); # store filehandle in hash of hash
    print { $Count_Group_OG{$category}->{'fh'} } "# OG#\tID\tgene\tEC_number\tproduct\torganism\n"; # block around fh needed if not a simple string variable (or glob)
    return 1;
}



### Print annotation of the OGs to the result files
sub print_og_results {
    my ($og, $genome, $fh) = @_;

    print $fh "$og"; # print OG number only for the first paralog
    foreach my $id (@{ $Ortho_Groups{$og}->{$genome} }) { # de-reference anonymous array for putative paralogs
        die "\n### Fatal error:\nID '$id' present in the Proteinortho matrix '$PO_Matrix_File' in OG '$og' but not present in the corresponding genome '$genome' multi-FASTA file. However, all IDs in the matrix and the genome files in directory '$Genome_Dir' need to be correspondent to each other. Make sure you chose the correct directory with the input genome files for the current Proteinortho matrix!\n" unless ($Annotation{$id}); # overlaps with CDS count check in multi-FASTA and PO matrix files, but here only IDs checked used for print out # the reverse check (an ID that's present in a genome file but not in the PO matrix) WAS during the parsing of the annotation from the multi-FASTA files and subroutine 'check_po_id($id)' (see sub why not)

        print $fh "\t$id\t";
        if ($Annotation{$id}->{'gene'}) {
            print $fh $Annotation{$id}->{'gene'}, "\t";
        } else {
            print $fh "\t";
        }
        if ($Annotation{$id}->{'ec'}) {
            print $fh $Annotation{$id}->{'ec'}, "\t";
        } else {
            print $fh "\t";
        }
        if ($Annotation{$id}->{'product'}) {
            print $fh $Annotation{$id}->{'product'}, "\t";
        } else { # very unlikely that a CDS has no product annotation
            print $fh "\t";
        }
        if ($Annotation{$id}->{'organism'}) {
            print $fh $Annotation{$id}->{'organism'}, "\n";
        } else {
            print $fh "\n";
        }
    }
    return 1;
}
