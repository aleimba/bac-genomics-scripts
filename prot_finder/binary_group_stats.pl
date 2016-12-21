#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<binary_group_stats.pl> - categorize binary matrix rows according to
column groups

=head1 SYNOPSIS

C<perl binary_group_stats.pl -i binary_matrix.tsv -g group_file.tsv
-p E<gt> overall_stats.tsv>

=head1 DESCRIPTION

Categorize rows of a delimited TEXT input B<binary> matrix (option
B<-i>) according to column group affiliations. All fields of the binary
matrix need to be filled with either a B<0> indicating absence or a
B<1> indicating presence, i.e. all rows need to have the same number of
columns. Use option B<-d> to set the delimiter of the input matrix,
default is set to tab-delimited/separated matrices.

The group affiliations of the columns are intended to get overall
presence/absence statistics for groups of columns and not simply
single columns of the matrix. Percentage inclusion (option
B<-cut_i>) and exclusion (option B<-cut_e>) cutoffs can be set to
define how strict the presence/absence of column groups within a row
are defined. Of course groups can also hold only single column
headers to get single column statistics. Group affiliations are
defined in a mandatory B<tab-delimited> group input file (option
B<-g>), including the column headers from the input binary matrix,
with B<minimal two> and B<maximal four> groups.

The script was designed to handle a presence/absence query protein
binary matrix from a C<prot_finder.pl> protein homolog search in
bacterial genomes with a subsequent C<prot_binary_matrix.pl> run to
create the matrix. But, any binary matrix can be used for
C<binary_group_stats.pl>, optionally beforehand transposed with
C<transpose_matrix.pl>. However, column headers in the first row and
row headers in the first column are B<mandatory> for the input
binary matrix. Only alphanumeric (a-z, A-Z, 0-9), underscore (_),
dash (-), and period (.) characters are allowed for the B<column
headers> and B<group names> in the group file (option B<-g>) to
avoid downstream problems with the operating/file system. As a
consequence, also no whitespaces are allowed in these! Additionally,
B<column headers>, B<row headers>, and B<group names> need to be
B<unique>.

<binary_group_stats.pl> is based upon
L<C<po2group_stats.pl>|/po2group_stats>, which does the same thing
for genomes in an ortholog/paralog output matrix from a
L<B<Proteinortho5>|http://www.bioinf.uni-leipzig.de/Software/proteinortho/>
calculation.

To explain the logic behind the categorization, the following
annotation for example groups will be used. A following '1'
exemplifies a group column presence count in a respective row E<gt>=
the rounded inclusion cutoff, a '0' a group column presence count
E<lt>= the rounded exclusion cutoff. The presence and absence of
rows for the column group affiliations are structured in different
categories depending on the number of groups. For B<two groups>
(e.g. A and B) there are five categories: 'A specific' (A:B = 1:0),
'B specific' (0:1), 'cutoff core' (1:1), 'underrepresented' (0:0),
and 'unspecific'. Unspecific rows have a column presence count for
at least B<one> group outside the cutoffs (exclusion cutoff E<lt>
column presence count E<lt> inclusion cutoff) and thus cannot be
categorized. The respective row headers of these 'unspecific' rows
will only be printed to a final result file with option B<-u>.
Overall stats for all categories are printed to C<STDOUT> in a final
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

The resulting B<group> presence/absence (according to the cutoffs)
can also be printed to a group binary matrix (option B<-b>) in the
result directory (option B<-r>), excluding the 'unspecific'
category. Since the categories are the logics underlying venn
diagrams, you can also plot the results in a venn diagram using the
group binary matrix (option B<-p>). The 'underrepresented' category
is exempt from the venn diagram, because it is outside of venn
diagram logics.

For an illustration of the logics have a look at the example venn
diagrams in the L<C<po2group_stats.pl>|/po2group_stats> folder
F</po2group_stats/pics/venn_diagram_logics.[svg|png]>.

There are two optional categories (which are only considered for the
final print outs and in the final stats matrix, not for the group
binary matrix and the venn diagram): 'strict core' (option B<-co>)
for rows where B<all> columns have a presence '1'. Of course all the
'strict core' rows are also included in the 'cutoff core' category
('strict core' is identical to 'cutoff core' with B<-cut_i> 1 and
B<-cut_e> 0). Option B<-s> activates the detection of
'singleton/column-specific' rows with a '1' in only B<one> column.
Depending on the cutoffs and number of columns in the groups,
category 'underrepresented' includes most of these singletons.

Each row header is printed in the respective category output file in
the result directory. The number of row headers in the category
result files are the same as listed in the venn diagram and the
final stats matrix. Groups with only a single column header will
have the same result as the respective 'singleton' category (with
option B<-s>).

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input delimited TEXT binary matrix (e.g. *.tsv, *.csv, or *.txt),
see option B<-d>

=item B<-g>=I<str>, B<-groups_file>=I<str>

Tab-delimited file with group affiliation for the columns from the
input binary matrix with B<minimal two> and B<maximal four> groups
(easiest to create in a spreadsheet software and save in
tab-separated format). B<All> column headers from the input binary
matrix need to be included. Column headers and group names can only
include alphanumeric (a-z, A-Z, 0-9), underscore (_), dash (-), and
period (.) characters (no whitespaces allowed either). Example
format with two column headers in group A, three in group B and D,
and one in group C:

group_A\tgroup_B\tgroup_C\tgroup_D
column_header1\tcolumn_header9\tcolumn_header3\tcolumn_header8
column_header7\tcolumn_header6\t\tcolumn_header5
\tcolumn_header4\t\tcolumn_header2

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-d>=I<str>, B<-delimiter>=I<str>

Set delimiter of input binary matrix (e.g. comma ',', single
space ' ' etc.) [default = tab-delimited/separated]

=item B<-r>=I<str>, B<-result_dir>=I<str>

Path to result folder [default = inclusion and exclusion percentage
cutoff, './results_i#_e#']

=item B<-cut_i>=I<float>, B<-cut_inclusion>=I<float>

Percentage inclusion cutoff for column presence counts in a group
per row, has to be E<gt> 0 and E<lt>= 1. Cutoff will be rounded
according to the column header number in each group and has to be
E<gt> the rounded exclusion cutoff in this group. [default = 0.9]

=item B<-cut_e>=I<float>, B<-cut_exclusion>=I<float>

Percentage exclusion cutoff, has to be E<gt>= 0 and E<lt> 1. Rounded
cutoff has to be E<lt> rounded inclusion cutoff. [default = 0.1]

=item B<-b>, B<-binary_group_matrix>

Print a group binary matrix with the presence/absence column group
results according to the cutoffs (excluding 'unspecific' category
rows)

=item B<-p>, B<-plot_venn>

Plot venn diagram from the group binary matrix (except 'unspecific'
and 'underrepresented' categories) with function C<venn> from B<R>
package B<gplots>, requires option B<-b>

=item B<-co>, B<-core_strict>

Include 'strict core' category in output for rows where B<all>
columns have a '1'

=item B<-s>, B<-singletons>

Include singleton/column-specific rows for each column header in the
output, activates also overall column header presence ('1') counts
in final stats matrix for columns with singletons

=item B<-u>, B<-unspecific>

Include 'unspecific' category in output

=item B<-a>, B<-all_column_presence_overall>

Report overall presence counts for all column headers (appended to
the final stats matrix), also those without singletons; will include
all overall column header presence counts without option B<-s>

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

=item F<./results_i#_e#/[*_specific|*_absent|cutoff_core|underrepresented]_rows.txt>

Files including the row headers for rows in non-optional categories

=item (F<./results_i#_e#/[*_singletons|strict_core|unspecific]_rows.txt>)

Optional category output files with the respective row headers

=item (F<./results_i#_e#/binary_matrix.tsv>)

Tab-delimited binary matrix of group presence/absence results
according to cutoffs (excluding 'unspecific' rows)

=item (F<./results_i#_e#/venn_diagram.pdf>)

Venn diagram for non-optional categories (except 'unspecific' and
'underrepresented' categories)

=back

=head1 EXAMPLES

=over

=item C<perl binary_group_stats.pl -i binary_matrix_transposed.csv -g group_file.tsv -d , -r result_dir -cut_i 0.7 -cut_e 0.2 -b -p -co -s -u -a E<gt> overall_stats.tsv>

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

 0.1                                                       06-06-2016

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
my $Binary_Matrix_File; # binary input matrix
my $Groups_File; # tab-separated file with group affiliation for the input binary matrix column headers
my $Delimiter = "\t"; # set delimiter/separator of input matrix
my $Result_Dir; # path to results folder; default is set below to '"./results_i".$Inclusion_Cutoff."c".$Exclusion_Cutoff'
my $Inclusion_Cutoff = 0.9; # inclusion percentage cutoff for column presence ('1') count per row (floating point number >0 && <=1)
my $Exclusion_Cutoff = 0.1; # exclusion percentage cutoff (floating point number >=0 && <1)
my $Opt_Group_Binary; # print resulting group binary matrix according to the group inclusion/exclusion results
my $Opt_Venn; # create venn diagram with R package 'gplots' function 'venn'; requires $Opt_Group_Binary
my $Opt_Strict_Core; # report also strict core, i.e. rows with presence in ALL columns of input binary matrix (identical to 'cutoff_core' with inclusion and exclusion cutoffs of 1 and 0 resp.)
my $Opt_Singletons; # report also singletons for each column header (outside of the group classifications)
my $Opt_Unspecific; # report also group-unspecific rows ('exclusion < row_group_column-presence_count < inclusion' for any group)
my $Opt_Columns_Overall; # report overall row presence counts also for column headers without singletons
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Binary_Matrix_File,
            'groups_file=s' => \$Groups_File,
            'delimiter=s' => \$Delimiter,
            'result_dir=s' => \$Result_Dir,
            'cut_inclusion=f' => \$Inclusion_Cutoff,
            'cut_exclusion=f' => \$Exclusion_Cutoff,
            'binary_group_matrix' => \$Opt_Group_Binary,
            'plot_venn' => \$Opt_Venn,
            'core_strict' => \$Opt_Strict_Core,
            'singletons' => \$Opt_Singletons,
            'unspecific' => \$Opt_Unspecific,
            'all_column_presence_overall' => \$Opt_Columns_Overall,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD and enforce mandatory options
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);

if (!$Binary_Matrix_File || !$Groups_File) {
    my $warning = "\n### Fatal error: Mandatory options '-i' or '-g' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}

if (($Inclusion_Cutoff <= 0 || $Inclusion_Cutoff > 1) || ($Exclusion_Cutoff < 0 || $Exclusion_Cutoff >= 1)) {
    my $warning = "\n### Fatal error:\nInclusion (0 < '-cut_i' <= 1) or exclusion (0 <= '-cut_e' < 1) cutoff(s) not chosen correctly!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}

print STDERR "Script call command: $Cmdline\n"; # print call command after '-h|-v'

if ($Opt_Venn && !$Opt_Group_Binary) {
    warn "Option '-p' to draw the venn diagram set, but not it's required option '-b' for the output group binary matrix. Forcing option '-b'!\n";
    $Opt_Group_Binary = 1;
}



### Parse column headers (from the input binary matrix) in groups file
print STDERR "Parsing group file '$Groups_File' "; # run status of script
open (my $Groups_File_Fh, "<", "$Groups_File");

my @Group_Names; # store group name column order
my %Column_Header_Groups; # hash of hash with group name as key, internal hash with resp. input binary matrix column header anonymous array (key 'column_headers'), or calculated integer inclusion/exclusion cutoffs (key 'inclusion' etc.)
while (<$Groups_File_Fh>) {
    chomp;
    next if (/^\s*$/); # skip empty lines

    # check groups file header and get group names
    if ($. == 1) { # header of groups file (first line)
        die "\n### Fatal error:\nGroups file '$Groups_File' does not have a correctly formatted header line with at least TWO and maximal FOUR tab-separated group name columns (without other whitespaces)!\n" if (!/^\S+\t\S+\t?\S*\t?\S*$/);
        foreach my $group_name (split(/\t/)) {
            check_file_system_conformity($group_name, 'group names'); # subroutine to check string only contains alphanumeric '0-9a-zA-Z', underscore '_', dash '-', and period '.' characters
            die "\n### Fatal error:\nGroup name '$group_name' is not unique in the group file '$Groups_File'! Please choose unique names for the groups in the group file!\n" if (grep {/$group_name/} @Group_Names);
            push (@Group_Names, $group_name); # store groups array
        }
        next; # skip header line for subsequent code
    }

    # parse input binary matrix column headers in groups file
    die "\n### Fatal error:\nLine '$.' in the groups file '$Groups_File' does not have the correct format with maximal four tab-separated column header names from the input binary matrix '$Binary_Matrix_File' according to the groups (without other whitespaces):\n$_\n" if (!/^\S*\t?\S*\t?\S*\t?\S*$/);
    my $column = -1; # group file column counter to associate input binary matrix column header to correct group (@Group_Names array zero-based, thus start with '-1')
    foreach my $column_header (split(/\t/)) {
        $column++;
        next if ($column_header =~ /^$/); # skip empty columns in group file, in case the groups have uneven column header numbers
        warn "\n### Warning:\nInput binary matrix column header '$column_header' in the group file has a binary ('0' or '1') instead of an alphanumeric string, sure this is a column header? The matrix needs to have mandatory row and column headers, see help with '-h'!\n" if ($column_header =~ /^(0|1)$/);
        check_file_system_conformity($column_header, 'column headers'); # subroutine
        die "\n### Fatal error:\nInput binary matrix column header '$column_header' is not unique in the group file '$Groups_File', but needs to be in the group file and the input binary matrix '$Binary_Matrix_File'! The group file needs to contain all of the column headers from the input binary matrix and vice versa. Make sure the strings correspond exactly (also lower/uppercase) to the header of the matrix!\n\n" if (check_column_header_in_groups($column_header)); # subroutine to check if a column header is present in %Column_Header_Groups
        push (@{ $Column_Header_Groups{$Group_Names[$column]}->{'column_headers'} }, $column_header); # push into anonymous array in hash of hash
    }
}
close $Groups_File_Fh;
print STDERR "with ", scalar map(@{ $Column_Header_Groups{$_}->{'column_headers'} }, keys %Column_Header_Groups), " input binary matrix column headers ...\n"; # run status of script

foreach (@Group_Names) {
    die "\n### Fatal error:\nGroup '$_' does not contain an element, please add a suitable input binary matrix column header or remove the group from the group file '$Groups_File'!\n" if (!$Column_Header_Groups{$_});
}



### Round inclusion and exclusion count cutoffs for each group
foreach (sort keys %Column_Header_Groups) {
    $Column_Header_Groups{$_}->{'inclusion'} = int((scalar @{ $Column_Header_Groups{$_}->{'column_headers'} } * $Inclusion_Cutoff) + 0.5); # round positive number half up to integer ((s)printf rounds half to even), see e.g. https://stackoverflow.com/questions/178539/how-do-you-round-a-floating-point-number-in-perl
    print STDERR "Group '$_' with ", scalar @{ $Column_Header_Groups{$_}->{'column_headers'} }, " input binary matrix column headers has a rounded inclusion cutoff of $Column_Header_Groups{$_}->{'inclusion'} column header(s) and ";

    $Column_Header_Groups{$_}->{'exclusion'} = int((scalar @{ $Column_Header_Groups{$_}->{'column_headers'} } * $Exclusion_Cutoff) + 0.5);
    print STDERR "an exclusion cutoff of $Column_Header_Groups{$_}->{'exclusion'} column header(s)!\n";
    die "\n### Fatal error:\nThe rounded inclusion cutoff has to be greater than the exclusion cutoff for all groups, please choose appropriate cutoffs!\n" if ($Column_Header_Groups{$_}->{'inclusion'} <= $Column_Header_Groups{$_}->{'exclusion'});
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



### Parse rows in input binary matrix
print STDERR "Parsing input binary matrix '$Binary_Matrix_File' ...\n"; # run status of script
open (my $Binary_Matrix_Fh, "<", "$Binary_Matrix_File");

my @Column_Headers; # store input binary matrix column header order
my %Row_Header_Groups; # hash of hash with row header as key (first column of each row) and internal hash with column header as key and '1' for presence ('0' absence is not stored)
while (<$Binary_Matrix_Fh>) {
    chomp;
    next if (/^\s*$/); # skip empty lines
    die "\n### Fatal error:\nSet delimiter/separator '$Delimiter' (option '-d') not found in line '$.' of the input binary matrix '$Binary_Matrix_File', sure the correct one is set?\n" if ($_ !~ /$Delimiter/);

    # check input binary matrix file header and get column header names
    if ($. == 1) { # header of input binary matrix file (first line)
        my (undef, @headers) = split(/$Delimiter/); # get rid of first line row header
        foreach my $column_header (@headers) {
            die "\n### Fatal error:\nColumn header '$column_header' is present in the input binary matrix '$Binary_Matrix_File' but not in the groups file '$Groups_File'. However, the group file needs to contain all of the column headers from the matrix and vice versa. Make sure the strings correspond exactly (also lower/uppercase) to the matrix header!\n" unless (check_column_header_in_groups($column_header)); # subroutine to check if a column header is present in %Column_Header_Groups
            push (@Column_Headers, $column_header);
        }
        die "\n### Fatal error:\nNot the same amount of input binary matrix column headers in the group file '$Groups_File' and matrix '$Binary_Matrix_File', but the group file needs to contain all of the column headers in the matrix and vice versa!\n" if (@Column_Headers != map(@{ $Column_Header_Groups{$_}->{'column_headers'} }, keys %Column_Header_Groups));
        next; # skip header line of input binary matrix for subsequent code
    }

    # parse input binary matrix
    my ($row_header, @presence_absence) = split(/$Delimiter/);
    die "\n### Fatal error:\nInput binary matrix '$Binary_Matrix_File' row '$row_header' doesn't have the same amount of columns as overall column headers, however the matrix needs to be filled/consistent in every row!\n" unless (@Column_Headers == @presence_absence);
    die "\n### Fatal error:\nInput binary matrix row header '$row_header' is not unique in the input binary matrix '$Binary_Matrix_File', however all row headers need to be unique!\n" if ($Row_Header_Groups{$row_header});
    warn "\n### Warning:\nLine number '$.' of the binary input matrix '$Binary_Matrix_File' has a binary ('0' or '1') row header instead of an alphanumeric string, sure this is a row header? The matrix needs to have mandatory row and column headers, see help with '-h'!\n\n" if ($row_header =~ /^(0|1)$/);
    my $column = -1; # column counter to associate presence/absence field to correct column header (@Column_Headers array zero-based, thus start with '-1')
    foreach (@presence_absence) {
        $column++;
        die "\n###Fatal error:\nThere's an empty field in the input binary matrix '$Binary_Matrix_File' in row '$row_header', column '$column+2'. This shouldn't happen, as all the fields of the BINARY matrix should either have a '0' indicating absence or a '1' indicating presence!\n" if (/^$/);
        die "\n### Fatal error:\nThe input binary matrix '$Binary_Matrix_File' contains a field with '$_' in row/line '$.', column '$column+2', instead of the mandatory binary '0' indicating absence or '1' indicating presence!\n" if (!/^(0|1)$/);
        next if (!$_); # skip absence columns in binary input matrix (indicated by '0'), have to skip after '$column++' otherwise column won't fit
        $Row_Header_Groups{$row_header}{$Column_Headers[$column]} = $_; # store presence ('1') in hash of hash
    }
}
close $Binary_Matrix_Fh;



### Calculate column stats from input binary matrix according to the group cutoffs, categorize the results and print to result files
print STDERR "Calculating column stats from the input binary matrix according to the group cutoffs ...\n"; # run status of script

# open result file for group binary matrix and print header
my $Group_Binary_Matrix_File; # filename needed also for venn diagram R script below
my $Group_Binary_Matrix_Fh;
if ($Opt_Group_Binary) {
    $Group_Binary_Matrix_File = "$Result_Dir/group_binary_matrix.tsv";
    open ($Group_Binary_Matrix_Fh, ">", "$Group_Binary_Matrix_File");
    print $Group_Binary_Matrix_Fh join("\t", @Group_Names), "\n"; # print header for group binary matrix
}

my %Count_Group_Pres_Abs; # hash of hash to store presence count and output filehandle for each category ('core', 'group-specific' etc.)
foreach my $row_header (sort {lc($a) cmp lc($b)} keys %Row_Header_Groups) {
    my %present_row_columns; # hash of array to store input binary matrix columns of a group with PRESENCE ('1') in the current row for cutoff tests and print outs

    # 'strict_core' category row, if ALL columns present in row (with option '-co')
    if ($Opt_Strict_Core && keys %{$Row_Header_Groups{$row_header}} == @Column_Headers) {
        $Count_Group_Pres_Abs{'strict_core'}->{'presence_count'}++;
        open_result_file("$Result_Dir/strict_core_rows.txt", 'strict_core') if (!$Count_Group_Pres_Abs{'strict_core'}->{'fh'}); # subroutine to open result file and store fh (if not opened yet)
        print { $Count_Group_Pres_Abs{'strict_core'}->{'fh'} } "$row_header\n"; # block around fh needed if not a simple string variable (or glob)
    }

    # count and save number of columns with presence ('1') for the current $row_header with their corresponding group affiliation
    foreach my $group (@Group_Names) {
        foreach my $column_header (@{ $Column_Header_Groups{$group}->{'column_headers'} }) {
            if ($Row_Header_Groups{$row_header}->{$column_header}) { # $column_header with presence in current $row_header
                push(@{ $present_row_columns{$group} }, $column_header); # push into anonymous array in hash

                # 'singleton' category row if only ONE column with presence (with option '-s')
                if ($Opt_Singletons && keys %{$Row_Header_Groups{$row_header}} == 1) {
                    $Count_Group_Pres_Abs{$column_header}->{'presence_count'}++; # here: category = $column_header
                    open_result_file("$Result_Dir/$column_header\_singletons.txt", $column_header) if (!$Count_Group_Pres_Abs{$column_header}->{'fh'}); # subroutine
                    print { $Count_Group_Pres_Abs{$column_header}->{'fh'} } "$row_header\n";
                }
            }
        }
    }

    # filter column presence counts of the groups in this row according to inclusion-/exclusion-cutoffs
    my @row_included_groups; # array for groups with number of present column in this row >= the inclusion cutoff
    my @row_excluded_groups; # array for groups with column presence <= the exclusion cutoff
    my @row_unspec_groups; # array for unspecific groups with 'exclusion < column_presence_count < inclusion'
    foreach my $group (@Group_Names) {
        if ($present_row_columns{$group}) { # only run if group has a column presence in this row (otherwise undefined ARRAY reference)
            if (@{ $present_row_columns{$group} } >= $Column_Header_Groups{$group}->{'inclusion'}) {
                push(@row_included_groups, $group);
            } elsif (@{ $present_row_columns{$group} } <= $Column_Header_Groups{$group}->{'exclusion'}) {
                push(@row_excluded_groups, $group);
            } elsif (@{ $present_row_columns{$group} } > $Column_Header_Groups{$group}->{'exclusion'} && @{ $present_row_columns{$group} } < $Column_Header_Groups{$group}->{'inclusion'}) {
                push(@row_unspec_groups, $group);
            }
        } else { # $present_row_columns{$group} == 0 for group (undef anonymous array), i.e. always <= exclusion
            push(@row_excluded_groups, $group);
        }
    }

    # 'unspecific' category row if at least ONE group present in '@row_unspec_groups'
    # no further categorization, i.e. skip rest of logics
    if (@row_unspec_groups) {
        $Count_Group_Pres_Abs{'unspecific'}->{'presence_count'}++;
        if ($Opt_Unspecific) { # print only to output file with option '-u'
            open_result_file("$Result_Dir/unspecific_rows.txt", 'unspecific') if (!$Count_Group_Pres_Abs{'unspecific'}->{'fh'}); # subroutine
            print { $Count_Group_Pres_Abs{'unspecific'}->{'fh'} } "$row_header\n";
        }
        next; # skip to next $row_header, because now unspec and not suitable for any further categories
    }

    # print group binary matrix (in contrast to rows in category 'unspecific' above, 'underrepresented' rows below not skipped because package 'venn' can handle a row with only zeros in the group binary matrix)
    if ($Opt_Group_Binary) {
        my $i = 0;
        foreach my $group (@Group_Names) {
            $i++;
            print $Group_Binary_Matrix_Fh "1" if (grep {/$group/} @row_included_groups);
            print $Group_Binary_Matrix_Fh "0" if (grep {/$group/} @row_excluded_groups);
            print $Group_Binary_Matrix_Fh "\t" if ($i < @Group_Names);
        }
        print $Group_Binary_Matrix_Fh "\n";
    }

    # 'cutoff_core' category row if ALL groups >= inclusion (includes 'strict_core' rows where ALL columns are present in the row)
    # 2 groups A/B [2 inclus :0 exclus], 3 groups A/B/C [3:0], 4 groups A/B/C/D [4:0]
    if (@row_included_groups == @Group_Names && @row_excluded_groups == 0) {
        $Count_Group_Pres_Abs{'cutoff_core'}->{'presence_count'}++;
        open_result_file("$Result_Dir/cutoff_core_rows.txt", 'cutoff_core') if (!$Count_Group_Pres_Abs{'cutoff_core'}->{'fh'}); # subroutine
        print { $Count_Group_Pres_Abs{'cutoff_core'}->{'fh'} } "$row_header\n";

    # 'underrepresented' category row if ALL groups <= exclusion_cutoff (depending on the cutoffs and number of columns in the groups will include most singletons)
    # 2 groups [0:2], 3 [0:3], 4 [0:4]
    } elsif (@row_included_groups == 0 && @row_excluded_groups == @Group_Names) {
        $Count_Group_Pres_Abs{'underrepresented'}->{'presence_count'}++;
        open_result_file("$Result_Dir/underrepresented_rows.txt", 'underrepresented') if (!$Count_Group_Pres_Abs{'underrepresented'}->{'fh'}); # subroutine
        print { $Count_Group_Pres_Abs{'underrepresented'}->{'fh'} } "$row_header\n";

    # 'group-specific' category row if only ONE group >= inclusion_cutoff
    # 2 groups [1:1], 3 [1:2], 4 [1:3]
    } elsif (@row_included_groups == 1 && @row_excluded_groups == @Group_Names - 1) {
        $Count_Group_Pres_Abs{"$row_included_groups[0]_specific"}->{'presence_count'}++; # only one element in array thus [0]
        open_result_file("$Result_Dir/$row_included_groups[0]_specific_rows.txt", "$row_included_groups[0]_specific") if (!$Count_Group_Pres_Abs{"$row_included_groups[0]_specific"}->{'fh'}); # subroutine
        print { $Count_Group_Pres_Abs{"$row_included_groups[0]_specific"}->{'fh'} } "$row_header\n";

    # 'group-absent' category row if ONE group <= exclusion_cutoff (only for 3 and 4 groups)
    # 3 groups (A+B not C [2:1], A+C not B, and B+C not A) and 4 groups (A+B+C not D [3:1], A+B+D not C, A+C+D not B, B+C+D not A)
    } elsif (@Group_Names >= 3 && (@row_included_groups == @Group_Names - 1 && @row_excluded_groups == 1)) {
        $Count_Group_Pres_Abs{"$row_excluded_groups[0]_absent"}->{'presence_count'}++;
        open_result_file("$Result_Dir/$row_excluded_groups[0]_absent_rows.txt", "$row_excluded_groups[0]_absent") if (!$Count_Group_Pres_Abs{"$row_excluded_groups[0]_absent"}->{'fh'}); # subroutine
        print { $Count_Group_Pres_Abs{"$row_excluded_groups[0]_absent"}->{'fh'} } "$row_header\n";

    # 'two group-specific' category row if TWO groups >= inclusion_cutoff (only for 4 groups)
    # 4 groups (B+D not A or C [2:2], B+C not A or D, A+D not B or C, A+C not B or D, A+B not C or D, C+D not A or B)
    } elsif (@Group_Names == 4 && (@row_included_groups == @row_excluded_groups)) {
        $Count_Group_Pres_Abs{"$row_included_groups[0]-$row_included_groups[1]_specific"}->{'presence_count'}++;
        open_result_file("$Result_Dir/$row_included_groups[0]-$row_included_groups[1]_specific_rows.txt", "$row_included_groups[0]-$row_included_groups[1]_specific") if (!$Count_Group_Pres_Abs{"$row_included_groups[0]-$row_included_groups[1]_specific"}->{'fh'}); # subroutine
        print { $Count_Group_Pres_Abs{"$row_included_groups[0]-$row_included_groups[1]_specific"}->{'fh'} } "$row_header\n";

    } else {
        die "\n### Fatal error:\nNo category logic test fits for row '$row_header', which shouldn't happen. Have you smuggled only one or five groups into the groups file '$Groups_File' (which also shouldn't be possible)? The script can handle only two, three, or four groups please edit the input!\n"; # shouldn't happen, since regex in group file parse checks for minimally two and maximally four tab-separated group names and all category logics should be covered
    }
}



### Close all output files
close $Group_Binary_Matrix_Fh if ($Opt_Group_Binary);
map { close $Count_Group_Pres_Abs{$_}->{'fh'} } keys %Count_Group_Pres_Abs;



### Plot the venn diagram with R package 'gplots' function 'venn'
if ($Opt_Venn) {
    my $tmp_r_script = "$Result_Dir/tmp_venn.R"; # filename of the R script
    my $venn_diagram_name = "$Result_Dir/venn_diagram.pdf"; # filename of the output PDF histogram
    print STDERR "Drawing venn diagram (option '-p') '$venn_diagram_name' ...\n"; # run status of script
    open (my $r_fh, ">", "$tmp_r_script");

    select $r_fh; # select fh for standard print/f output
    print "#!/usr/bin/Rscript --vanilla --slave\n"; # header of R script
    print "suppressMessages(library(gplots))\n"; # load library 'gplots' for function 'venn' and suppress STDOUT message from loading
    print "group_binary_matrix <- read.table(\"$Group_Binary_Matrix_File\", header=TRUE, sep=\"\\t\")\n";
    print "pdf(\"$venn_diagram_name\")\n";
    print "venn(group_binary_matrix)\n"; # create the venn diagram
    print "out <- dev.off()\n"; # write histogram to PDF and suppress STDOUT output by diverting it
    select STDOUT; # change standard print/f output back to STDOUT
    close $r_fh;

    # execute tmp R script with 'Rscript'
    system("Rscript $tmp_r_script") == 0 or die "### Fatal error:\nPlotting the venn diagram didn't work. Either statistical programming language 'R' or the required R package 'gplots' is not installed, not in global \$PATH, or something wrong with the group binary matrix '$Group_Binary_Matrix_File' or the temporary R script '$tmp_r_script', which runs the drawing of the venn diagram. Function 'venn' of the 'gplots' package is required and needs to be installed!\n";
    unlink $tmp_r_script;
}



### Count overall category rows, and plausibility check for above logics
my $Total_Category_Rows = 0;
foreach my $category (keys %Count_Group_Pres_Abs) { # can't use a map in case of optionally strict_core, singletons categories (which overlap non-optional categories)
    next if ($category =~ /strict_core/ || grep {$category eq $_} @Column_Headers); # skip optional 'strict_core' (option '-co') or singleton ('-s') categories

    $Total_Category_Rows += $Count_Group_Pres_Abs{$category}->{'presence_count'};
}
die "\n### Fatal error:\nResulting categories don't include all present rows, this shouldn't happen!\n" if (keys %Row_Header_Groups != $Total_Category_Rows); # is this really always the case?



### Print final output stats matrix
print STDERR "Finished printing output files to result directory '$Result_Dir', printing final stats matrix to STDOUT ...\n"; # run status of script
print "# category\tpresence_count\n"; # header for matrix
print "total_rows\t", scalar keys %Row_Header_Groups, "\n";
foreach my $category (sort {lc($a) cmp lc($b)} keys %Count_Group_Pres_Abs) {
    print "$category";
    if (grep {$category eq $_} @Column_Headers) { # print overall presence counts for column headers with singletons (option '-s')
        print " overall_presence\t", scalar grep($Row_Header_Groups{$_}->{$category} , keys %Row_Header_Groups), "\n";
        print "$category singletons"; # include category "key word" singleton for category stats below
    }
    print "\t$Count_Group_Pres_Abs{$category}->{'presence_count'}\n";
}

# print out overall presence counts for all columns even without singletons (option '-a')
if ($Opt_Columns_Overall) {
    foreach my $column_header (sort {lc($a) cmp lc($b)} @Column_Headers) {
        next if (grep {$column_header eq $_} keys %Count_Group_Pres_Abs); # skip column headers with singletons (option '-s'), already printed above
        print "$column_header overall_presence\t", scalar grep($Row_Header_Groups{$_}->{$column_header} , keys %Row_Header_Groups), "\n";
    }
}

exit;


#############
#Subroutines#
#############

### Subroutine to check if a column header is present in hash %Column_Header_Groups (hash filled from the $Groups_File)
sub check_column_header_in_groups {
    my $column_header = shift;
    foreach my $group (keys %Column_Header_Groups) {
        return 1 if (grep {$_ eq $column_header} @{ $Column_Header_Groups{$group}->{'column_headers'} });
    }
    return; # return 'undef'
}



### Subroutine to check if a string adheres to file system conformity, esp. for filenames
sub check_file_system_conformity {
    my ($string, $type) = @_;
    die "\n\n### Fatal error:\nOnly alphanumeric '0-9a-zA-Z', underscore '_', dash '-', and period '.' characters allowed for $type to avoid operating/file system problems. Thus, please adapt '$string' in the group file '$Groups_File' and if necessary also in the input binary matrix '$Binary_Matrix_File' accordingly!\n" if ($string !~ /^[\w_\-\.]+$/);
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
    open ($Count_Group_Pres_Abs{$category}->{'fh'}, ">", "$filename"); # store filehandle in hash of hash
    print { $Count_Group_Pres_Abs{$category}->{'fh'} } "# row_header\n"; # block around fh needed if not a simple string variable (or glob)
    return 1;
}
