#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<transpose_matrix.pl> - transpose a delimited TEXT matrix

=head1 SYNOPSIS

C<perl transpose_matrix.pl input_matrix.tsv E<gt>
input_matrix_transposed.tsv>

B<or>

C<perl prot_binary_matrix.pl blast_hits.tsv | perl
transpose_matrix.pl E<gt> binary_matrix_transposed.tsv>

=head1 DESCRIPTION

This script transposes a delimited TEXT input matrix, i.e. rows will
become columns and columns rows. Use option B<-d> to set the
delimiter of the input and output matrix, default is set to
tab-delimited/separated matrices. Input matrices can be given
directly via C<STDIN> or as a file. The script is intended for the
resulting presence/absence binary matrices of
C<prot_binary_matrix.pl>, but can be used for any TEXT matrix.

The binary matrix of C<prot_binary_matrix.pl> has the query protein
IDs as column headers and the subject genomes as row headers. Thus,
C<transpose_matrix.pl> is very useful to transpose the
C<prot_binary_matrix.pl> matrix for the usage with
C<binary_group_stats.pl> to calculate presence/absence statistics
for groups of columns/genomes (and not simply single columns of the
matrix). C<binary_group_stats.pl> also has a comprehensive manual
with its option B<-h>.

Additionally, option B<-e> can be used to fill empty cells of the
input matrix with a value in the transposed matrix (e.g. 'NA', '0'
etc.).

=head1 OPTIONS

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-d>=I<str>, B<-delimiter>=I<str>

Set delimiter of input and output matrix (e.g. comma ',', single
space ' ' etc.) [default = tab-delimited/separated]

=item B<-e>=I<str>, B<-empty>=I<str>

Fill empty cells of the input matrix with a value in the transposed
matrix (e.g. 'NA', '0' etc.)

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The transposed matrix is printed to C<STDOUT>. Redirect or pipe into
another tool as needed.

=back

=head1 EXAMPLES

=over

=item C<perl transpose_matrix.pl -d ' ' -e NA input_matrix_space-delimit.txt E<gt> input_matrix_space-delimit_transposed.txt>

=back

B<or>

=over

=item C<for matrix in *.tsv; do perl transpose_matrix.pl "$matrix" E<gt> "${matrix%.*}_transposed.tsv"; done>

=back

B<or>

=over

=item C<perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl -l -c | perl transpose_matrix.pl -d , E<gt> binary_matrix_transposed.csv>

=back

B<or>

=over

=item C<mkdir result_dir && ./prot_finder_pipe.sh -q query.faa -s subject.faa -d result_dir -m | tee result_dir/blast_hits.tsv | perl prot_binary_matrix.pl | tee result_dir/binary_matrix.tsv | perl transpose_matrix.pl E<gt> result_dir/binary_matrix_transposed.tsv>

=back

=head1 VERSION

 0.1                                                       12-04-2016

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

=head1 ACKNOWLEDGEMENT

The Perl implementation for transposing a matrix on Stack Overflow
was very useful:
L<https://stackoverflow.com/questions/1729824/transpose-a-file-in-bash>

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

### Get the options with Getopt::Long
my $Delimiter = "\t"; # set separator/delimiter of input/output matrix
my $Empty; # optionally, fill empty cells with a value
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('delimiter=s' => \$Delimiter,
            'empty=s' => \$Empty,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);


### Run perldoc on POD and set option defaults
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);


### Check input
if (-t STDIN && ! @ARGV) {
    my $warning = "\n### Fatal error: No STDIN and no input file given as argument, please supply one of them and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
} elsif (!-t STDIN && @ARGV) {
    my $warning = "\n### Fatal error: Both STDIN and an input file given as argument, please supply only either one and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
}
die "\n### Fatal error: Too many arguments given, only STDIN or one input file allowed as argument! Please see the usage with option '-h' if unclear!\n\n" if (@ARGV > 1);
die "\n### Fatal error: File '$ARGV[0]' does not exist!\n\n" if (@ARGV && $ARGV[0] ne '-' && !-e $ARGV[0]);


### Parse input matrix
my %Input_Matrix; # hash of hash to store the input matrix
my $Max_Columns = 0; # maximum number of columns, needed in case not every row of input matrix has the same number of columns
my $Row_Num = 0; # count input matrix number of rows
while (<>) {
    chomp;
    warn "### Warning: Set separator/delimiter '$Delimiter' (option '-d') not found in the following first line/header of input matrix, sure the correct one is set?\n$_\n\n" if ($_ !~ /$Delimiter/ && $. == 1);

    my $col_num = 0; # count number of columns for each row
    foreach my $cell (split(/$Delimiter/)) { # split each row for the cells
        $cell = $Empty if ($cell =~ /^$/); # needed for empty cells in between cells with values, because for these $cell is defined in print out below
        $Input_Matrix{$Row_Num}{$col_num++} = $cell;
    }

    $Max_Columns = $col_num if ($col_num > $Max_Columns);
    $Row_Num++;
}


### Print out transposed matrix
my $Max_Rows = $Row_Num;
for (my $col_num = 0; $col_num < $Max_Columns; $col_num++) {
    for ($Row_Num = 0; $Row_Num < $Max_Rows; $Row_Num++) { # repurposing $Row_Num
        print "$Delimiter" if ($Row_Num > 0); # separator only after the first transposed column
        if (defined $Input_Matrix{$Row_Num}{$col_num}) { # 'defined' needed, in case $cell has '0' as value
            print $Input_Matrix{$Row_Num}{$col_num};
        } elsif (defined $Empty) { # for rows of the input matrix with columns < $Max_Columns; 'defined' needed, in case $Empty is set to '0'
            print $Empty;
        }
    }
    print "\n";
}

exit;
