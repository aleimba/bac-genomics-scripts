#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<prot_binary_matrix.pl> - create a presence/absence matrix from
C<prot_finder.pl> output

=head1 SYNOPSIS

C<perl prot_binary_matrix.pl blast_hits.tsv E<gt> binary_matrix.tsv>

B<or>

C<perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl E<gt> binary_matrix.tsv>

=head1 DESCRIPTION

This script is intended to create a presence/absence matrix from the
significant C<prot_finder.pl>
L<B<BLASTP>|http://blast.ncbi.nlm.nih.gov/Blast.cgi>) hits (or the
companion bash pipe C<prot_finder_pipe.sh>). The tab-separated
C<prot_finder.pl> output can be given directly via C<STDIN> or as a
file. By default a tab-delimited binary presence/absence matrix for
query hits per subject organism will be printed to C<STDOUT>. Use
option B<-t> to count all query hits per subject organism, not just
the binary presence/absence. You can transpose the presence/absence
binary matrix with the script C<transpose_matrix.pl> (see its help
with B<-h>).

The resulting matrix can be used to associate the presence/absence
data with a phylogenetic tree, e.g. use the Interactive Tree Of Life
website (L<B<iTOL>|http://itol.embl.de/>). B<iTOL> likes individual
comma-separated input files, thus use options B<-s -c> for this
purpose.

For B<iTOL> the organism names have to have identical names to the
leaves of the phylogenetic tree, thus manual adaptation, e.g. in a
spreadsheet software, might be needed. B<Careful>, subject organisms
without a significant B<BLASTP> hit won't be included in the
tab-separated C<prot_finder.pl> result table and hence can't be
included by C<prot_binary_matrix.pl>. If needed add them manually to
the result matri(x|ces).

Additionally, you can give the presence/absence binary matrix to
C<binary_group_stats.pl> to calculate presence/absence statistics
for groups of columns and not simply single columns of the matrix.
C<binary_group_stats.pl> also has a comprehensive manual with its
option B<-h>.

=head1 OPTIONS

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-s>, B<-separate>

Separate presence/absence files for each query protein printed to
the result directory [default without B<-s> = C<STDOUT> matrix for
all query proteins combined]

=item B<-d>=I<str>, B<-dir_result>=I<str>

Path to result folder, requires option B<-s> [default =
'./binary_matrix_results']

=item B<-t>, B<-total>

Count total occurrences of query proteins, not just binary
presence/absence

=item B<-c>, B<-csv>

Output matri(x|ces) in comma-separated format (*.csv) instead of
tab-delimited format (*.tsv)

=item B<-l>, B<-locus_tag>

Use the locus_tag B<prefixes> in the subject_ID column of the
C<prot_finder.pl> output (instead of the subject_organism columns) as
organism IDs to associate query hits to organisms. The subject_ID
column will include locus_tags if they're annotated for a genome
(see the L<C<cds_extractor.pl>|/cds_extractor> format description).
Useful if the L<C<cds_extractor.pl>|/cds_extractor> output doesn't
include strain names for 'o=' in the FASTA IDs, because the prefix
of a locus_tag should be unique for a genome (see
L<http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation>).

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 17

=item C<STDOUT>

The resulting presence/absence matrix is printed to C<STDOUT>
without option B<-s>. Redirect or pipe into another tool as needed.

=item (F<./binary_matrix_results>)

Separate query presence/absence files are stored in a result folder
with option B<-s>

=item (F<./binary_matrix_results/query-ID_binary_matrix.(tsv|csv)>)

Separate query presence/absence files with option B<-s>

=back

=head1 EXAMPLES

=over

=item C<perl prot_binary_matrix.pl -s -d result_dir -t blast_hits.tsv>

=back

B<or>

=over

=item C<perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl -l -c E<gt> binary_matrix.csv>

=back

B<or>

=over

=item C<mkdir result_dir && ./prot_finder_pipe.sh -q query.faa -s subject.faa -d result_dir -m | tee result_dir/blast_hits.tsv | perl prot_binary_matrix.pl E<gt> binary_matrix.tsv>

=back

=head1 VERSION

 0.6                                               update: 23-11-2015
 0.1                                                       25-10-2012

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

### Get options with Getopt::Long
my $Opt_Separate; # separate presence/absence files for each query printed to result_dir (default: single presence/absence file for all queries printed to STDOUT)
my $Result_Dir; # path to result folder, requires option '-s'; default set below to 'binary_matrix_results'
my $Opt_Total; # count total occurrences of query proteins not just presence/absence binary
my $Opt_Csv; # output in csv format (default: tsv)
my $Opt_Locus_Tag; # use locus_tag prefixes (from subject_ID column, see cds_exractor) instead of subject_organism as ID to count query hits
my $VERSION = 0.6;
my ($Opt_Version, $Opt_Help);
GetOptions ('separate' => \$Opt_Separate,
            'dir_result=s' => \$Result_Dir,
            'total' => \$Opt_Total,
            'csv' => \$Opt_Csv,
            'locus_tag' => \$Opt_Locus_Tag,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD and set option defaults
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if ($Result_Dir && !$Opt_Separate) {
    warn "### Warning: Option '-d' given but not its required option '-s', forcing option '-s'!\n";
    $Opt_Separate = 1;
}

my $Separator = "\t";
$Separator = "," if ($Opt_Csv); # optional csv output format



### Check input
if (-t STDIN && ! @ARGV) {
    my $warning = "\n### Fatal error: No STDIN and no input file given as argument, please supply one of them and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
} elsif (!-t STDIN && @ARGV) {
    my $warning = "\n### Fatal error: Both STDIN and an input file given as argument, please supply only either one and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
}
die "\n### Fatal error: Too many arguments given, only STDIN or one input file allowed as argument! Please see the usage with option '-h' if unclear!\n" if (@ARGV > 1);
die "\n### Fatal error: File '$ARGV[0]' does not exist!\n" if (@ARGV && $ARGV[0] ne '-' && !-e $ARGV[0]);



### Create result folder, only for option '-s'
if ($Opt_Separate) {
    $Result_Dir = 'binary_matrix_results' if (!$Result_Dir);
    $Result_Dir =~ s/\/$//; # get rid of a potential '/' at the end of $Result_Dir path
    if (-e $Result_Dir) {
        empty_dir($Result_Dir); # subroutine to empty a directory with user interaction
    } else {
        mkdir $Result_Dir;
    }
}



### Parse the input from 'prot_finder.pl' to associate organism with query hit
my @Queries; # store all query proteins
my %Hits; # hash of hash to associate subject_organism/subject_ID with query hit

while (<>) { # read STDIN or file input
    chomp;
    if ($. == 1) { # $. = check only first line of input (works with STDIN and file input)
        die "\n### Fatal error: Input doesn't have the correct format, it has to be the output of 'prot_finder.pl' with the following header:\n# subject_organism\tsubject_ID\tsubject_gene\tsubject_protein_desc\tquery_ID\tquery_desc\tquery_coverage [%]\tquery_identities [%]\tsubject/hit_coverage [%]\te-value of best HSP\tbit-score of best HSP\n" if (!/# subject_organism\tsubject_ID\tsubject_gene\tsubject_protein_desc\tquery_ID\tquery_desc/);
        next; # skip header line
    }

    my @line = split (/\t/, $_); # $line[0] = subject_organism; $line[1] = subject_ID (mostly locus_tag, see cds_extractor); $line[4] = query_ID
    my $query = $line[4];
    my $id;
    if ($Opt_Locus_Tag) { # use subject_ID
        die "\n### Fatal error: The subject_ID of the following line doesn't look like an NCBI locus tag (see: http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation). Column subject_ID needs to include only locus_tags for option '-l'!\n$_\n" if ($line[1] !~ /^([a-zA-Z][0-9a-zA-Z]{2,11})_[0-9a-zA-Z]+$/); # check if subject_ID is locus_tag ('\w' not used because allows alphanumeric and '_')
        # excerpt: The locus_tag prefix must be 3-12 alphanumeric characters and the first character may not be a digit. All chromosomes and plasmids of an individual genome must use the exactly same locus_tag prefix followed by an underscore and then an alphanumeric identification number that is unique within the given genome. Other than the single underscore used to separate the prefix from the identification number no other special characters can be used in the locus_tag.
        $id = $1; # locus_tag prefix, unique for each genome
    } else { # use subject_organism as ID
        $id = $line[0];
    }

    if ($Opt_Total) { # count total occurrences of query proteins
        $Hits{$id}{$query}++;

    } else { # only binary output (0 or 1)
        $Hits{$id}{$query} = 1;
    }

    push (@Queries, $query) if (!grep($_ eq $query, @Queries)); # push each query only once in @Queries
}



### Print binary data to a joined or separate (for each query; as needed by iTOL) file(s)
if (!$Opt_Separate) { # joined output
    # print header
    if ($Opt_Locus_Tag) {
        print "locus_tag";
    } else {
        print "organism";
    }
    print "$Separator";
    print join("$Separator", sort @Queries), "\n";

    # print data to STDOUT
    foreach my $id (sort keys %Hits) {
        print "$id";
        foreach my $query (sort @Queries) {
            if ($Hits{$id}->{$query}) {
                print "$Separator", "$Hits{$id}->{$query}";
            } else {
                print "$Separator", "0";
            }
        }
        print "\n";
    }

} elsif ($Opt_Separate) { # separated output for each query
    foreach my $query (sort @Queries) {
        my $file = "$Result_Dir/$query\_binary\_matrix.";
        if ($Opt_Csv) {
            $file .= "csv";
        } else {
            $file .= "tsv";
        }
        open (my $binary_matrix_fh, ">", "$file");
        foreach my $id (sort keys %Hits) {
            print $binary_matrix_fh "$id";
            if ($Hits{$id}->{$query}) {
                print $binary_matrix_fh "$Separator", "$Hits{$id}->{$query}\n";
            } else {
                print $binary_matrix_fh "$Separator", "0\n";
            }
        }
        close $binary_matrix_fh;
    }
}

exit;


###############
# Subroutines #
###############

### Subroutine to empty a directory with user interaction
sub empty_dir {
    my $dir = shift;
    print STDERR "\nDirectory '$dir' already exists! You can use either option '-d' to set a different output result directory name, or do you want to replace the directory and all its contents [y|n]? ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
        unlink glob "$dir/*"; # remove all files in results directory
    } else {
        die "\nScript abborted!\n";
    }
    return 1;
}
