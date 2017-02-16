#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<cdd2cog.pl> - assign COG categories to protein sequences

=head1 SYNOPSIS

C<perl cdd2cog.pl -r rps-blast.out -c cddid.tbl -f fun.txt -w whog>

=head1 DESCRIPTION

The script assigns COG (L<cluster of orthologous
groups|http://www.ncbi.nlm.nih.gov/COG/>) categories to proteins.
For this purpose, the query proteins need to be blasted with
RPS-BLAST+ (L<Reverse Position-Specific BLAST|http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>)
against NCBI's Conserved Domain Database
(L<CDD|http://www.ncbi.nlm.nih.gov/cdd>). Use
L<C<cds_extractor.pl>|/cds_extractor> beforehand to extract multi-fasta
protein files from GENBANK or EMBL files.

Both tab-delimited RPS-BLAST+ outformats, B<-outfmt 6> and B<-outfmt
7>, can be processed by C<cdd2cog.pl>. By default, RPS-BLAST+ hits
for each query protein are filtered for the best hit (lowest
e-value). Use option B<-a|all_hits> to assign COGs to all BLAST hits
and e.g. do a downstream filtering in a spreadsheet application.
Results are written to tab-delimited files in the F<./results>
folder, overall assignment statistics are printed to C<STDOUT>.

Several files are needed from NCBI's FTP server to run the RPS-BLAST+
and C<cdd2cog.pl>:

=over

=item 1.) L<CDD|ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/>

More information about the files in the CDD FTP archive can be found
in the respective F<README> file.

=item 1.1.) F<cddid.tbl.gz>

The file needs to be unpacked:

C<gunzip cddid.tbl.gz>

Contains summary information about the CD models in a tab-delimited
format. The columns are: PSSM-Id, CD accession (e.g. COG#), CD short
name, CD description, and PSSM (position-specific scoring matrices)
length.

=item 1.2.) F<./little_endian/Cog_LE.tar.gz>

Unpack and untar via:

C<tar xvfz Cog_LE.tar.gz>

Preformatted RPS-BLAST+ database of the CDD COG distribution for
Intel CPUs and Unix/Windows architectures.

=item 2.) L<COG|ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/>

Read F<readme> for more information about the respective files in
the COG FTP archive.

=item 2.1.) F<fun.txt>

One-letter functional classification used in the COG database.

=item 2.2.) F<whog>

Name, description, and corresponding functional classification of
each COG.

=back

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-r>=I<str>, B<-rps_report>=I<str>

Path to RPS-BLAST+ report/output, outfmt 6 or 7

=item B<-c>=I<str>, B<-cddid>=I<str>

Path to CDD's F<cddid.tbl> file

=item B<-f>=I<str>, B<-fun>=I<str>

Path to COG's F<fun.txt> file

=item B<-w>=I<str>, B<-whog>=I<str>

Path to COG's F<whog> file

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-a>, B<-all_hits>

Don't filter RPS-BLAST+ output for the best hit, rather assign COGs
to all hits

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

Overall assignment statistics

=item F<./results>

All tab-delimited output files are stored in this result folder

=item F<rps-blast_cog.txt>

COG assignments concatenated to the RPS-BLAST+ results for filtering

=item F<protein-id_cog.txt>

Slimmed down F<rps-blast_cog.txt> only including query id (first
BLAST report column), COGs, and functional categories

=item F<cog_stats.txt>

Assignment counts for each used COG

=item F<func_stats.txt>

Assignment counts for single-letter functional categories

=back

=head1 EXAMPLES

=head2 RPS-BLAST+

=over

=item C<rpsblast -query protein.fasta -db Cog -out rps-blast.out
-evalue 1e-2 -outfmt 6>

=item C<rpsblast -query protein.fasta -db Cog -out rps-blast.out
-evalue 1e-2 -outfmt '7 qseqid sseqid pident length mismatch gapopen
qstart qend sstart send evalue bitscore qcovs'>

=back

=head2 C<cdd2cog.pl>

=over

=item C<perl cdd2cog.pl -r rps-blast.out -c cddid.tbl -f fun.txt
-w whog -a>

=back

=head1 VERSION

 0.2                                               update: 2017-02-16
 0.1                                                       2013-08-01

=head1 AUTHOR

 Andreas Leimbach                         aleimba[at]gmx[dot]de

=head1 ACKNOWLEDGEMENTS

I got the idea for using NCBI's CDD PSSMs for COG assignment from JGI's L<IMG/ER annotation
system|http://img.jgi.doe.gov/>, which employes the same technique.


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
my $Rps_Report; # path to the rps-blast report/output
my $CDDid_File; # path to the CDD 'cddid.tbl' file
my $Fun_File; # path to the COG 'fun' file
my $Whog_File; # path to the COG 'whog' file
my $Opt_All_Hits; # give all blast hits for a query, not just the best (lowest evalue)
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('rps_report=s' => \$Rps_Report,
            'cddid=s' => \$CDDid_File,
            'fun=s' => \$Fun_File,
            'whog=s' => \$Whog_File,
            'all_hits' => \$Opt_All_Hits,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Rps_Report || !$CDDid_File || !$Fun_File || !$Whog_File) {
    my $warning = "\n### Fatal error: Option(s) or arguments for '-r', '-c', '-f', or '-w' are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Parse the 'cddid.tbl', 'fun.txt' and 'whog' file contents and store info in hash structures
my (%CDDid, %Fun, %Whog); # global hashes
parse_cdd_cog(); # subroutine



### Create results directory for output files
my $Out_Dir = './results/';
if (-e $Out_Dir) {
    print "###Directory '$Out_Dir' already exists! Replace the directory and all its contents [y|n]? ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
        unlink glob "$Out_Dir*"; # remove all files in results directory
        rmdir $Out_Dir; # remove the empty directory
    } else {
        die "Script abborted!\n";
    }
}
mkdir $Out_Dir or die "Can't create directory \"$Out_Dir\": $!\n";



### Parse the rps-blast report/output file and assign COGs
my %Cog_Stats; # store the total number of query protein hits for each COG, written to '$Cogstats_Out' below

my $Blast_Out = 'rps-blast_cog.txt'; # output file for COG assignments appended to RPS-BLAST results
open (my $Blast_Out_Fh, ">", "$Out_Dir"."$Blast_Out");
print $Blast_Out_Fh "query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tCOG#\tfunctional categories\t\t\t\t\tCOG protein description\n"; # header for $Blast_Out

my $Locus_Cog = "protein-id_cog.txt"; # slimmed down $Blast_Out only including locus_tags, COGs, and functional categories
open (my $Locus_Cog_Fh, ">", "$Out_Dir"."$Locus_Cog");

print "Parsing RPS-BLAST report ...\n"; # status message
my $Skip = ''; # only keep best blast hit per query (lowest e-value), except option 'all_hits' is given
open (my $Rps_Report_Fh, "<", "$Rps_Report");
while (<$Rps_Report_Fh>) {
    if (/^#/) { # skip comment lines in blast report for BLAST+ "outfmt 7"
        next;
    }
    chomp;

    my @line = split(/\t/, $_); # split tab-separated RPS-BLAST report

    if ($Skip eq $line[0] && !$Opt_All_Hits) {
        # only keep best blast hit per query, only if option 'all_hits' is NOT set
        # $line[0] is query id, and should be locus_tag or specific ID from multi-fasta protein query file
        next;
    }
    $Skip = $line[0];

    my $pssm_id = $1 if $line[1] =~ /^CDD\:(\d+)/; # get PSSM-Id from the subject hit
    my $cog = $CDDid{$pssm_id}; # get the COG# according to the PSSM-Id as listed in 'cddid.tbl'
    $Cog_Stats{$cog}++; # increment hit-number for specific COG

    ### Collect functional categories stats
    my @functions = split('', $Whog{$cog}->{'function'}); # split the single-letter functional categories to count them and join them as tab-separated below
    foreach (@functions) {
        $Fun{$_}->{'count'}++; # increment hit-number for specific functional category
    }

    ### Print to result files
    my $functions = join("\t", @functions); # join functional categories tab-separated
    print $Locus_Cog_Fh "$line[0]\t$cog\t$functions\n"; # locus_tag\tCOG\tfunctional categories
    $functions .= "\t" x (5 - @functions); # add additional tabs for COGs with fewer than five functions (which is the maximum number)
    print $Blast_Out_Fh "$_\t$cog\t$functions\t$Whog{$cog}->{'desc'}\n"; # $_ = RPS-BLAST line
}

close $Rps_Report_Fh;
close $Blast_Out_Fh;
close $Locus_Cog_Fh;



### Total COG and functional categories stats
print "Writing assignment statistic files in '$Out_Dir' folder ...\n"; # status message

my $Cogstats_Out = 'cog_stats.txt'; # output file for assignment numbers for each COG
open (my $Cog_Stats_Fh, ">", "$Out_Dir"."$Cogstats_Out");
my $prot_stats = 0; # store total number of query proteins, which have a COG assignment
foreach my $cog (sort keys %Cog_Stats) {
    print $Cog_Stats_Fh "$cog\t$Whog{$cog}->{'desc'}\t$Cog_Stats{$cog}\n"; # COG protein descriptions stored in %Whog
    $prot_stats += $Cog_Stats{$cog}; # sum up total COG assignments
}
close $Cog_Stats_Fh;

my $Funcstats_Out = 'func_stats.txt'; # output file for assignment numbers for each functional category
open (my $Func_Stats_Fh, ">", "$Out_Dir"."$Funcstats_Out");
my $func_cats = 0; # store total number of assigned functional categories
foreach my $func (sort keys %Fun) {
    print $Func_Stats_Fh "$func\t$Fun{$func}->{'desc'}\t$Fun{$func}->{'count'}\n";
    $func_cats += $Fun{$func}->{'count'}; # sum up total functional category assignments
}
close $Func_Stats_Fh;



### State which files were created and print overall statistics
print "\n############################################################################\n";
print "The following tab-delimited files were created in the '$Out_Dir' directory:\n";
print "- $Blast_Out: COG assignments concatenated to the RPS-BLAST results for filtering\n";
print "- $Locus_Cog: Slimmed down '$Blast_Out' only including query id (first BLAST report column), COGs, and functional categories\n";
print "- $Cogstats_Out: COG assignment counts\n";
print "- $Funcstats_Out: Functional category assignment counts\n";
print "##############################################################################\n";
print "Overall assignment statistics:\n";
print "~ Total query proteins categorized into COGs: $prot_stats\n";
print "~ Total COGs used for the query proteins [of ", scalar keys %CDDid, " overall]: ", scalar keys %Cog_Stats, "\n";
print "~ Total number of assigned functional categories: $func_cats\n";
print "~ Total functional categories used for the query proteins [of ", scalar keys %Fun, " overall]: ", scalar grep ($Fun{$_}->{'count'} > 0, keys %Fun), "\n\n"; # grep for functional categories with a count > 0 to get the ones with assigned query proteins

exit;



###############
# Subroutines #
###############

### Subroutine to parse the 'cddid.tbl', 'fun' and 'whog' file contents and store in hash structures
sub parse_cdd_cog {

    ### 'cddid.tbl'
    open (my $cddid_fh, "<", "$CDDid_File");
    print "\nParsing CDDs '$CDDid_File' file ...\n"; # status message
    while (<$cddid_fh>) {
        chomp;
        my @line = split(/\t/, $_); # split line at the tabs
        if ($line[1] =~ /^COG\d{4}$/) { # search for COG CD accessions in cddid
            $CDDid{$line[0]} = $line[1]; # hash to store info; $line[0] = PSSM-Id
        }
    }
    close $cddid_fh;

    ### 'fun.txt'
    open (my $fun_fh, "<", "$Fun_File");
    print "Parsing COGs '$Fun_File' file ...\n"; # status message
    while (<$fun_fh>) {
        chomp;
        $_ =~ s/^\s*|\s+$//g; # get rid of all leading and trailing whitespaces
        if (/^\[(\w)\]\s*(.+)$/) {
            $Fun{$1} = {'desc' => $2, 'count' => 0}; # anonymous hash in hash
            # $1 = single-letter functional category, $2 = description of functional category
            # count used to find functional categories not present in the query proteins for final overall assignment statistics
        }
    }
    close $fun_fh;

    ### 'whog'
    open (my $whog_fh, "<", "$Whog_File");
    print "Parsing COGs '$Whog_File' file ...\n"; # status message
    while (<$whog_fh>) {
        chomp;
        $_ =~ s/^\s*|\s+$//g; # get rid of all leading and trailing whitespaces
        if (/^\[(\w+)\]\s*(COG\d{4})\s+(.+)$/) {
            $Whog{$2} = {'function' => $1, 'desc' => $3}; # anonymous hash in hash
            # $1 = single-letter functional categories, maximal five per COG (only COG5032 with five)
            # $2 = COG#, $3 = COG protein description
        }
    }
    close $whog_fh;

    return 1;
}
