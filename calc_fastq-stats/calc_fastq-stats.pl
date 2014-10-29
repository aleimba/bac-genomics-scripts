#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<calc_fastq-stats.pl> - basic statistics for bases and reads in a FASTQ file

=head1 SYNOPSIS

C<perl calc_fastq-stats.pl -i reads.fastq>

B<or>

C<gzip -dc reads.fastq.gz | perl calc_fastq-stats.pl -i ->

=head1 DESCRIPTION

Calculates some simple statistics, like individual and total base
counts, GC content, and basic stats for the read lengths, and
read/base qualities in a FASTQ file. The GC content calculation does
not include 'N's. Stats are printed to C<STDOUT> and optionally to an
output file.

Because the quality of a read degrades over its length with all NGS
machines, it is advisable to also plot the quality for each cycle as
implemented in tools like
L<FastQC|http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
or the L<fastx-toolkit|http://hannonlab.cshl.edu/fastx_toolkit/>.

If the sequence and the quality values are interrupted by line
breaks (i.e. a read is B<not> represented by four lines), please fix
with Heng Li's L<seqtk|https://github.com/lh3/seqtk>:

C<seqtk seq -l 0 infile.fastq E<gt> outfile.fastq>

An alternative tool, which is a lot faster, is B<fastq-stats> from
L<ea-utils|https://code.google.com/p/ea-utils/>.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input FASTQ file or piped STDIN (-) from a gzipped file

=item B<-q>=I<int>, B<-qual_offset>=I<int>

ASCII quality offset of the Phred (Sanger) quality values [default 33]

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-c>=I<int>, B<-coverage_limit>=I<int>

Number of bases to sample from the top of the file

=item B<-n>=I<int>, B<-num_read>=I<int>

Number of reads to sample from the top of the file

=item B<-o>=I<str>, B<-output>=I<str>

Print stats in addition to C<STDOUT> to the specified output file

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

Calculated stats are printed to C<STDOUT>

=item F<(outfile)>

Optional outfile for stats

=back

=head1 EXAMPLES

=over

=item C<zcat reads.fastq.gz | perl calc_fastq-stats.pl -i - -q 64 -c 175000000 -n 3000000>

=back

=head1 DEPENDENCIES

If the following modules are not installed get them from
L<CPAN|http://www.cpan.org/>:

=over 20

=item C<Statistics::Descriptive>

Perl module to calculate basic descriptive statistics

=item C<Statistics::Descriptive::Discrete>

Perl module to calculate descriptive statistics for discrete data sets

=item C<Statistics::Descriptive::Weighted>

Perl module to calculate descriptive statistics for weighted variates

=back

=head1 VERSION

 0.1                                                       28-10-2014

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

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

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Pod::Usage;
eval { # check if modules are installed
    require Statistics::Descriptive; # module of basic descriptive statistics
    Statistics::Descriptive->import();
    require Statistics::Descriptive::Discrete; # module for statistics with 'discrete' values
    Statistics::Descriptive::Discrete->import();
    require Statistics::Descriptive::Weighted; # module for weighted statistics (includes approximations for quantiles in contrast to 'Statistics::Descriptive::Discrete')
    Statistics::Descriptive::Weighted->import();
}; # semi-colon needed
if ($@) { # if module(s) not installed die with error message
    die "\n### Fatal error: One or several of the required statistical Perl modules 'Statistics::Descriptive', 'Statistics::Descriptive::Weighted', or 'Statistics::Descriptive::Discrete' are not installed! Please install from CPAN!\n\n";
}

### Get the options with Getopt::Long
my $Input_File; # either input FASTQ file or STDIN
my $Output_File; # stats are printed to STDOUT, optionally to this output file
my $Qual_Ascii_Offset = 33; # ASCII quality offset for Phred (Sanger) quality scores (see below)
my $Coverage_Limit; # number of bases to sample
my $Num_Read; # number of reads to sample
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_File,
            'output:s' => \$Output_File,
            'qual_offset=i' => \$Qual_Ascii_Offset,
            'coverage_limit:i' => \$Coverage_Limit,
            'num_read:i' => \$Num_Read,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_File) {
    my $warning = "\n### Fatal error: Option '-i' or its argument is missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Open FASTQ file or accept STDIN, open optional output file
my $Input_Fh;
if ($Input_File eq '-') { # file input via STDIN
    $Input_Fh = *STDIN; # capture typeglob of STDIN
} else { # input via FASTQ file
    open ($Input_Fh, "<", "$Input_File");
}
my $Output_Fh;
if ($Output_File) {
    file_exist($Output_File); # subroutine
    open ($Output_Fh, ">", "$Output_File");
}



### Parse read data in FASTQ file
my ($A, $C, $G, $T, $U, $N);
$A = $C = $G = $T = $U = $N = 0;
my $Total_Bases = 0; # total number of bases

my $Line_Count = 0; # number of lines in the FASTQ file
my $Read_Count = 0; # number of reads
my @Read_Lengths; # length of each individual read

my %Base_Quality_Counts; # quality counts for each individual base
my @Mean_Read_Quals; # average qualities for each read

my $I = 1; # mio. counter for processing status message
while (<$Input_Fh>) { # FASTQ "usually" format uses four lines per sequence
    chomp;
    $Read_Count++;
    $Line_Count++;

    # status message
    if ($Read_Count/1000000 == $I) {
        print STDERR "$I Mio reads processed ...\r"; # carriage return to overwrite messages and not clutter STDERR
        $I++;
    }

    # sequence identifier/read name (line 1)
    # e.g.: @H108:287:D0M79ACXX:7:1101:5248:1997 1:N:0:TAGGCATGAGAGTAGA
    # @ <instrument‐name>:<run ID>:<flowcell ID>:<lane‐number>:<tile‐number>:<x‐pos>:<y‐pos> <read number>:<is filtered>:<control number>:<barcode sequence>
    # [Illumina, CASAVA v1.8 Changes, 05.01.2011]
    my $seq_id = $_;
    die "\n### Fatal error:\nThis read doesn't have a sequence identifier/read name according to FASTQ specs, it should begin with a '\@':\n$seq_id\n" if ($seq_id !~ /^@/);
    $seq_id =~ s/^@(.+)/$1/; # remove '@' to make comparable to $qual_id

    # DNA sequence of the read (line 2)
    my $seq = read_line(); # subroutine
    die "\n### Fatal error:\nRead '$seq_id' has a whitespace in its sequence, which is not allowed according to FASTQ specs:\n$seq\n" if ($seq =~ /\s+/);
    base_count($seq, $seq_id); # subroutine
    push(@Read_Lengths, length $seq);
    $Total_Bases += length $seq;

    # optional sequence ID of quality (line 3)
    my $qual_id = read_line(); # subroutine
    die "\n### Fatal error:\nThe optional sequence identifier/read name for the quality line of read '$seq_id' is not according to FASTQ specs, it should begin with a '+':\n$qual_id\n" if ($qual_id !~ /^\+/);
    $qual_id =~ s/^\+(.*)/$1/; # if optional ID is present check if equal to $seq_id in line 1
    die "\n### Fatal error:\nThe sequence identifier/read name of read '$seq_id' doesn't fit to the optional ID in the quality line:\n$qual_id\n" if ($qual_id && $qual_id ne $seq_id);

    # ASCII quality values (line 4)
    # The quality scores (Q; from Illumina 1.3 onwards: Phred [i.e. Sanger] Q = -10 log10(Pe); with Pe being estimated probability for base calling error) are transformed from integer to a single ASCII character for brevity so that a string can represent all of the quality scores within a read. ASCII offset of 33 (Illumina 1.8+; from Illumina 1.3+ was offset of 64). Q + 33 = ASCII. [Illumina, CASAVA v1.8 Changes, 05.01.2011]
    my $qual = read_line(); # subroutine
    die "\n### Fatal error:\nRead '$seq_id' has a whitespace in its quality values, which is not allowed according to FASTQ specs:\n$qual\n" if ($qual =~ /\s+/);
    die "\n### Fatal error:\nRead '$seq_id' has a non-ASCII character in its quality values, which is not allowed according to FASTQ specs:\n$qual\n" if ($qual =~ /[^[:ascii]]/);
    die "\n### Fatal error:\nThe quality line of read '$seq_id' doesn't have the same number of symbols as letters in the sequence:\n$seq\n$qual\n" if (length $qual != length $seq);

    my $total_read_qual = 0;
    foreach (split('', $qual)) {
        my $phred_qual = ord($_)-$Qual_Ascii_Offset; # ord converts a character to its ASCII value (an integer), which has to be offset to get the Phred quality value
        $Base_Quality_Counts{$phred_qual}++; # count individual base qualities
        $total_read_qual += $phred_qual;
    }

    my $mean_read_qual = sprintf("%.2f", $total_read_qual/length $qual);
    push(@Mean_Read_Quals, $mean_read_qual);

    if ($Coverage_Limit && $Total_Bases >= $Coverage_Limit) {
        print STDERR "\nReached coverage limit of '$Coverage_Limit bp' at read number '$Read_Count' with '$Total_Bases bp'!\n";
        last;
    }
    last if ($Num_Read && $Read_Count == $Num_Read); # skip rest of reads if $num_read is reached
}
print "\n"; # get rid of the carriage return for status messages
close $Input_Fh;



### Final sanity checks
die "\n### Fatal error:\nFASTQ file doesn't contain four lines for each read, which is an accepted convention in the FASTQ specs! If the sequence and the quality values are interrupted by line breaks please fix with Heng Li's 'seqtk' (https://github.com/lh3/seqtk):\nseqtk seq -l 0 infile.fastq > outfile.fastq\n" if ($Line_Count/4 != $Read_Count); # should be already covered as the parser assumes four lines per read (see above)
my $Total_Nucs = $A + $C + $G + $T + $U + $N;
die "\n### Fatal error:\nThe total number of bases '$Total_Nucs bp' is not equal to the read length sum '$Total_Bases bp'!\nAre there degenerate IUPAC bases (except for 'N') or non-nucleotide characters in the sequence (which is not allowed according to FASTQ specs)?\n" if ($Total_Nucs != $Total_Bases); # should be already covered in subroutine 'base_count' and checks in the parser above



### Read length stats
my %Read_Lengths_Stats; # store stats for read lengths
stats_full(\@Read_Lengths, \%Read_Lengths_Stats); # subroutine for 'Statistics::Descriptive'
print "#Read stats:\n";
print_out_std("number_of_reads\t$Read_Count\n"); # subroutine
print_out_std("min_read_length\t$Read_Lengths_Stats{'min'}\n");
print_out_std("max_read_length\t$Read_Lengths_Stats{'max'}\n");
print_out_std("mean_read_length\t$Read_Lengths_Stats{'mean'}\n");
print_out_std("q1_read_length\t$Read_Lengths_Stats{'q1'}\n");
print_out_std("median_read_length\t$Read_Lengths_Stats{'median'}\n");
print_out_std("q3_read_length\t$Read_Lengths_Stats{'q3'}\n");
print_out_std("stdev_read_length\t$Read_Lengths_Stats{'sd'}\n");



### Mean read Phred quality stats
my %Mean_Read_Quals_Stats; # store mean quality stats for reads
stats_full(\@Mean_Read_Quals, \%Mean_Read_Quals_Stats); # subroutine for 'Statistics::Descriptive'
print "#Read quality stats:\n";
print_out_std("min_read_qual\t$Mean_Read_Quals_Stats{'min'}\n"); # subroutine
print_out_std("max_read_qual\t$Mean_Read_Quals_Stats{'max'}\n");
print_out_std("mean_read_qual\t$Mean_Read_Quals_Stats{'mean'}\n");
print_out_std("q1_read_qual\t$Mean_Read_Quals_Stats{'q1'}\n");
print_out_std("median_read_qual\t$Mean_Read_Quals_Stats{'median'}\n");
print_out_std("q3_read_qual\t$Mean_Read_Quals_Stats{'q3'}\n");
print_out_std("stdev_read_qual\t$Mean_Read_Quals_Stats{'sd'}\n");



### Total bases Phred quality stats
my %Base_Quals_Stats; # store quality stats for individual bases
stats_discrete_full(\%Base_Quality_Counts, \%Base_Quals_Stats); # subroutine for 'Statistics::Descriptive::Discrete' and 'Statistics::Descriptive::Weighted'
print "#Base quality stats:\n";
print_out_std("min_base_qual\t$Base_Quals_Stats{'min'}\n"); # subroutine
print_out_std("max_base_qual\t$Base_Quals_Stats{'max'}\n");
print_out_std("mean_base_qual\t$Base_Quals_Stats{'mean'}\n");
print_out_std("q1_base_qual\t$Base_Quals_Stats{'q1'}\n");
print_out_std("median_base_qual\t$Base_Quals_Stats{'median'}\n");
print_out_std("q3_base_qual\t$Base_Quals_Stats{'q3'}\n");
print_out_std("stdev_base_qual\t$Base_Quals_Stats{'sd'}\n");



### Base stats
print "#Base count stats:\n";
print_out_std("number_of_bases\t$Total_Bases\n"); # subroutine
print_out_std("A_count\t$A\n");
print_out_std("C_count\t$C\n");
print_out_std("G_count\t$G\n");
print_out_std("T_count\t$T\n");
if ($U) {
    print STDERR "Nucleotide 'U' found in the sequences, RNA?\n" if ($U);
    print_out_std("U_count\t$U\n") if ($U);
}
print_out_std("N_count\t$N\n") if ($N);
my $GC_Content = sprintf ("%.2f", (($C + $G)/$Total_Bases)*100);
print_out_std("GC_content_[%]\t$GC_Content\n");
close Output_Fh if ($Output_File);


exit;

#############
#Subroutines#
#############

### Subroutine to count bases
sub base_count {
    my ($seq, $seq_id) = @_;
    die "\n### Fatal error:\nRead '$seq_id' has a IUPAC degenerate base (except for 'N') or non-nucleotide character in its sequence, which is not allowed according to FASTQ specs:\n$seq\n" if ($seq =~ /[^acgtun]/i);
    $A += ($seq =~ tr/[aA]//);
    $C += ($seq =~ tr/[cC]//);
    $G += ($seq =~ tr/[gG]//);
    $T += ($seq =~ tr/[tT]//);
    $U += ($seq =~ tr/[uU]//);
    $N += ($seq =~ tr/[nN]//);
    return 1;
}



### Subroutine to test for file existence and give warning to STDERR
sub file_exist {
    my $file = shift;
    if (-e $file) {
        warn "\nThe result file '$file' exists already and will be overwritten!\n\n";
        return 1;
    }
    return 0;
}



### Subroutine to read in a line from input, chomp, and count it
sub read_line {
    my $line = <$Input_Fh>;
    $Line_Count++;
    chomp $line;
    return $line;
}



### Subroutine to calculate stats with modules 'Statistics::Descriptive::Discrete' and 'Statistics::Descriptive::Weighted' for large data sets, which overburden RAM and/or module 'Statistics::Descriptive'
sub stats_discrete_full {
    my ($data_hash_ref, $hash_stat_ref) = @_;

    # convert data to needed formats
    my @discrete_tuple; # for 'Statistics::Descriptive::Discrete'
    my (@values, @weights); # for 'Statistics::Descriptive::Weighted'
    foreach (sort {$a <=> $b} keys %$data_hash_ref) { # de-reference hash ref
        push(@discrete_tuple, ($_, $data_hash_ref->{$_})); # expects first value then weight, de-reference
        push(@values, $_);
        push(@weights, $data_hash_ref->{$_});
    }

    my $stat_discrete = Statistics::Descriptive::Discrete->new();
    $stat_discrete->add_data_tuple(@discrete_tuple);
    $hash_stat_ref->{'mean'} = sprintf("%.2f", $stat_discrete->mean()); # rounded to two decimal places
    $hash_stat_ref->{'median'} = sprintf("%.2f", $stat_discrete->median());
    $hash_stat_ref->{'sd'} = sprintf("%.2f", $stat_discrete->standard_deviation());
    $hash_stat_ref->{'min'} = $stat_discrete->min();
    $hash_stat_ref->{'max'} = $stat_discrete->max();

    # quantile approximations only implemented in 'Statistics::Descriptive::Weighted' not in 'Statistics::Descriptive::Discrete'
    my $stat_weighted = Statistics::Descriptive::Weighted::Full->new();
    $stat_weighted->add_data(\@values, \@weights); # add data as array refs
    $hash_stat_ref->{'q1'} = sprintf("%.2f", $stat_weighted->quantile(.25)); # Q1, first quartile (25th percentile)
    $hash_stat_ref->{'q3'} = sprintf("%.2f", $stat_weighted->quantile(.75)); # Q3, third quartile (75th percentile)
    return 1;
}



### Subroutine to calculate stats with module 'Statistics::Descriptive'
sub stats_full {
    my ($data_array_ref, $hash_stat_ref) = @_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data_array_ref); # de-reference array ref
    $hash_stat_ref->{'mean'} = sprintf("%.2f", $stat->mean());
    $hash_stat_ref->{'median'} = sprintf("%.2f", $stat->median());
    $hash_stat_ref->{'sd'} = sprintf("%.2f", $stat->standard_deviation());
    $hash_stat_ref->{'min'} = $stat->min();
    $hash_stat_ref->{'max'} = $stat->max();
    $hash_stat_ref->{'q1'} = sprintf("%.2f", $stat->quantile(1)); # Q1, first quartile (25th percentile)
    $hash_stat_ref->{'q3'} = sprintf("%.2f", $stat->quantile(3)); # Q3, third quartile (75th percentile)
    # $stat->clear(); # remove stats from the module for next round; not needed because of 'new()' initialisation at begin of subroutine
    return 1;
}



### Subroutine to print to STDOUT and optionally to output file
sub print_out_std {
    print "@_"; # print line to STDOUT
    print $Output_Fh "@_" if ($Output_File); # additionally, print to optional output file
    return 1;
}
