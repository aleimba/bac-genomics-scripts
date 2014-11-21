#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<sam_insert-size.pl> - insert size and read length statistics for paired-end reads

=head1 SYNOPSIS

C<perl sam_insert-size.pl -i file.sam>

B<or>

C<samtools view -h file.bam | perl sam_insert-size.pl -i ->

=head1 DESCRIPTION

Calculate insert size and read length statistics for paired-end reads
in SAM/BAM alignment format. The program gives the arithmetic mean,
median, and standard deviation (stdev) among other statistical values.

Insert size is defined as the total length of the original fragment
put into sequencing, i.e. the sequenced DNA fragment between the
adaptors. The 16-bit FLAG of the SAM/BAM file is used to filter reads
(see the L<SAM
specifications|http://samtools.sourceforge.net/SAM1.pdf>).

B<Read length> statistics are calculated for all mapped reads
(irrespective of their pairing).

B<Insert size> statistics are calculated only for B<paired reads>.
Typically, the insert size is perturbed by artifacts, like chimeras,
structural re-arrangements or alignment errors, which result in a
very high maximum insert size measure. As a consequence the mean and
stdev can be strongly misleading regarding the real distribution. To
avoid this, two methods are implemented that first trim the insert
size distribution to a 'core' to calculate the respective statistics.
Additionally, secondary alignments for multiple mapping reads and
supplementary alignments for chimeric reads, as well as insert sizes
of zero are not considered (option B<-min_ins_cutoff> is set to
B<one> by default).

The B<-a|-align> method includes only proper/concordant paired reads
in the statistical calculations (as determined by the mapper and the
options for insert size minimum and maximum used for mapping). This
is the B<default> method.

The B<-p|-percentile> method first calculates insert size statistics
for all read pairs, where the read and the mate are mapped ('raw
data'). Subsequently, the 10th and the 90th percentile are discarded
to calculate the 10% truncated mean and stdev. Discarding the lowest
and highest 10% of insert sizes gives the advantage of robustness
(insensitivity to outliers) and higher efficiency in heavy-tailed
distributions.

Alternative tools, which are a lot faster, are
L<C<CollectInsertSizeMetrics>|https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics>
from L<Picard Tools|https://broadinstitute.github.io/picard/> and
L<C<sam-stats>|https://code.google.com/p/ea-utils/wiki/SamStats> from
L<ea-utils|https://code.google.com/p/ea-utils/>.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input SAM file or piped C<STDIN> (-) from a BAM file e.g. with
L<C<samtools view>|http://www.htslib.org/doc/samtools-1.1.html> from
L<Samtools|http://www.htslib.org/>

=item B<-a>, B<-align>

B<Default method:> Align method to calculate insert size statistics,
includes only reads which are mapped in a proper/concordant pair (as
determined by the mapper). Excludes option B<-p>.

B<or>

=item B<-p>, B<-percentile>

Percentile method to calculate insert size statistics, includes only
read pairs with an insert size within the 10th and the 90th
percentile range of all mapped read pairs. However, the frequency
distribution as well as the histogram will be plotted with the 'raw'
insert size data before percentile filtering. Excludes option B<-a>.

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-d>, B<-distro>

Create distribution histograms for the insert sizes and read lengths
with L<R|http://www.r-project.org/>. The calculated median and mean
(that are printed to C<STDOUT>) are plotted as vertical lines into
the histograms. Use it to control the correctness of the statistical
calculations.

=item B<-f>, B<-frequencies>

Print the frequencies of the insert sizes and read lengths to
tab-delimited files 'ins_frequency.txt' and 'read_frequency.txt',
respectively.

=item B<-max>=I<int>, B<-max_ins_cutoff>=I<int>

Set a maximal insert size cutoff, all insert sizes above this cutoff
will be discarded (doesn't affect read length). With B<-min> and
B<-max> you can basically run both methods, by first running the
script with B<-p> and then using the 10th and 90th percentile of the
'raw data' as B<-min> and B<-max> for option B<-a>.

=item B<-min>=I<int>, B<-min_ins_cutoff>=I<int>

Set a minimal insert size cutoff [default = 1]

=item B<-n>=I<int>, B<-num_read>=I<int>

Number of reads to sample for the calculations from the start of the
SAM/BAM file. Significant statistics can usually be calculated from a
fraction of the total SAM/BAM alignment file.

=item B<-xlim_i>=I<int>, B<-xlim_ins>=I<int>

Set an upper limit for the x-axis of the B<'R'> B<insert size>
histogram, overriding automatic truncation of the histogram tail.
The default cutoff is one and a half times the third quartile Q3
(75th percentile) value. The minimal cutoff is set to the lowest
insert size automatically. Forces option B<-d>.

=item B<-xlim_r>=I<int>, B<-xlim_read>=I<int>

Set an upper limit for the x-axis of the optional B<'R'> B<read
length> histogram. Default value is as in B<-xlim_i>. Forces option
B<-d>.

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

Calculated stats are printed to C<STDOUT>

=item F<./results>

All optional output files are stored in this results folder

=item (F<./results/ins_frequency.txt>)

Frequencies of insert size 'raw data', tab-delimited

=item (F<./results/ins_histo.pdf>)

Distribution histogram for the insert size 'raw data'

=item (F<./results/read_frequency.txt>)

Frequencies of read lengths, tab-delimited

=item (F<./results/read_histo.pdf>)

Distribution histogram for the read lengths. Not informative if
there's no variation in the read lengths.

=back

=head1 EXAMPLES

=over

=item C<perl sam_insert-size.pl -i file.sam -a -d -f>

=item C<samtools view -h file.bam | perl sam_insert-size.pl -i - -p
-max 500 -n 2000000 -xlim_i 350>

=back

=head1 DEPENDENCIES

=over

=item B<Statistics::Descriptive>

Perl module to calculate descriptive statistics, if not installed
already get it from L<CPAN|http://www.cpan.org/>

=item B<Statistical computing language L<R|http://www.r-project.org/>>

C<Rscript> is needed to plot the histograms with option B<-d>

=back

=head1 VERSION

 0.2                                               update: 29-10-2014
 0.1                                                       27-11-2013

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

=head1 ACKNOWLEDGEMENTS

References/thanks go to:

- Tobias Rausch's online courses/workshops (EMBL Heidelberg) on the
introduction to SAM files and flags L<http://www.embl.de/~rausch/>

- The CBS NGS Analysis course for the percentile filtering idea:
L<http://www.cbs.dtu.dk/courses/27626/programme.php>

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
eval { # check if module is installed
    require Statistics::Descriptive; # module of basic descriptive statistics
    Statistics::Descriptive->import();
}; # semi-colon needed
if ($@) { # if module not installed die with error message
    die "### Fatal error: Perl module 'Statistics::Descriptive' is not installed, but required for this script! Please install from CPAN!\n\n";
}


### Get the options with Getopt::Long
my $Input_Sam; # input SAM file
my $Opt_Align; # option '-a' to include only properly/concordantly mapped read pairs for insert size stats, excludes '-p'
my $Opt_Percentile; # option '-p' to calculate insert size stats within the 10th percentile (P10) and P90, excludes '-a'
my $Max_Ins_Cutoff; # maximum insert size cutoff
my $Min_Ins_Cutoff = 1; # minimum insert size cutoff
my $Num_Read; # number of reads to sample
my $Opt_Freq; # print insert size/read length counts to files 'ins_frequency.txt' and 'read_frequency.txt'
my $Opt_Histo; # draw distribution histograms with 'Rscript'
my $Ins_Xlim; # set x-axis upper limit for the insert size R histogram
my $Read_Xlim; # set x-axis upper limit for the read length R histogram
my $VERSION = 0.2;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_Sam,
            'align' => \$Opt_Align,
            'percentile' => \$Opt_Percentile,
            'max_ins_cutoff=i' => \$Max_Ins_Cutoff,
            'min_ins_cutoff=i' => \$Min_Ins_Cutoff,
            'num_read=i' => \$Num_Read,
            'frequencies' => \$Opt_Freq,
            'distro' => \$Opt_Histo,
            'xlim_ins=i' => \$Ins_Xlim,
            'xlim_read=i' => \$Read_Xlim,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_Sam) {
    my $warning = "\n### Fatal error: Option '-i' or its argument is missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Enforce mandatory or optional options
if (!$Opt_Align && !$Opt_Percentile) {
    warn "### Warning: None of the mandatory options ('-a' or '-p') given. Forcing default option '-a'!\n";
    $Opt_Align = 1;
} elsif ($Opt_Align && $Opt_Percentile) {
    die "\n### Fatal error: Both mandatory options ('-a' and '-p') given! Choose only one of the filter methods!\n";
}
if (($Ins_Xlim || $Read_Xlim) && !$Opt_Histo) { # force option '-d' to create the histograms for $Ins_Xlim and $Read_Xlim
    warn "### Warning: One 'xlim'-option given, but not option '-d' to create histograms. Forcing option '-d'!\n";
    $Opt_Histo = 1;
}



### Open SAM file or accept STDIN
my $Sam_Fh;
if ($Input_Sam eq '-') { # file input via STDIN
    $Sam_Fh = *STDIN; # capture typeglob of STDIN
} else { # input via SAM file
    open ($Sam_Fh, "<", "$Input_Sam");
}



### Variables for basic alignment stats
my ($Ref_Name, $Ref_Size), # name of the reference and ref sequence length
my $Ref_Count = 0; # number of ref sequences (e.g. if a draft genome is used)
my $Total_Reads = 0; # total reads in SAM file
my @Mapq; # Phred/Sanger-based mapping quality
my $Qc_Fail = 0; # reads failing platform/vendor quality checks
my $Unmapped_Reads = 0; # reads with no match in the alignment
my $Second_Aln_Reads = 0; # secondary alignments for multiple mapping reads (typically best alignment is primary for that read)
my $Supp_Aln_Reads = 0; # supplementary alignments for chimeric reads (only one partial alignment is primary for that read)
my $Paired_Reads = 0; # reads paired in sequencing
my $Mapped_Pair_Reads = 0; # paired reads, both reads (read and mate) mapped
my $Proper_Pair_Reads = 0; # reads mapped in a proper pair



### Read in SAM/BAM file, store insert sizes, read lengths, and counts for basic alignment stats
my %Insert_Size_Counts; # store insert size counts
my @Ins_Sizes; # ALL insert sizes for stats and histograms
my %Read_Length_Counts; # store read length counts
my @Read_Lengths; # ALL read lengths
my $I = 1; # mio. counter for processing status message
while (<$Sam_Fh>) {
    chomp;
    next if (/^\s*$/); # skip empty lines
    if (/^\@SQ/) { # parse ref name and ref size out of SAM header
        $Ref_Count++;
        if ($Ref_Count == 1) { # only keep name and size of first reference
            $Ref_Name = $1 if (/SN\:(\S+)/); # ref name
            $Ref_Size = $1 if (/LN\:(\d+)/); # ref seq length
        }
    }
    next if (/^@/); # skip remaining SAM header lines

    my @fields = split(/\t/, $_); # split tab-separated SAM lines
    # $fields[1]= 16-bit FLAG, [4]= Phred/Sanger-based mapping quality, [5]= CIGAR string, [8]= insert size
    # FLAG bits:
    # Bits      Decimal     Hexadecimal    Description
    # 2^0       1           0x0001         read paired in sequencing, independent of pair-mapping
    # 2^1       2           0x0002         read mapped in proper pair
    # 2^2       4           0x0004         read unmapped
    # 2^3       8           0x0008         mate is unmapped
    # 2^4       16          0x0010         strand of read mapping (0 for fw; 1 for rev)
    # 2^5       32          0x0020         strand of mate mapping
    # 2^6       64          0x0040         read first in pair
    # 2^7       128         0x0080         read second in pair
    # 2^8       256         0x0100         alignment not primary (read maps several times on reference)
    # 2^9       512         0x0200         read fails plaftorm/vendor quality checks
    # 2^10      1024        0x0400         read either PCR or optical duplicate
    # 2^11      2048        0x0800         supplementary alignment of chimeric read
    die "\n### Fatal error: Less than 11 tab-separated fields in SAM/BAM file! Sure this is a SAM/BAM file?\n" if (@fields < 11); # SAM format has 11 mandatory tab-separated fields per line

    # discard secondary and supplementary alignments (has to be before total read count)
    $Second_Aln_Reads++ if ($fields[1] & 0x100); # count secondary alignments
    if ($fields[1] & 0x800) { # count and skip supplementary alignments
        $Supp_Aln_Reads++;
        next;
    }
    next if ($fields[1] & 0x100); # skip secondary alignments

    # status message to STDERR with number of reads processed
    if ($Total_Reads/1000000 == $I) {
        print STDERR "$I Mio reads processed ...\r"; # carriage return to overwrite messages and not clutter STDERR
        $I++;
    }
    last if ($Num_Read && ($Total_Reads == $Num_Read)); # skip rest of reads if $Num_Read is reached
    $Total_Reads++;


    # read length counts
    $Qc_Fail++ if ($fields[1] & 0x200); # count reads that didn't pass platform/vendor QC checks
    if ($fields[1] & 0x4) { # count reads not mapped on the reference
        $Unmapped_Reads++;
    } else { # if read is mapped include in read length stats
        push(@Mapq, $fields[4]); # store mapping quality of mapped read
        $Read_Length_Counts{read_len($fields[5])}++; # subroutine to calculate the read length from the CIGAR string
        push(@Read_Lengths, read_len($fields[5])); # subroutine
    }
    $fields[8] = abs $fields[8]; # absolute value of insert size (script only stores info from first read in a pair, see below, which can be on either strand)

    # mapped read pairs
    if ($fields[1] & 0x1) { # read paired in sequencing
        $Paired_Reads++;
        if (!($fields[1] & 0x4) && !($fields[1] & 0x8)) { # pair (read and mate) mapped
            $Mapped_Pair_Reads++;
            if ($Opt_Percentile && ($fields[8] >= $Min_Ins_Cutoff && (!$Max_Ins_Cutoff || $fields[8] <= $Max_Ins_Cutoff))) { # option '-p', only include reads inside the '-min' (default 1) and '-max' cutoffs
                if ($fields[1] & 0x40) { # only include first read in pair, otherwise insert sizes will be double counted
                    $Insert_Size_Counts{$fields[8]}++;
                    push(@Ins_Sizes, $fields[8]);
                }
            }
        }
    }

    # reads mapped in proper pair
    if ($fields[1] & 0x2) {
        $Proper_Pair_Reads++;
        if ($Opt_Align && ($fields[8] >= $Min_Ins_Cutoff && (!$Max_Ins_Cutoff || $fields[8] <= $Max_Ins_Cutoff))) { # option '-a'
            if ($fields[1] & 0x40) { # first read in pair
                $Insert_Size_Counts{$fields[8]}++;
                push(@Ins_Sizes, $fields[8]);
            }
        }
    }
}
warn "\n"; # get rid of the carriage return for status messages



### Basic alignment stats
if ($Ref_Name && $Ref_Size) { # print only if reference info given in the SAM/BAM header
    print "\nStatistics for the short-read mapping on reference '$Ref_Name' ($Ref_Size bp)";
    if ($Ref_Count > 1) { # if several references, print only first name
        print " and ", $Ref_Count-1, " additional references:\n";
    } else {
        print ":\n";
    }
}
print "Overall read statistics:\n";
print "\tDiscarded secondary alignments of multiple mapping reads: $Second_Aln_Reads\n";
print "\tDiscarded supplementary alignments of chimeric reads: $Supp_Aln_Reads\n";
print "\tTotal read count: $Total_Reads\n";
print "\tReads failing platform/vendor quality checks: $Qc_Fail ("; perc($Qc_Fail, $Total_Reads); print "%)\n"; # subroutine to print the percentage with printf
print "\tReads paired in sequencing: $Paired_Reads ("; perc($Paired_Reads, $Total_Reads); print "%)\n"; # subroutine
print "\tReads mapped on reference: ", scalar @Read_Lengths, " ("; perc(scalar @Read_Lengths, $Total_Reads); print "%)\n"; # subroutine
print "\tUnmapped reads: $Unmapped_Reads ("; perc($Unmapped_Reads, $Total_Reads); print "%)\n"; # subroutine
print "\tMean Phred/Sanger-based mapping quality: ";
my %Mapq_Stats; # store stats for mapping qualities
stats_full(\@Mapq, \%Mapq_Stats); # subroutine to calculate the stats for the data in the array with 'Statistics::Descriptive' and store result stats in the given hash
print "$Mapq_Stats{'mean'}\n";
print "\tPaired reads mapped on reference ('raw data' used for option '-p'): $Mapped_Pair_Reads ("; perc($Mapped_Pair_Reads, $Total_Reads); print "%)\n"; # subroutine
print "\tReads mapped in a proper pair (used for option '-a'): $Proper_Pair_Reads ("; perc($Proper_Pair_Reads, $Total_Reads); print "%)\n"; # subroutine



### Read lengths stats
my %Read_Length_Stats; # store stats for read lengths
if (@Read_Lengths > 0) {
    stats_full(\@Read_Lengths, \%Read_Length_Stats); # subroutine for 'Statistics::Descriptive'
    print "\nRead length statistics of all mapped reads:\n";
    print "\tRead length min value: $Read_Length_Stats{'min'}\tmax value: $Read_Length_Stats{'max'}\n";
    print "\tRead length quantiles, Q1 (25th percentile): $Read_Length_Stats{'q1'}\tQ3 (75th percentile): $Read_Length_Stats{'q3'}\n";
    print "\tRead length median: $Read_Length_Stats{'median'}\n";
    print "\tRead length mean: $Read_Length_Stats{'mean'}\n";
    print "\tRead length standard deviation: $Read_Length_Stats{'sd'}\n";
} else {
    die "\n### Fatal error: No CIGAR strings found in the SAM file, sure this is a SAM file from a read mapping?\n";
}



### Insert sizes stats
my %Ins_Size_Stats; # store stats for insert sizes
if (@Ins_Sizes > 0) {
    print "\nInsert size statistics of mapped reads ";
    if ($Opt_Percentile) { # for option '-p' filter insert sizes by the 10th and 90th percentile of the original 'raw data'
        print "with option '-p':\n";
        print "\t\"Raw data\" insert sizes before percentile filtering (pairs mapped on reference divided by two, excluding insert sizes of zero): ", scalar @Ins_Sizes, "\n";
        print "\t\"Raw data\" statistics for subsequent percentile filtering: ";
        stats_full(\@Ins_Sizes, \%Ins_Size_Stats); # subroutine to calculate stats with the "raw data" to subsequently filter by 10th and 90th percentile with grep
        print "10th_percentile=$Ins_Size_Stats{'p10'} 90th_percentile=$Ins_Size_Stats{'p90'} min=$Ins_Size_Stats{'min'} max=$Ins_Size_Stats{'max'} Q1=$Ins_Size_Stats{'q1'} median=$Ins_Size_Stats{'median'} Q3=$Ins_Size_Stats{'q3'} mean=$Ins_Size_Stats{'mean'} sd=$Ins_Size_Stats{'sd'}\n";
        my @filter_ins_size; # new array with filtered insert sizes
        @filter_ins_size = grep ($_ >= $Ins_Size_Stats{'p10'} && $_ <= $Ins_Size_Stats{'p90'}, @Ins_Sizes); # grep elements according to calculated 10th and 90th percentile; retain 'raw data' @Ins_Sizes for frequency file and histogram
        print "\tInsert sizes used to calculate statistics (after percentile filtering): ", scalar @filter_ins_size, "\n";
        stats_full(\@filter_ins_size, \%Ins_Size_Stats); # subroutine for 'Statistics::Descriptive'
    } elsif ($Opt_Align) { # option '-a'
        print "with option '-a':\n";
        print "\tInsert sizes used to calculate statistics (reads mapped in proper pair divided by two, excluding insert sizes of zero): ", scalar @Ins_Sizes, "\n";
        stats_full(\@Ins_Sizes, \%Ins_Size_Stats); # subroutine
    }
    print "\tInsert size min value: $Ins_Size_Stats{'min'}\tmax value: $Ins_Size_Stats{'max'}\n";
    print "\tInsert size quantiles, Q1: $Ins_Size_Stats{'q1'}\tQ3: $Ins_Size_Stats{'q3'}\n";
    print "\tInsert size median: $Ins_Size_Stats{'median'}\n";
    print "\tInsert size mean: $Ins_Size_Stats{'mean'}\n";
    print "\tInsert size standard deviation: $Ins_Size_Stats{'sd'}\n\n";
} else {
    warn "\n### Warning: No insert sizes detected in the SAM file, sure this is paired-end data?\n\n";
}



### Weighted stats; only needed if arrays @Read_Lengths and @Ins_Sizes get too big
### then the distribution hashes more useful as they are smaller
### CPAN modules 'Statistics::Descriptive::Discrete' and 'Statistics::Descriptive::Weighted' can be used then (see 'calc_fastq-stats.pl')
#my $Size_Weighted = 0;
#my $Elements = 0; # insert size count
## weighted mean
#foreach my $size (keys %Insert_Size_Counts) {
    #$Size_Weighted = $Size_Weighted + ($size*$Insert_Size_Counts{$size});
    #$Elements += $Insert_Size_Counts{$size};
#}
#my $Mean_Weighted = $Size_Weighted/$Elements;
#print "\tWeighted mean: ", $Mean_Weighted, "\n";
#$Size_Weighted = 0; # set back to zero to calculate weighted variance
## weighted stdev
#foreach my $size (keys %Insert_Size_Counts) {
    #$Size_Weighted += ($size - $Mean_Weighted)*($size - $Mean_Weighted)*$Insert_Size_Counts{$size}; # variance
#}
#print "\tWeighted stdev: ", sqrt($Size_Weighted/$Elements), "\n"; # square root for stdev



### Optionally, create results directory for output files
my $Results_Dir = './results';
if($Opt_Freq || $Opt_Histo) {
    mkdir $Results_Dir if (!-e $Results_Dir);
}



### Optionally, create insert size and read length distribution files
if ($Opt_Freq) {
    warn "The insert size and read length distribution files are created in the '$Results_Dir' directory:\n";

    my $ins_freq_out = 'ins_frequency.txt';
    file_exist("$Results_Dir/$ins_freq_out");
    open (my $ins_freq_fh, ">", "$Results_Dir"."/$ins_freq_out");
    warn "\t$ins_freq_out\n";
    print $ins_freq_fh "insert_size\tcount\n"; # header for the file
    foreach my $size (sort {$a <=> $b} keys %Insert_Size_Counts) {
        print $ins_freq_fh "$size\t$Insert_Size_Counts{$size}\n";
    }

    my $read_freq_out = 'read_frequency.txt';
    file_exist("$Results_Dir/$read_freq_out");
    open (my $read_freq_fh, ">", "$Results_Dir"."/$read_freq_out");
    warn "\t$read_freq_out\n\n";
    print $read_freq_fh "read_length\tcount\n"; # header
    foreach my $length (sort {$a <=> $b} keys %Read_Length_Counts) {
        print $read_freq_fh "$length\t$Read_Length_Counts{$length}\n";
    }
    close $ins_freq_fh;
    close $read_freq_fh;
}



### Optionally, create histogram pdfs for the distribution of insert size and read length
if ($Opt_Histo) {
    warn "The following histogram pdfs are created in the '$Results_Dir' directory:\n";
    r_histo('ins', \@Ins_Sizes, \%Ins_Size_Stats, $Ins_Xlim); # subroutine to create the insert size histogram with 'Rscript' of 'R'
    r_histo('read', \@Read_Lengths, \%Read_Length_Stats, $Read_Xlim); # subroutine for read length histogram
    warn "\n";
}


close $Sam_Fh; # closing this filehandle to early will result in a warning message for STDIN input?!
exit;


###############
# Subroutines #
###############

### Subroutine to test for file existence and give warning to STDERR
sub file_exist {
    my $file = shift;
    if (-e $file) {
        warn "\nThe result file '$file' exists already and will be overwritten!\n\n";
        return 1;
    }
    return 0;
}



### Subroutine to print the precentage rounded to two decimal places
sub perc {
    my ($numerator, $denominator) = @_;
    printf("%.2f", ($numerator/$denominator)*100);
    return 1;
}



### Subroutine to plot a histogram by creating an R script and executing it with 'Rscript'
sub r_histo {
    my ($prefix, $data_array_ref, $hash_stat_ref, $xlim) = @_; # prefix is either 'ins' for insert size or 'read' for read length
    $xlim = 1.5 * $hash_stat_ref->{'q3'} if (!$xlim); # set xlim upper limit to one and a half times Q3 if not given as option

    # label for the histogram
    my $label;
    if ($prefix eq 'ins') {
        $label = 'Insert size';
    } elsif ($prefix eq 'read') {
        $label = 'Read length';
    }

    # create temporary input file for R with input sizes/read lengths
    my $tmp_file = "tmp_"."$prefix".".txt";
    open (my $tmp_file_fh, ">", "$tmp_file");
    foreach (@$data_array_ref) { # de-reference array-ref for iteration
        print $tmp_file_fh "$_\n" if ($_ <= $xlim); # only insert sizes/read lengths below xlim needed
    }
    close $tmp_file_fh;

    # create R script
    my $tmp_r_script = "tmp_"."$prefix"."_hist.r"; # filename of the R script
    my $histo_name = "$prefix"."_histo.pdf"; # filename of the output pdf histogram
    file_exist("$Results_Dir/$histo_name");
    open (my $r_fh, ">", "$tmp_r_script");

    select $r_fh; # select fh for standard print/f output
    print "#!/usr/bin/Rscript --vanilla --slave\n"; # header of R script
    print "$prefix = scan(file=\"$tmp_file\", quiet=T)\n"; # scan in data and suppress sdtout output
    print "pdf(\"$Results_Dir/$histo_name\")\n"; # filename of the pdf
    print "hist($prefix, breaks=$xlim, xlim=c(min($prefix),$xlim), main=NULL, xlab=\"$label [bp]\", ylab=\"$label count\")\n"; # create the histogram with labels, but without MAIN-title (see below)
    if ($Ref_Name && $Ref_Count == 1) { # title for histogram plot
        print "title(\"$label distribution for mapping on\\n$Ref_Name\", adj=0)\n"; # align title left to have room for mtext below
    } else {
        print "title(\"$label distribution\")\n";
    }
    print "abline(v=$hash_stat_ref->{'median'}, col=\"blue\", lwd=1)\n"; # plot the calculated median into the histogram
    print "abline(v=$hash_stat_ref->{'mean'}, col=\"green\", lwd=1)\n";  # mean
    print "legend(\"topright\", c(\"median\", \"mean\"), cex=0.8, col=c(\"blue\", \"green\"), lwd=1)\n";
    print "mtext(\"median=$hash_stat_ref->{'median'}\\nmean=$hash_stat_ref->{'mean'}\\nstdev=$hash_stat_ref->{'sd'}\", side=3, adj=1, cex=0.8, col=\"red\")\n"; # add calculated stats to the plot margin
    print "out <- dev.off()\n"; # write histogram to the pdf and suppress STDOUT output by diverting it
    close $r_fh;
    select STDOUT; # change standard print/f output to STDOUT

    # execute R script with Rscript
    system("Rscript $tmp_r_script") == 0 or die "### Fatal error: Statistical programming language 'R' either not installed, not in \$PATH, or something wrong with '$tmp_r_script'! Install 'R' to create the histograms (required is 'Rscript')!\n";
    warn "\t$histo_name\n"; # print to STDERR which file has been created
    unlink ($tmp_file, $tmp_r_script); # remove the tmp files
    return 1;
}



### Subroutine to get the read length from the CIGAR string
sub read_len {
    my $cigar = shift;
    $cigar =~ s/\d+(D|N|H|P)//g; # SAM specs: sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ, hence delete all other CIGAR operations with global modifier
    my @cigar_split = split(/\D/, $cigar); # split CIGAR string via non-digital characters
    my $read_length = 0;
    foreach (@cigar_split) { # sum up the CIGAR operations
        $read_length += $_;
    }
    return $read_length;
}



### Subroutine to calculate stats with module 'Statistics::Descriptive'
sub stats_full {
    my ($data_array_ref, $hash_stat_ref) = @_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data_array_ref); # de-reference array ref
    $hash_stat_ref->{'mean'} = sprintf("%.2f", $stat->mean()); # rounded to two decimal places
    $hash_stat_ref->{'median'} = $stat->median();
    $hash_stat_ref->{'sd'} = sprintf("%.2f", $stat->standard_deviation());
    $hash_stat_ref->{'min'} = $stat->min();
    $hash_stat_ref->{'max'} = $stat->max();
    $hash_stat_ref->{'q1'} = $stat->quantile(1); # Q1, first quartile (25th percentile)
    $hash_stat_ref->{'q3'} = $stat->quantile(3); # Q3, third quartile (75th percentile)
    my ($percentile, $index) = $stat->percentile(10); # 10th percentile (in list context also returns the index of the percentile, which is not needed)
    $hash_stat_ref->{'p10'} = $percentile;
    ($percentile, $index) = $stat->percentile(90); # 90th percentile
    $hash_stat_ref->{'p90'} = $percentile;
    # $stat->clear(); # remove stats from the module for next round; not needed because of 'new()' initialisation at begin of subroutine?!
    return 1;
}
