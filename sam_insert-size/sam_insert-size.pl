#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

sam_insert-size.pl                                          27-11-2013

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
(for SAM specifications see: L<http://samtools.sourceforge.net/SAM1.pdf>).

B<Read length> statistics are calculated for all mapped reads
(irrespective of their pairing).

B<Insert size> statistics are calculated only for B<paired reads>.
Typically, the insert size is perturbed by artifacts, like chimeras,
structural re-arrangements or alignment errors, which result in a very
high maximum insert size measure. As a consequence the mean and stdev
can be strongly misleading regarding the real distribution. To avoid
this, two methods are implemented that first trim the insert
size distribution to a 'core' to calculate the respective statistics.
Additionally, secondary alignments for multiple mapping reads and
supplementary alignments for chimeric reads, as well as insert sizes
of zero are not considered (option '-min_ins_cutoff' is set to
B<one> by default).

The B<'-a|-align'> method includes only proper/concordant paired reads
in the statistical calculations (as determined by the mapper and the
options for insert size minimum and maximum used for mapping). This is
the B<default> method.

The B<'-p|-percentile'> method first calculates insert size statistics
for all read pairs, where the read and the mate are mapped ('raw
data'). Subsequently, the 10th and the 90th percentile are discarded
to calculate the 10% truncated mean and stdev. Discarding the lowest
and highest 10% of insert sizes gives the advantage of robustness
(insensitivity to outliers) and higher efficiency in heavy-tailed
distributions.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input SAM file or piped STDIN (-) from a BAM file with C<samtools
view>.

=item B<-a>, B<-align>

B<Default method:> Align method to calculate insert size statistics,
includes only reads which are mapped in a proper/concordant pair (as
determined by the mapper). If option '-p' is not given, '-a' will be
forced as default option.

=item B<-p>, B<-percentile>

Percentile method to calculate insert size statistics, includes only
read pairs with an insert size within the 10th and the 90th percentile
range of all mapped read pairs. However, the frequency distribution as
well as the histogram will be plotted with the 'raw' insert size data
before percentile filtering. Excludes option '-a'.

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD).

=item B<-d>, B<-distro>

Create distribution histograms for the insert sizes and read lengths
with B<'R'>. The calculated median and mean (that are printed to stdout)
are plotted as vertical lines into the histograms. Use it to control
the correctness of the statistical calculations.

=item B<-f>, B<-frequencies>

Print the frequencies of the insert sizes and read lengths to the
tab-delimited files 'ins_frequency.txt' and 'read_frequency.txt',
respectively.

=item B<-max>=I<int>, B<-max_ins_cutoff>=I<int>

Set a maximal insert size cutoff, all insert sizes above this cutoff
will be discarded (doesn't affect read length). With '-min' and '-max'
you can basically run both methods, by first running the script with
B<'-p'> and then using the 10th and 90th percentile of the 'raw data'
as '-min' and '-max' for option B<'-a'>.

=item B<-min>=I<int>, B<-min_ins_cutoff>=I<int>

Set a minimal insert size cutoff [default = 1].

=item B<-n>=I<int>, B<-num_read>=I<int>

Number of reads to sample for the calculations from the start of the
SAM/BAM file. Significant statistics can usually be calculated from a
fraction of the total SAM/BAM alignment file.

=item B<-xlim_i>=I<int>, B<-xlim_ins>=I<int>

Set an upper limit for the x-axis of the B<'R'> B<insert size> histogram,
overriding automatic truncation of the histogram tail. The default
cutoff is one and a half times the third quartile Q3 (75th percentile)
value. The minimal cutoff is set to the lowest insert size
automatically. Forces option B<'-d'>.

=item B<-xlim_r>=I<int>, B<-xlim_read>=I<int>

Set an upper limit for the x-axis of the optional B<'R'> B<read length>
histogram. Default value is as in B<'-xlim_i>'. Forces option B<'-d'>.

=item B<-v>, B<-version>

Print version number to STDERR.

=back

=head1 OUTPUT

=over 20

=item F<./results>

All optional output files are stored in this results folder. Already
existing output files will be overwritten without warning!

=item (F<ins_frequency.txt>)

Frequencies of insert size 'raw data', tab-delimited.

=item (F<ins_histo.pdf>)

Distribution histogram for the insert size 'raw data'.

=item (F<read_frequency.txt>)

Frequencies of read lengths, tab-delimited.

=item (F<read_histo.pdf>)

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

The Perl module is needed to calculate statistics, if not already in
the Perl core packages install from CPAN L<http://search.cpan.org/>.

=item B<Statistical computing language I<R>>

B<Rscript> is needed to plot the histograms with option B<'-d'>.

=back

=head1 VERSION

0.1

=head1 CHANGELOG

=over

=item C<0.x>

Nothing yet ;-).

=back

=head1 ACKNOWLEDGEMENTS

References/thanks go to:

- Tobias Rauschs, EMBL Heidelberg, online courses/workshops on
the introduction to SAM files and flags:
http://www.embl.de/~rausch/

- The CBS NGS Analysis course for the percentile filtering idea:
http://www.cbs.dtu.dk/courses/27626/programme.php

=head1 AUTHOR

A Leimbach                                     aleimba[at]gmx[dot]de

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
eval { # check if module is installed
    require Statistics::Descriptive; # module of basic descriptive statistics
    Statistics::Descriptive->import();
}; # semi-colon needed
if ($@) { # if module not installed die with error message
    die "### Fatal error: Perl module \'Statistics::Descriptive\' is not installed, but required for this script! Please install from CPAN!\n\n";
}


### Get the options with Getopt::Long
my $input_sam; # input SAM file
my $opt_align; # option '-a', excludes '-p'
my $opt_percentile; # option '-p', excludes '-a'
my $max_ins_cutoff; # maximum insert size cutoff
my $min_ins_cutoff = 1; # minimum insert size cutoff
my $num_read; # option '-n'
my $opt_freq; # option '-f'
my $opt_histo; # option '-d'
my $ins_xlim; # set x-axis upper limit for the insert size R histogram
my $read_xlim; # set x-axis upper limit for the read length R histogram
my $version = 0.1;
my ($opt_version, $opt_help);
GetOptions ('input=s' => \$input_sam,
            'align' => \$opt_align,
            'percentile' => \$opt_percentile,
            'max_ins_cutoff:i' => \$max_ins_cutoff,
            'min_ins_cutoff:i' => \$min_ins_cutoff,
            'num_read:i' => \$num_read,
            'frequencies' => \$opt_freq,
            'distro' => \$opt_histo,
            'xlim_ins:i' => \$ins_xlim,
            'xlim_read:i' => \$read_xlim,
            'version' => \$opt_version,
            'help|?' => \$opt_help);


### Run perldoc on POD
pod2usage(-verbose => 2) if ($opt_help);
if ($opt_version) {
    die "$0 $version\n";
} elsif (!$input_sam) {
    my $warning = "\n### Fatal error: Option \'-i\' or its argument are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}


### Enforce mandatory or optional options
if (!$opt_align && !$opt_percentile) {
    warn "### Warning: None of the mandatory options (\'-a\' or \'-p\') given. Forcing default option \'-a\'!\n";
    $opt_align = 1;
} elsif ($opt_align && $opt_percentile) {
    die "\n### Fatal error: Both mandatory options (\'-a\' and \'-p\') given! Choose only one of the filter methods!\n";
}
if (($ins_xlim || $read_xlim) && !$opt_histo) { # force option '-d' to create the histograms for $ins_xlim and $read_xlim
    warn "### Warning: One 'xlim'-option given, but not option \'-d\' to create histograms. Forcing option \'-d\'!\n";
    $opt_histo = 1;
}


### Open SAM file or accept STDIN
my $sam_fh; # filehandle for SAM file
if ($input_sam eq '-') { # file input via STDIN
    $sam_fh = *STDIN; # capture typeglob of STDIN
} else { # input via SAM file
    open ($sam_fh, "<", "$input_sam");
}


### Variables for basic alignment stats
my ($ref_name, $ref_size), # name of the reference and ref sequence length
my $ref_count = 0; # number of ref sequences (e.g. if a draft genome is used)
my $total_reads = 0; # total reads in SAM file
my @mapq; # phred-based mapping quality
my $qc_fail = 0; # reads failing platform/vendor quality checks
my $unmapped_reads = 0; # reads with no match in the alignment
my $second_aln_reads = 0; # secondary alignments for multiple mapping reads (typically best alignment is primary for that read)
my $supp_aln_reads = 0; # supplementary alignments for chimeric reads (only one partial alignment is primary for that read)
my $paired_reads = 0; # reads paired in sequencing
my $mapped_pair_reads = 0; # paired reads, both reads (read and mate) mapped
my $proper_pair_reads = 0; # reads mapped in a proper pair


### Read in SAM/BAM file, store insert sizes, read lengths, and counts for basic alignment stats
my %insert_sizes; # store insert size distribution
my @ins_size_array; # ALL insert sizes for stats and histograms
my %read_lengths; # store read length distribution
my @read_length_array; # ALL read lengths
my $i = 1; # mio. counter for processing status message
while (<$sam_fh>) {
    chomp;
    next if (/^\s*$/); # skip empty lines
    if (/^\@SQ/) { # parse ref name and ref size out of SAM header
        $ref_count++;
        if ($ref_count == 1) { # only keep name and size of first reference
            $ref_name = $1 if (/SN\:(\S+)/); # ref name
            $ref_size = $1 if (/LN\:(\d+)/); # ref seq length
        }
    }
    next if (/^@/); # skip remaining SAM header lines

    my @fields = split(/\t/, $_); # split tab-separated SAM lines
    # $fields[1]= 16-bit FLAG, [4]= Phred-based mapping quality, [5]= CIGAR string, [8]= insert size
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
    $second_aln_reads++ if ($fields[1] & 0x100); # count secondary alignments
    if ($fields[1] & 0x800) { # count and skip supplementary alignments
        $supp_aln_reads++;
        next;
    }
    next if ($fields[1] & 0x100); # skip secondary alignments

    # count total reads
    last if (defined $num_read && ($total_reads == $num_read)); # skip rest of reads if $num_read is reached
    $total_reads++;

    # status message to stderr with number of reads processed
    if ($total_reads/1000000 == $i) {
        print STDERR "$i Mio reads processed ...\r"; # carriage return to overwrite messages and not clutter stderr
        $i++;
    }

    # read length counts
    $qc_fail++ if ($fields[1] & 0x200); # count reads that didn't pass platform/vendor QC checks
    if ($fields[1] & 0x4) { # count reads not mapped on the reference
        $unmapped_reads++;
    } else { # if read is mapped include in read length stats
        push(@mapq, $fields[4]); # store mapping quality of mapped read
        $read_lengths{read_len($fields[5])}++; # subroutine to calculate the read length from the CIGAR string
        push(@read_length_array, read_len($fields[5])); # subroutine
    }
    $fields[8] = abs $fields[8]; # absolute value of insert size (script only stores info from first read in a pair, see below, which can be on either strand)

    # mapped read pairs
    if ($fields[1] & 0x1) { # read paired in sequencing
        $paired_reads++;
        if (!($fields[1] & 0x4) && !($fields[1] & 0x8)) { # pair (read and mate) mapped
            $mapped_pair_reads++;
            if ($fields[8] >= $min_ins_cutoff || (defined $max_ins_cutoff && $fields[8] <= $max_ins_cutoff)) { # only include reads inside the '-min' (default 1) and '-max' cutoffs
                if ($opt_percentile) { # option '-p'
                    if ($fields[1] & 0x40) { # only include first read in pair, otherwise insert sizes will be double counted
                        $insert_sizes{$fields[8]}++;
                        push(@ins_size_array, $fields[8]);
                    }
                }
            }
        }
    }

    # reads mapped in proper pair
    if ($fields[1] & 0x2) {
        $proper_pair_reads++;
        if ($fields[8] >= $min_ins_cutoff || (defined $max_ins_cutoff && $fields[8] <= $max_ins_cutoff)) {
            if ($opt_align) { # option '-a'
                if ($fields[1] & 0x40) { # first read in pair
                    $insert_sizes{$fields[8]}++;
                    push(@ins_size_array, $fields[8]);
                }
            }
        }
    }
}
warn "\n"; # get rid of the carriage return for status messages


### Basic alignment stats
if ($ref_name && $ref_size) { # print only if reference info given in the SAM/BAM header
    print "\nStatistics for the short-read mapping on reference \'$ref_name\' ($ref_size bp)";
    if ($ref_count > 1) { # if several references, print only first name
        print " and ", $ref_count-1, " additional references:\n";
    } else {
        print ":\n";
    }
}
print "Overall read statistics:\n";
print "\tDiscarded secondary alignments of multiple mapping reads: $second_aln_reads\n";
print "\tDiscarded supplementary alignments of chimeric reads: $supp_aln_reads\n";
print "\tTotal read count: $total_reads\n";
print "\tReads failing platform/vendor quality checks: $qc_fail ("; perc($qc_fail, $total_reads); print "%)\n"; # subroutine to print the percentage with printf
print "\tReads paired in sequencing: $paired_reads ("; perc($paired_reads, $total_reads); print "%)\n"; # subroutine
print "\tReads mapped on reference: ", scalar @read_length_array, " ("; perc(scalar @read_length_array, $total_reads); print "%)\n"; # subroutine
print "\tUnmapped reads: $unmapped_reads ("; perc($unmapped_reads, $total_reads); print "%)\n"; # subroutine
print "\tMean phred-based mapping quality: ";
my %mapq_stats; # store stats for mapping qualities
stats_full(\@mapq, \%mapq_stats); # subroutine to calculate the stats for the data in the array with 'Statistics::Descriptive' and store result stats in the given hash
print "$mapq_stats{'mean'}\n";
print "\tPaired reads mapped on reference (\'raw data\' used for option \'-p\'): $mapped_pair_reads ("; perc($mapped_pair_reads, $total_reads); print "%)\n"; # subroutine
print "\tReads mapped in a proper pair (used for option \'-a\'): $proper_pair_reads ("; perc($proper_pair_reads, $total_reads); print "%)\n"; # subroutine


### Read lengths stats
my %read_stats; # store stats for read lengths
if (@read_length_array > 0) {
    stats_full(\@read_length_array, \%read_stats); # subroutine for 'Statistics::Descriptive'
    print "\nRead length statistics of all mapped reads:\n";
    print "\tRead length min value: $read_stats{'min'}\tmax value: $read_stats{'max'}\n";
    print "\tRead length quantiles, Q1 (25th percentile): $read_stats{'q1'}\tQ3 (75th percentile): $read_stats{'q3'}\n";
    print "\tRead length median: $read_stats{'median'}\n";
    print "\tRead length mean: $read_stats{'mean'}\n";
    print "\tRead length standard deviation: $read_stats{'sd'}\n";
    } else {
    die "\n### Fatal error: No CIGAR strings found in the SAM file, sure this is a SAM file from a read mapping?\n";
}


### Insert sizes stats
my %ins_stats; # store stats for insert sizes
if (@ins_size_array > 0) {
    print "\nInsert size statistics of mapped reads ";
    if ($opt_percentile) { # for option '-p' filter insert sizes by the 10th and 90th percentile of the original 'raw data'
        print "with option \'-p\':\n";
        print "\t\"Raw data\" insert sizes before percentile filtering (pairs mapped on reference divided by two, excluding insert sizes of zero): ", scalar @ins_size_array, "\n";
        print "\t\"Raw data\" statistics for subsequent percentile filtering: ";
        stats_full(\@ins_size_array, \%ins_stats); # subroutine to calculate stats with the "raw data" to subsequently filter by 10th and 90th percentile with grep
        print "10th_percentile=$ins_stats{'p10'} 90th_percentile=$ins_stats{'p90'} min=$ins_stats{'min'} max=$ins_stats{'max'} Q1=$ins_stats{'q1'} Q3=$ins_stats{'q3'} mean=$ins_stats{'mean'} sd=$ins_stats{'sd'}\n";
        my @filter_ins_size; # new array with filtered insert sizes
        @filter_ins_size = grep ($_ >= $ins_stats{'p10'} && $_ <= $ins_stats{'p90'}, @ins_size_array); # grep elements according to calculated 10th and 90th percentile; retain 'raw data' @ins_size_array for frequency file and histogram
        print "\tInsert sizes used to calculate statistics (after percentile filtering): ", scalar @filter_ins_size, "\n";
        stats_full(\@filter_ins_size, \%ins_stats); # subroutine for 'Statistics::Descriptive'
    } elsif ($opt_align) { # option '-a'
        print "with option \'-a\':\n";
        print "\tInsert sizes used to calculate statistics (reads mapped in proper pair divided by two, excluding insert sizes of zero): ", scalar @ins_size_array, "\n";
        stats_full(\@ins_size_array, \%ins_stats); # subroutine
    }
    print "\tInsert size min value: $ins_stats{'min'}\tmax value: $ins_stats{'max'}\n";
    print "\tInsert size quantiles, Q1: $ins_stats{'q1'}\tQ3: $ins_stats{'q3'}\n";
    print "\tInsert size median: $ins_stats{'median'}\n";
    print "\tInsert size mean: $ins_stats{'mean'}\n";
    print "\tInsert size standard deviation: $ins_stats{'sd'}\n\n";
} else {
    warn "\n### Warning: No insert sizes detected in the SAM file, sure this is paired-end data?\n\n";
}


### Weighted stats; only needed if arrays @read_length_array and @ins_size_array get too big
### then the distribution hashes more useful as they are smaller (however percentiles etc. could not be calculated with 'Statistics::Descriptive' anymore)
# my $size_weighted = 0;
# my $elements = 0; # insert size count
## weighted mean
# foreach my $size (keys %insert_sizes) {
    # $size_weighted = $size_weighted + ($size*$insert_sizes{$size});
    # $elements = $elements + $insert_sizes{$size};
# }
# my $mean_weighted = $size_weighted/$elements;
# print "\tWeighted mean: ", $mean_weighted, "\n";
# $size_weighted = 0; # set back to zero to calculate weighted variance, save one variable declaration
## weighted stdev
# foreach my $size (keys %insert_sizes) {
    # $size_weighted = $size_weighted + (($size - $mean_weighted)*($size - $mean_weighted)*$insert_sizes{$size}); # variance
# }
# print "\tWeighted stdev: ", sqrt($size_weighted/$elements), "\n"; # square root for stdev


### Optionally, create results directory for output files
my $results_dir = './results';
if($opt_freq || $opt_histo) {
    mkdir $results_dir if (!-e $results_dir);
}


### Optionally, create insert size and read length distribution files
if ($opt_freq) {
    warn "The insert size and read length distribution files are created in the \'$results_dir\' directory:\n";
    my $ins_freq_out = 'ins_frequency.txt';
    open (my $ins_freq_fh, ">", "$results_dir"."/$ins_freq_out");
    warn "\t$ins_freq_out\n";
    print $ins_freq_fh "sizes\tcounts\n";
    foreach my $size (sort {$a <=> $b} keys %insert_sizes) {
        print $ins_freq_fh "$size\t$insert_sizes{$size}\n";
    }
    my $read_freq_out = 'read_frequency.txt';
    open (my $read_freq_fh, ">", "$results_dir"."/$read_freq_out");
    warn "\t$read_freq_out\n\n";
    print $read_freq_fh "lengths\tcounts\n";
    foreach my $length (sort {$a <=> $b} keys %read_lengths) {
        print $read_freq_fh "$length\t$read_lengths{$length}\n";
    }
    close $ins_freq_fh;
    close $read_freq_fh;
}


### Optionally, create histogram pdfs for the distribution of insert size and read length
if ($opt_histo) {
    warn "The following histogram pdfs are created in the \'$results_dir\' directory:\n";
    r_histo('ins', \@ins_size_array, \%ins_stats, $ref_name, $results_dir, $ins_xlim); # subroutine to create the insert size histogram with 'Rscript' of 'R'
    r_histo('read', \@read_length_array, \%read_stats, $ref_name, $results_dir, $read_xlim); # subroutine for read length histogram
    warn "\n";
}


close $sam_fh; # closing this filehandle to early will result in stdin a warning message
exit;


###############
# Subroutines #
###############

### Subroutine to print the precentage rounded to two decimal places
sub perc {
    my ($numerator, $denominator) = @_;
    printf("%.2f", ($numerator/$denominator)*100);
    return 1;
}


### Subroutine to plot a histogram by creating an R script and executing it with 'Rscript'
sub r_histo {
    my ($prefix, $data_array_ref, $hash_stat_ref, $ref_name, $results_dir, $xlim) = @_; # prefix is either 'ins' for insert size or 'read' for read length
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
    open (my $r_fh, ">", "$tmp_r_script");
    print $r_fh "#!/usr/bin/Rscript --vanilla --slave\n"; # header of R script
    print $r_fh "$prefix = scan(file=\"$tmp_file\", quiet=T)\n"; # scan in data and suppress sdtout output
    print $r_fh "pdf(\"$results_dir/$histo_name\")\n"; # filename of the pdf
    print $r_fh "hist($prefix, breaks=$xlim, xlim=c(min($prefix),$xlim), main=NULL, xlab=\"$label [bp]\", ylab=\"$label count\")\n"; # create the histogram with labels, but without MAIN-title (see below)
    if ($ref_name && $ref_count == 1) { # title for histogram plot
        print $r_fh "title(\"$label distribution for mapping on\\n$ref_name\", adj=0)\n"; # align title left to have room for mtext below
    } else {
        print $r_fh "title(\"$label distribution\")\n";
    }
    print $r_fh "abline(v=$hash_stat_ref->{'median'}, col=\"blue\", lwd=1)\n"; # plot the calculated median into the histogram
    print $r_fh "abline(v=$hash_stat_ref->{'mean'}, col=\"green\", lwd=1)\n";  # mean
    print $r_fh "legend(\"topright\", c(\"median\", \"mean\"), cex=0.8, col=c(\"blue\", \"green\"), lwd=1)\n";
    print $r_fh "mtext(\"median=$hash_stat_ref->{'median'}\\nmean=$hash_stat_ref->{'mean'}\\nstdev=$hash_stat_ref->{'sd'}\", side=3, adj=1, cex=0.8, col=\"red\")\n"; # add calculated stats to the plot margin
    print $r_fh "out <- dev.off()\n"; # write histogram to the pdf and suppress stdout output by diverting it
    close $r_fh;

    # execute R script with Rscript
    system("Rscript $tmp_r_script") == 0 or die "### Fatal error: Statistical programming language \'R\' either not installed, not in \$PATH, or something wrong with \'$tmp_r_script\'! Install \'R\' to create the histograms (required is \'Rscript\')!\n";
    warn "\t$histo_name\n"; # print to stderr which file has been created
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
    # $stat->clear(); # remove stats from the module for next round; not needed because of 'new()' initialisation at begin of subroutine
    return 1;
}
