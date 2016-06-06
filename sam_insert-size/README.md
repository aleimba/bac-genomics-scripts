sam_insert-size
===============

`sam_insert-size.pl` is a script to calculate insert size and read length statistics for paired-end reads in SAM/BAM format.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Dependencies](#dependencies)
* [Author - contact](#author---contact)
* [Acknowledgements](#acknowledgements)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl sam_insert-size.pl -i file.sam

**or**

    samtools view -h file.bam | perl sam_insert-size.pl -i -

## Description

Calculate insert size and read length statistics for paired-end reads
in SAM/BAM alignment format. The program gives the arithmetic mean,
median, and standard deviation (stdev) among other statistical values.

Insert size is defined as the total length of the original fragment
put into sequencing, i.e. the sequenced DNA fragment between the
adaptors. The 16-bit FLAG of the SAM/BAM file is used to filter reads
(see the [SAM specifications](http://samtools.sourceforge.net/SAM1.pdf)).

**Read length** statistics are calculated for all mapped reads
(irrespective of their pairing).

**Insert size** statistics are calculated only for **paired reads**.
Typically, the insert size is perturbed by artifacts, like chimeras,
structural re-arrangements or alignment errors, which result in a
very high maximum insert size measure. As a consequence the mean and
stdev can be strongly misleading regarding the real distribution. To
avoid this, two methods are implemented that first trim the insert
size distribution to a 'core' to calculate the respective statistics.
Additionally, secondary alignments for multiple mapping reads and
supplementary alignments for chimeric reads, as well as insert sizes
of zero are not considered (option **-min_ins_cutoff** is set to
**one** by default).

The **-a|-align** method includes only proper/concordant paired reads
in the statistical calculations (as determined by the mapper and the
options for insert size minimum and maximum used for mapping). This
is the **default** method.

The **-p|-percentile** method first calculates insert size statistics
for all read pairs, where the read and the mate are mapped ('raw
data'). Subsequently, the 10th and the 90th percentile are discarded
to calculate the 10% truncated mean and stdev. Discarding the lowest
and highest 10% of insert sizes gives the advantage of robustness
(insensitivity to outliers) and higher efficiency in heavy-tailed
distributions.

Alternative tools, which are a lot faster, are [`CollectInsertSizeMetrics`](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics)
from [Picard Tools](https://broadinstitute.github.io/picard/) and
[`sam-stats`](https://code.google.com/p/ea-utils/wiki/SamStats) from
[ea-utils](https://code.google.com/p/ea-utils/).

## Usage

    samtools view -h file.bam | perl sam_insert-size.pl -i - -p -d -f -min 50 -max 500 -n 2000000 -xlim_i 350 -xlim_r 200

## Options

### Mandatory options

- -i, -input

    Input SAM file or piped *STDIN* (-) from a BAM file e.g. with [`samtools view`](http://www.htslib.org/doc/samtools-1.1.html) from [Samtools](http://www.htslib.org/)

- -a, -align

    **Default method:** Align method to calculate insert size statistics, includes only reads which are mapped in a proper/concordant pair (as determined by the mapper). Excludes option **-p**.

**or**

- -p, -percentile

    Percentile method to calculate insert size statistics, includes only read pairs with an insert size within the 10th and the 90th percentile range of all mapped read pairs. However, the frequency distribution as well as the histogram will be plotted with the 'raw' insert size data before percentile filtering. Excludes option **-a**.

### Optional options

- -h, -help

    Help (perldoc POD)

- -d, -distro

    Create distribution histograms for the insert sizes and read lengths with [R](http://www.r-project.org/). The calculated median and mean (that are printed to *STDOUT*) are plotted as vertical lines into the histograms. Use it to control the correctness of the statistical calculations.

- -f, -frequencies

    Print the frequencies of the insert sizes and read lengths to tab-delimited files 'ins_frequency.txt' and 'read_frequency.txt', respectively.

- -max, -max_ins_cutoff

    Set a maximal insert size cutoff, all insert sizes above this cutoff will be discarded (doesn't affect read length). With **-min** and **-max** you can basically run both methods, by first running the script with **-p** and then using the 10th and 90th percentile of the 'raw data' as **-min** and **-max** for option **-a**.

- -min, -min_ins_cutoff

    Set a minimal insert size cutoff [default = 1]

- -n, -num_read

    Number of reads to sample for the calculations from the start of the SAM/BAM file. Significant statistics can usually be calculated from a fraction of the total SAM/BAM alignment file.

- -xlim_i, -xlim_ins

    Set an upper limit for the x-axis of the **'R' insert size** histogram, overriding automatic truncation of the histogram tail. The default cutoff is one and a half times the third quartile Q3 (75th percentile) value. The minimal cutoff is set to the lowest insert size automatically. Forces option **-d**.

- -xlim_r, -xlim_read

    Set an upper limit for the x-axis of the optional **'R' read length** histogram. Default value is as in **-xlim_i**. Forces option **-d**.

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    Calculated stats are printed to *STDOUT*

- ./results

    All **optional** output files are stored in this results folder

- (./results/ins_frequency.txt)

    Frequencies of insert size 'raw data', tab-delimited

- (./results/ins_histo.pdf)

    Distribution histogram for the insert size 'raw data'

- (./results/read_frequency.txt)

    Frequencies of read lengths, tab-delimited

- (./results/read_histo.pdf)

    Distribution histogram for the read lengths. Not informative if there's no variation in the read lengths.

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies

- `Statistics::Descriptive`

    Perl module to calculate descriptive statistics, if not installed already get it from [CPAN](http://www.cpan.org/)

- Statistical computing language [R](http://www.r-project.org/)

    `Rscript` is needed to plot the histograms with option **-d**

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Acknowledgements

References/thanks go to:

- Tobias Rausch's online courses/workshops (EMBL Heidelberg) on the introduction to SAM files and flags (http://www.embl.de/~rausch/)

- The CBS NGS Analysis course for the percentile filtering idea (http://www.cbs.dtu.dk/courses/27626/programme.php)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.2 (29.10.2014)
    - Fixed bug for options '-min_ins_size' and '-max_ins_size'
    - warn if result files already exist
    - simplify prints to R script with Perl function 'select'
    - minor Perl syntax changes so all Perl scripts conform to the same syntax
    - minor changes to POD
    - finally included README.md
- v0.1 (27.11.2013)
