calc_fastq-stats
================

`calc_fastq-stats.pl` is a script to calculate basic statistics for bases and reads in a FASTQ file.

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
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl calc_fastq-stats.pl -i reads.fastq

**or**

    gzip -dc reads.fastq.gz | perl calc_fastq-stats.pl -i -

## Description

The script calculates some simple statistics, like individual and total base
counts, GC content, and basic stats for the read lengths, and
read/base qualities in a FASTQ file. The GC content calculation does
not include 'N's. Stats are printed to *STDOUT* and optionally to an
output file.

Because the quality of a read degrades over its length with all NGS
machines, it is advisable to also plot the quality for each cycle as
implemented in tools like
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
or the [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/).

If the sequence and the quality values are interrupted by line
breaks (i.e. a read is **not** represented by four lines), please fix
with Heng Li's [seqtk](https://github.com/lh3/seqtk):

    seqtk seq -l 0 infile.fastq > outfile.fastq

An alternative tool, which is a lot faster, is **fastq-stats** from
[ea-utils](https://code.google.com/p/ea-utils/).

## Usage

    zcat reads.fastq.gz | perl calc_fastq-stats.pl -i - -q 64 -c 175000000 -n 3000000

## Options

### Mandatory options

- -i, -input

Input FASTQ file or piped STDIN (-) from a gzipped file

- -q, -qual_offset

ASCII quality offset of the Phred (Sanger) quality values [default 33]

### Optional options

- -h, -help:

Help (perldoc POD)

- -c, -coverage_limit

Number of bases to sample from the top of the file

- -n, -num_read

Number of reads to sample from the top of the file

- -o, -output

Print stats in addition to *STDOUT* to the specified output file

- -v, -version

Print version number to *STDERR*

## Output

- *STDOUT*

Calculated stats are printed to *STDOUT*

- (outfile)

Optional outfile for stats

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies

If the following modules are not installed get them from
[CPAN](http://www.cpan.org/):

- `Statistics::Descriptive`

Perl module to calculate basic descriptive statistics

- `Statistics::Descriptive::Discrete`

Perl module to calculate descriptive statistics for discrete data sets

- `Statistics::Descriptive::Weighted`

Perl module to calculate descriptive statistics for weighted variates

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.1 (28.10.2014)
