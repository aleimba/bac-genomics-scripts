sample_fastx-txt
================

`sample_fastx-txt.pl` is a script to randomly subsample FASTA, FASTQ, or TEXT files.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [Subsample paired-end read data and retain pairing](#subsample-paired-end-read-data-and-retain-pairing)
  * [Subsample TEXT file and skip three header lines during subsampling](#subsample-text-file-and-skip-three-header-lines-during-subsampling)
  * [Subsample TEXT file and remove two header lines for final output](#subsample-text-file-and-remove-two-header-lines-for-final-output)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Acknowledgements](#acknowledgements)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl sample_fastx-txt.pl -i infile.fasta -n 100 > subsample.fasta

**or**

    zcat reads.fastq.gz | perl sample_fastx-txt.pl -i - -n 100000 > subsample.fastq

## Description

Randomly subsample FASTA, FASTQ, and TEXT files.

Empty lines in the input files will be skipped and not included in
sampling. Format TEXT assumes one entry per single line. FASTQ
format assumes **four** lines per read, if this is not the case run
the FASTQ file through [`fastx_fix.pl`](/fastx_fix) or use Heng
Li's [`seqtk seq`](https://github.com/lh3/seqtk):

    seqtk seq -l 0 infile.fq > outfile.fq

The file type is detected automatically. However, if automatic
detection fails, TEXT format is assumed. As a last resort, you can
set the file type manually with option **-f**.

This script is an implementation of the *reservoir sampling*
algorithm (or *Algorithm R (3.4.2)*) described in Donald Knuth's
[*The Art of Computer Programming*](https://en.wikipedia.org/wiki/The_Art_of_Computer_Programming).
It is designed to randomly pull a small sample size from a
(potential) huge input file of indeterminate size, which
(potentially) doesn't fit into main memory. The beauty of reservoir
sampling is that it requires only one pass through the input file.
The memory consumption of the algorithm is proportional to the
sample size, thus large sample sizes will consume lots of memory as
the whole sample will be held in memory. On the other hand, the size
of the initial file is irrelevant.

An alternative tool, which is a lot faster, is `seqtk sample` from
the [*seqtk toolkit*](https://github.com/lh3/seqtk>).

## Usage

### Subsample paired-end read data and retain pairing

    perl sample_fastx-txt.pl -i read-pair_1.fq -n 1000000 -s 123 > sub-pair_1.fq

    perl sample_fastx-txt.pl -i read-pair_2.fq -n 1000000 -s 123 > sub-pair_2.fq

### Subsample TEXT file and skip three header lines during subsampling

    perl sample_fastx-txt.pl -i infile.txt -n 100 -f text -t 3 > subsample.txt

### Subsample TEXT file and remove two header lines for final output

    perl sample_fastx-txt.pl -i infile.txt -n 350 -t 2 | sed '1,2d' > sub_no-header.txt

## Options

### Mandatory options

- -i, -input

    Input FASTA/Q or TEXT file, or piped *STDIN* (-)

- -n, -num

    Number of entries/reads to subsample

### Optional options

- -h, -help

    Help (perldoc POD)

- -f, -file_type

    Set the file type manually [fasta|fastq|text]

- -s, -seed

    Set starting random seed. For **paired-end** read data use the **same random seed** for both FASTQ files with option **-s** to retain pairing (see [Subsample paired-end read data and retain pairing](#subsample-paired-end-read-data-and-retain-pairing) above).

- -t, -title_skip

    Skip the specified number of header lines in TEXT files before subsampling and append them again afterwards. If you want to get rid of the header as well, pipe the subsample output to [`sed`](https://www.gnu.org/software/sed/manual/sed.html) (see `man sed` and [Subsample TEXT file and remove two header lines for final output](#subsample-text-file-and-remove-two-header-lines-for-final-output) above).

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    The subsample of the input file is printed to *STDOUT*. Redirect or pipe into another tool as needed.

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Acknowledgements

I got the idea for reservoir sampling from Sean Eddy's keynote at
the Janelia meeting on [*High Throughput Sequencing for Neuroscience*](http://cryptogenomicon.wordpress.com/2014/11/01/high-throughput-sequencing-for-neuroscience/)
which he posted in his blog
[*Cryptogenomicon*](http://cryptogenomicon.wordpress.com/). The [*Wikipedia article*](https://en.wikipedia.org/wiki/Reservoir_sampling) and the
[*PerlMonks*](http://www.perlmonks.org/index.pl?node_id=177092) implementation helped a lot, as well.

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.1 (18.11.2014)
