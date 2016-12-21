revcom_seq
==========

`revcom_seq.pl` is a script to reverse complement (multi-)sequence files.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
* [Options](#options)
* [Output](#output)
* [Run environment](#run-environment)
* [Dependencies](#dependencies)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl revcom_seq.pl seq-file.embl > seq-file_revcom.embl

**or**

    perl cat_seq.pl multi-seq_file.embl | perl revcom_seq.pl -i embl > seq_file_cat_revcom.embl

## Description

This script reverse complements (multi-)sequence files. The
features/annotations in RichSeq files (e.g. EMBL or GENBANK format)
will also be adapted accordingly. Use option **-o** to specify a
different output sequence format. Input files can be given directly via
*STDIN* or as a file. If *STDIN* is used, the input sequence file
format has to be given with option **-i**. Be careful to set the
correct input format.

## Usage

    perl revcom_seq.pl -o gbk seq-file.embl > seq-file_revcom.gbk

**or** reverse complement all sequence files in the current working directory:

    for file in *.embl; do perl revcom_seq.pl -o fasta "$file" > "${file%.embl}"_revcom.fasta; done

## Options

- **-h**, **-help**

    Help (perldoc POD)

- **-o**=*str*, **-outformat**=*str*

    Specify different sequence format for the output [fasta, embl, or gbk]

- **-i**=*str*, **-informat**=*str*

    Specify the input sequence file format, only needed for *STDIN* input

- **-v**, **-version**

    Print version number to *STDOUT*

## Output

- *STDOUT*

    The reverse complemented sequence file is printed to *STDOUT*.
    Redirect or pipe into another tool as needed.

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies

- [**BioPerl**](http://www.bioperl.org) (tested version 1.007001)

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2 (2015-12-10)
    * included a POD instead of a simple usage text
    * included `pod2usage` with Pod::Usage
    * included 'use autodie' pragma
    * options with Getopt::Long
    * output format now specified with option **-o**
    * included version switch, **-v**
    * allowed file and *STDIN* input, instead of only file; thus new option **-i** for input format
    * output printed to *STDOUT* now, instead of output file
    * fixed bug, that only first sequence in multi-sequence file is reverse complemented. Now all sequences in a multi-seq file are reverse complemented.
* v0.1 (2013-02-08)
