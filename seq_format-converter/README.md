seq_format-converter
====================

A script to convert a sequence file to another format.

## Synopsis

    perl seq_format-converter.pl -i seq_file.gbk -f gbk -o embl

## Description

This script converts a (multi-)sequence file of a specific format to a differently formatted output file. The most common sequence formats are: **embl**, **fasta**, and **gbk** (genbank).

Since sequence formats change from time to time, BioPerl is not always up to date. For all available BioPerl sequence formats see: http://www.bioperl.org/wiki/HOWTO:SeqIO#Formats. **Warning**: The *bioperl-ext* package and the *io_lib* library from the **Staden** package (http://staden.sourceforge.net/) need to be installed in order to read the scf, abi, alf, pln, exp, ctf, ztr formats.

## Usage

    perl seq_format-converter.pl -i seq_file -f in_format -o out_format

### UNIX loop to reformat all sequence files in the current working directory

    for i in *.[embl|gbk]; do perl seq_format-converter.pl -i $i -f [embl|gbk] -o [embl|fasta|gbk]; done

## Options for *seq_format-converter.pl*

### Mandatory options

* -i, -input

Input sequence file

* -f, -format

Input sequence format (e.g. 'embl' or 'gbk)

* -o, -out_format

Output sequence format (e.g. 'embl', 'fasta' or 'gbk)

### Optional options

* -h, -help

Print usage

* -v, -version

Print version number

## Output

* seq_file.[embl|fasta|gbk]

Output sequence file in the specified format

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies (not in the core Perl modules)

* BioPerl (tested with version 1.006901)

## Author/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2 (03.02.2014)
    - allow short 'gbk' format instead of 'genbank'
    - also short 'gbk' file-extension for output file
    - included 'use autodie'
    - usage as HERE document
    - options with Getopt::Long
    - version switch
* v0.1 (10.11.2011)
