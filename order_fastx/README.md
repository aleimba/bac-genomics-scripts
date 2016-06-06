order_fastx
===========

`order_fastx.pl` is a script to order sequences in FASTA or FASTQ files.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)


## Synopsis

    perl order_fastx.pl -i infile.fasta -l order_id_list.txt > ordered.fasta

## Description

Order sequence entries in FASTA or FASTQ sequence files according to
an ID list with a given order. Beware, the IDs in the order list
have to be **identical** to the entire IDs in the sequence file.

However, the ">" or "@" ID identifiers of FASTA or FASTQ files,
respectively, can be omitted in the ID list.

The file type is detected automatically. But, you can set the file
type manually with option **-f**. FASTQ format assumes **four** lines
per read, if this is not the case run the FASTQ file through
[`fastx_fix.pl`](/fastx_fix) or use Heng Li's [`seqtk
seq`](https://github.com/lh3/seqtk):

    seqtk seq -l 0 infile.fq > outfile.fq

The script can also be used to pull a subset of sequences in the ID
list from the sequence file. Probably best to set option flag **-s**
in this case, see [Optional options](#optional-options) below. But, rather use
[`filter_fastx.pl`](/filter_fastx).

## Usage

    perl order_fastx.pl -i infile.fq -l order_id_list.txt -s -f fastq > ordered.fq

    perl order_fastx.pl -i infile.fasta -l order_id_list.txt -e > ordered.fasta

## Options

### Mandatory options

- -i, -input

    Input FASTA or FASTQ file

- -l, -list

    List with sequence IDs in specified order

### Optional options

- -h, -help

    Help (perldoc POD)

- -f, -file_type

    Set the file type manually [fasta|fastq]

- -e, -error_files

    Write missing IDs in the seq file or the order ID list without an equivalent in the other to error files instead of *STDERR* (see [Output](#output) below)

- -s, -skip_errors

    Skip missing ID error statements, excludes option **-e**

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    The newly ordered sequences are printed to *STDOUT*. Redirect or pipe into another tool as needed.

- (order_ids_missing.txt)

    If IDs in the order list are missing in the sequence file with option **-e**

- (seq_ids_missing.txt)

    If IDs in the sequence file are missing in the order ID list with option **-e**

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.1 (20.11.2014)
