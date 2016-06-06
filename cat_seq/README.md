cat_seq
=======

A script to merge multi-sequence RichSeq files into one single-entry 'artificial' sequence file.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [Merge multi-sequence file](#merge-multi-sequence-file)
  * [Merge multi-sequence file and specify different output format](#merge-multi-sequence-file-and-specify-different-output-format)
  * [UNIX loop to concatenate each multi-sequence file in the current working directory](#unix-loop-to-concatenate-each-multi-sequence-file-in-the-current-working-directory)
  * [Concatenate multi-sequence fasta files faster with UNIX's `grep`](#concatenate-multi-sequence-fasta-files-faster-with-unixs-grep)
* [Output](#output)
* [Dependencies](#dependencies)
* [Run environment](#run-environment)
* [Alternative software](#alternative-software)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl cat_seq.pl multi-seq_file.embl

## Description

This script concatenates multiple sequences in a RichSeq file (embl or genbank, but also fasta) to a single artificial sequence. The first sequence in the file is used as a foundation to add the subsequent sequences, along with all features and annotations.

Optionally, a different output file format can be specified (fasta/embl/genbank).

## Usage

### Merge multi-sequence file

    perl cat_seq.pl multi-seq_file.gbk

### Merge multi-sequence file and specify different output format

    perl cat_seq.pl multi-seq_file.embl [fasta|genbank]

### UNIX loop to concatenate each multi-sequence file in the current working directory

    for i in *.[embl|fasta|gbk]; do perl cat_seq.pl $i [embl|fasta|genbank]; done

### Concatenate multi-sequence fasta files faster with UNIXs *grep*
If you're working only with fasta files UNIX's `grep` is a faster choice to concatenate sequences.

    grep -v ">" seq.fasta > seq_artificial.fasta

Subsequently add as a first line a fasta ID (starting with '>') with an editor.

## Output

* *\_artificial.[embl|fasta|genbank]

Concatenated artificial sequence in the input format, or optionally the specified output sequence format.

## Dependencies

* BioPerl (tested with version 1.006901)

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Alternative software

The EMBOSS (The European Molecular Biology Open Software Suite) application ***union*** can also be used for this task (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/union.html).

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.1 (08.02.2013)
