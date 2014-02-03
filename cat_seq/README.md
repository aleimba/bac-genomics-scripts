cat_seq
=======

A script to merge multi-sequence RichSeq files into one single-entry 'artificial' sequence file.

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
If you're working only with fasta files UNIXs *grep* is a faster choice to concatenate sequences.

    grep -v ">" seq.fasta > seq_artificial.fasta

Subsequently add as a first line a fasta ID (starting with '>') with an editor.

## Output

* *\_artificial.[embl|fasta|genbank]

Concatenated artificial sequence in the input format, or optionally the specified output sequence format.

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies (not in the core Perl modules)

* BioPerl (tested with version 1.006901)

## Alternative software

The EMBOSS (The European Molecular Biology Open Software Suite) application ***union*** can also be used for this task (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/union.html).

## Author/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

* v0.1 (08.02.2013)