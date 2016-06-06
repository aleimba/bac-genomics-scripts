genomes_feature_table
=====================

`genomes_feature_table.pl` is a script to create a feature table for genomes in EMBL and GENBANK format.

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

    perl genomes_feature_table.pl path/to/genome_dir > feature_table.tsv

## Description

A genome feature table lists basic stats/info (e.g. genome size, GC
content, coding percentage, accession number(s)) and the numbers of
annotated primary features (e.g. CDS, genes, RNAs) of genomes. It
can be used to have an overview of these features in different
genomes, e.g. in comparative genomics publications.

`genomes_feature_table.pl` is designed to extract (or calculate)
these basic stats and **all** annotated primary features from RichSeq
files (**EMBL** or **GENBANK** format) in a specified directory (with the
correct file extension, see option **-e**). The **default** directory
is the current working directory. The primary features are
counted and the results for each genome printed in tab-separated
format. It is a requirement that each file contains **only one**
genome (complete or draft, with or without plasmids).

The most important features will be listed first, like genome
description, genome size, GC content, coding percentage (calculated
based on non-pseudo CDS annotation), CDS and gene numbers, accession
number(s) (first..last in the sequence file), RNAs (rRNA, tRNA,
tmRNA, ncRNA), and unresolved bases (IUPAC code 'N'). If plasmids are
annotated in a sequence file, the number of plasmids are
counted and listed as well (needs a */plasmid="plasmid_name"* tag in the
*source* primary tag, see e.g. Genbank accession number
[CP009167](http://www.ncbi.nlm.nih.gov/nuccore/CP009167)). Use option **-p**
to list plasmids as separate entries (lines) in the feature table.

For draft genomes the number of contigs/scaffolds are counted. All
contigs/scaffolds of draft genomes should be marked with the *WGS*
keyword (see e.g. draft NCBI Genbank entry
[JSAY00000000](http://www.ncbi.nlm.nih.gov/nuccore/JSAY00000000)). If this is
not the case for your file(s) you can add those keywords to each
sequence entry with the following Perl one-liners (will
edit files in place). For files in **GENBANK** format if 'KEYWORDS&nbsp;&nbsp;&nbsp;&nbsp;.' is present

    perl -i -pe 's/^KEYWORDS(\s+)\./KEYWORDS$1WGS\./' file

or if 'KEYWORDS' isn't present at all

    perl -i -ne 'if(/^ACCESSION/){ print; print "KEYWORDS    WGS.\n";} else{ print;}' file

For files in **EMBL** format if 'KW&nbsp;&nbsp;&nbsp;.' is present

    perl -i -pe 's/^KW(\s+)\./KW$1WGS\./' file

or if 'KW' isn't present at all

    perl -i -ne 'if(/^DE/){ $dw=1; print;} elsif(/^XX/ && $dw){ print; $dw=0; print "KW   WGS.\n";} else{ print;}' file

## Usage

    perl genomes_feature_table.pl -p -e gb,gbk > feature_table_plasmids.tsv

    perl genomes_feature_table.pl path/to/genome_dir/ -e gbf -e embl > feature_table.tsv

## Options

- -h, -help

    Help (perldoc POD)

- -e, -extensions

    File extensions to include in the analysis (EMBL or GENBANK format),
    either comma-separated list or multiple occurences of the option
    [default = ebl,emb,embl,gb,gbf,gbff,gbank,gbk,genbank]

- -p, -plasmids

    Optionally list plasmids as extra entries in the feature table, if
    they are annotated with a */plasmid="plasmid_name"* tag in the
    *source* primary tag

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    The resulting feature table is printed to *STDOUT*. Redirect or
    pipe into another tool as needed (e.g. `cut`, `grep`, or `head`).

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies

- [BioPerl](http://www.bioperl.org) (tested version 1.006923)

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.5 (14.09.2015)
    - changed script name to `genomes_feature_table.pl`
    - included a POD
    - options with Getopt::Long
    - included `pod2usage` with Pod::Usage
    - major code overhaul with restructuring (removing code redundancy, print out without temp file etc.) and Perl syntax changes
    - changed input options to get folder path from STDIN
    - as a consequence new option **-e|-extensions**
    - accession numbers not essential anymore, changed hash key to filename; but requires now only one genome per file
    - draft genomes should include 'WGS' keyword (warning if not)
    - option **-p|-plasmids** works now correctly with complete and draft genomes
    - count plasmids without option **-p**
- v0.4 (11.08.2013)
    - included 'use autodie;' pragma
    - included version switch
- v0.3 (05.11.2012)
    - new option **p** to report plasmid features in multi-sequence draft files separately
- v0.2 (19.09.2012)
- v0.1 (25.11.2011)
    - **original** script name: `get_genome_features.pl`
