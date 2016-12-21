trunc_seq
=========

`trunc_seq.pl` is a script to truncate sequence files.

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

    perl trunc_seq.pl 20 3500 seq-file.embl > seq-file_trunc_20_3500.embl

**or**

    perl trunc_seq.pl file_of_filenames_and_coords.tsv

## Description

This script truncates sequence files according to the given
coordinates. The features/annotations in RichSeq files (e.g. EMBL or
GENBANK format) will also be adapted accordingly. Use option **-o** to
specify a different output sequence format. Input can be given directly
as a file and truncation coordinates to the script, with the start
position as the first argument, stop as the second and (the path to)
the sequence file as the third. In this case the truncated sequence
entry is printed to *STDOUT*. Input sequence files should contain only
one sequence entry, if a multi-sequence file is used as input only the
**first** sequence entry is truncated.

Alternatively, a file of filenames (fof) with respective coordinates
and sequence files in the following **tab-separated** format can be
given to the script (the header is optional):

\#start&emsp;stop&emsp;seq-file<br>
300&emsp;9000&emsp;(path/to/)seq-file<br>
50&emsp;1300&emsp;(path/to/)seq-file2<br>

With a fof the resulting truncated sequence files are printed into a
results directory. Use option **-r** to specify a different results
directory than the default.

It is also possible to truncate a RichSeq sequence file loaded into the
[Artemis](http://www.sanger.ac.uk/science/tools/artemis) genome browser
from the Sanger Institute: Select a subsequence and then go to Edit ->
Subsequence (and Features)

## Usage

    perl trunc_seq.pl -o gbk 120 30000 seq-file.embl > seq-file_trunc_120_3000.gbk

**or**

    perl trunc_seq.pl -o fasta 5300 18500 seq-file.gbk | perl revcom_seq.pl -i fasta > seq-file_trunc_revcom.fasta

**or**

    perl trunc_seq.pl -r path/to/trunc_embl_dir -o embl file_of_filenames_and_coords.tsv

## Options

- **-h**, **-help**

    Help (perldoc POD)

- **-o**=*str*, **-outformat**=*str*

    Specify different sequence format for the output (files) [fasta, embl, or gbk]

- **-r**=*str*, **-result\_dir**=*str*

    Path to result folder for fof input \[default = './trunc\_seq\_results'\]

- **-v**, **-version**

    Print version number to *STDOUT*

## Output

- *STDOUT*

    If a single sequence file is given to the script the truncated sequence
    file is printed to *STDOUT*. Redirect or pipe into another tool as
    needed.

**or**

- ./trunc_seq_results

    If a fof is given to the script, all output files are stored in a
    results folder

- ./trunc_seq_results/seq-file_trunc_start_stop.format

    Truncated output sequence files are named appended with 'trunc' and the
    corresponding start and stop positions

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies

- [**BioPerl**](http://www.bioperl.org) (tested version 1.007001)

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2 (2015-12-07)
    * Merged funtionality of `trunc_seq.pl` and `run_trunc_seq.pl` in one single script
        * Allows now single file and file of filenames (fof) with coordinates input
        * output for single file input printed to *STDOUT* now
        * output for fof input printed into files in a result directory, new option **-r** to specify result directory
    * included a POD instead of a simple usage text
    * included `pod2usage` with Pod::Usage
    * included 'use autodie' pragma
    * options with Getopt::Long
    * output format now specified with option **-o**
    * included version switch, **-v**
    * fixed bug to remove input filepaths from fof input for output files
    * skip empty or comment lines (/^#/) in fof input
    * check and warn if input seq file has more than one seq entries
* v0.1 (2013-02-08)
    * In v0.1 `trunc_seq.pl` only for single sequence input, but included additional wrapper script `run_trunc_seq.pl` for a fof input
