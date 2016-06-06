rename_fasta_id
===============

`rename_fasta_id.pl` is a script to rename fasta IDs according to regular expressions.

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

    perl rename_fasta_id.pl -i file.fasta -p "NODE_.+$" -r "K-12_" -n -a c > out.fasta

**or**

    zcat file.fasta.gz | perl rename_fasta_id.pl -i - -p "coli" -r "" -o > out.fasta

## Description

This script uses the built-in Perl substitution operator `s///` to
replace strings in FASTA IDs. To do this, a **pattern** and a
**replacement** have to be provided (Perl regular expression syntax
can be used). The leading '>' character for the FASTA ID will be
removed before the substitution and added again afterwards. FASTA
IDs will be searched for matches with the **pattern**, and if found
the **pattern** will be replaced by the **replacement**.

**IMPORTANT**: Enclose the **pattern** and the **replacement** in
quotation marks (' or ") if they contain characters that would be
interpreted by the shell (e.g. pipes '|', brackets etc.).

For substitutions without any appendices in a UNIX OS you can of
course just use the great
[`sed`](https://www.gnu.org/software/sed/manual/sed.html) (see
`man sed`), e.g.:

    sed 's/^>pattern/>replacement/' file.fasta

## Usage

    perl rename_fasta_id.pl -i file.fasta -p "T" -r "a" -c -g -o

## Options

### Mandatory options

- -i, -input

Input FASTA file or piped STDIN (-) from a gzipped file

- -p, -pattern

Pattern to be replaced in FASTA ID

- -r, -replacement

Replacement to replace the pattern with. To entirely remove the pattern use '' or "" as input for **-r**.

### Optional options

- -h, -help

Help (perldoc POD)

- -c, -case-insensitive

Match pattern case-insensitive

- -g, -global

Replace pattern globally in the string

- -n, -numerate

Append a numeration/the count of the pattern hits to the replacement. This is e.g. useful to number contigs consecutively in a draft genome.

- -a, -append

Append a string after the numeration, e.g. 'c' for chromosome

- -o, -output

Verbose output of the substitutions that were carried out, printed to *STDERR*

- -v, -version

Print version number to *STDERR*

## Output

- *STDOUT*

The FASTA file with substituted ID lines is printed to *STDOUT*. Redirect or pipe into another tool as needed.

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

- v0.1 (09.11.2014)
