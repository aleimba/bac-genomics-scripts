cds_extractor
=============

`cds_extractor.pl` is a script to extract amino acid or nucleotide sequences from coding sequence (CDS) features in annotated genomes.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [Extract amino acid sequences](#extract-amino-acid-sequences)
  * [Extract nucleotide sequences](#extract-nucleotide-sequences)
  * [UNIX loop to extract sequences from all files in the current working directory](#unix-loop-to-extract-sequences-from-all-files-in-the-current-working-directory)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Dependencies](#dependencies)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl cds_extractor.pl -i seq_file.[embl|gbk] -p

## Description

This script extracts protein or DNA sequences of CDS features from a (multi)-RichSeq file (e.g. EMBL or GENBANK format) and writes them to a multi-FASTA file. The FASTA headers for each CDS include either the locus tag, if that's not available, protein ID, gene, or an internal CDS counter as identifier (in this order). The organism info includes also possible plasmid names. Pseudogenes (tagged by **/pseudo**) are not included (except in the CDS counter).

In addition to the identifier, FASTA headers include gene (**g=**), product (**p=**), organism (**o=**), and EC numbers (**ec=**), if these are present for a CDS. Individual EC numbers are separated by **semicolons**. The location/position (**l=** start..stop) of a CDS will always be included. If gene is used as FASTA header ID '**g=** gene' will only be included with option **-f**.

Fuzzy locations in the feature table of a sequence file are not taken into consideration for **l=**. If you set options **-u** and/or **-d** and the feature location overlaps a **circular** replicon boundary, positions are marked with '<' or '>' in the direction of the exceeded boundary. Features with overlapping locations in **linear** sequences (e.g. contigs) will be skipped and are **not** included in the output! A CDS feature is on the lagging strand if start > stop in the location. In the special case of overlapping circular sequence boundaries this is reversed.

Of course, the **l=** positions are separate for each sequence in a multi- sequence file. Thus, if you want continuous positions for the CDSs run these files first through [`cat_seq.pl`](/cat_seq).

Optionally, a file with locus tags can be given to extract only these CDS features with option **-l** (each locus tag in a new line).

## Usage

### Extract amino acid sequences

    perl cds_extractor.pl -i Ecoli_MG1655.gbk -p [-l locus_tags.txt -c MG1655 -f]

### Extract nucleotide sequences

    perl cds_extractor.pl -i Banthracis_Ames.embl -n [-l locus_tags.txt -u 100 -d 20 -c Ames -f]

### UNIX loop to extract sequences from all files in the current working directory

    for file in *.embl; do perl cds_extractor.pl -i "$file" -p [-l locus_tags.txt]; done

## Options

### Mandatory options

* **-i**=_str_, **-input**=_str_

    Input RichSeq sequence file including CDS annotation (e.g. EMBL or GENBANK format)

* **-p**, **-protein**

    Extract **protein** sequence for each CDS feature, excludes option **-n**

**or**

* **-n**, **-nucleotide**

    Extract **nucleotide** sequence for each CDS feature, excludes option **-p**

### Optional options

* **-h**, **-help**

    Help (perldoc POD)

* **-u**=_int_, **-upstream**=_int_

    Include given number of flanking nucleotides upstream of each CDS feature, forces option **-n**

* **-d**=_int_, **-downstream**=_int_

    Include given number of flanking nucleotides downstream of each CDS feature, forces option **-n**

* **-c**=_str_, **-cds_prefix**=_str_

    Prefix for the internal CDS counter [default = 'CDS']

* **-l**=_str_, **-locustag_list**=_str_

    List of locus tags to extract only those (each locus tag on a new line)

* **-f**, **-full_header**

    If gene is used as ID include additionally '**g=** gene' in FASTA headers, so downstream analyses can recognize the gene tag (e.g. [`prot_finder.pl`](/prot_finder)).

* **-v**, **-version**

    Print version number to *STDERR*

## Output

* \*.faa

    Multi-FASTA file of CDS protein sequences

**or**

* \*.ffn

    Multi-FASTA file of CDS DNA sequences

* (no_annotation_err.txt)

    Lists input files missing CDS annotation, script exited with **fatal error** i.e. no FASTA output file

* (double_id_err.txt)

    Lists input files with ambiguous FASTA IDs, script exited with **fatal error** i.e. no FASTA output file

* (locus_tag_missing_err.txt)

    Lists CDS features without locus tags

* (linear_seq_cds_overlap_err.txt)

    Lists CDS features overlapping sequence border of a **linear** molecule, which are **not** included in the result multi-FASTA file

## Dependencies

* [BioPerl](http://www.bioperl.org) (tested with version 1.006923)

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.7.1 (26.10.2015)
    - changed output file extensions from **\_cds\_aa.fasta* or **\_cds\_nuc.fasta* to **.faa* or **.ffn*, respectively
    - minor syntax changes in README, included TOC
    - minor syntax changes in POD
* v0.7 (31.03.2014)
    - location (l=) and EC numbers (ec=) for CDS features are included in the FASTA header
    - 'ec=', 'g=', 'p=', and 'o=' only included in FASTA header if these tags are present for a CDS feature, or additionally for 'g=' with option **-f**
    - if, with options '-u' and/or '-d', the location of a CDS feature overlaps a sequence boundary, the positions are marked with '<' or '>' in 'l='
    - additionally, CDS features whose location overlaps the sequence boundary of a linear molecule will not be included in the output, but IDs written to an error file
    - new option **-c** to chose prefix for internal CDS counter
    - /product feature value will not be used as FASTA ID anymore, skip directly to internal CDS counter, if /locus_tag, /protein_id, or /gene is missing for a CDS (too many 'hypothetical proteins')
    - internal CDS counter counts all CDSs of multi-sequence files sequential (doesn't start new with each new sequence in the multi-sequence file)
    - 'control_double' subroutine also called if /gene is used as FASTA ID
    - fixed bug introduced in v0.6 to exit if no CDS primary features found, because a draft multi-sequence file might have unannotated small contigs
    - new error files: no_annotation_err.txt, double_id_err.txt, linear_seq_cds_overlap_err.txt (the first two come in handy if you run `cds_extractor.pl` in a UNIX loop with many files)
    - included 'use autodie'
    - included version switch
    - included pod2usage with Pod::Usage
    - reorganized code into more subroutines to remove useless double codings (which contained also some bugs) and to make the script more concise
    - minor changes to Perl syntax
* v0.6 (06.06.2013)
    - exit with error if no CDS primary features present in input file, as /translation feature only present in CDS features (some GENBANK files are only annotated with 'gene')
    - included Bio::SeqFeatureI's method *spliced-seq* for CDS with split nucleotide sequences (CDS position indicated by 'join')
    - minor changes how the optional list of locus tags is handled
* v0.5 (03.06.2013)
    - included a POD
    - options with Getopt::Long
    - option **-n** to alternatively extract nucleotide sequences for CDS features (optionally with upstream and downstream sequences)
    - option to include full FASTA ID header for downstream [`prot_finder.pl`](/prot_finder) analysis
    - exit with error if the values for two (or more) /locus_tag or /protein_id tags are not unambiguous
    - print message to *STDOUT* if and which locus tags were not found in a given locus tag list (option **-l**)
* v0.4 (06.02.2013)
    - replace whitespaces of /product values with underscores
* v0.3 (06.09.2012)
    - internal CDS counter to use in FASTA ID for CDS features without a /locus_tag, /protein_id, /gene, or /product tag
    - include also organism (and possible plasmid) information in FASTA ID lines
    - give a warning to *STDOUT* if a CDS feature without a /locus_tag is found (but only for the first occurence)
    - additionally, *locus_tag_errors.txt* error file to list all CDSs without locus tags
    - catch errors with *eval* if a tag is missing
* v0.2 (04.09.2012)
    - if a CDS feature does not have a /locus_tag, then use the value for /protein_id, /gene, or /product (in this order) in the FASTA ID lines of the result file
    - optional extract only CDSs with locus tags given in a file
* v0.1 (24.05.2012)
