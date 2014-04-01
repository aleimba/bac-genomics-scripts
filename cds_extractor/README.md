cds_extractor
=============

A script to extract amino acid or nucleotide sequences from coding sequence (CDS) features in annotated genomes.

## Synopsis

    perl cds_extractor.pl -i seq_file.[embl|gbk] -p

## Description

This script extracts protein or DNA sequences of CDS features from a (multi)-RichSeq file (e.g. embl, genbank) and writes them to a multi-fasta file. The fasta headers for each CDS include either the locus tag, if that's not available, protein id, gene, or an internal CDS counter as identifier (in this order). The organism info includes also possible plasmid names. Pseudogenes (tagged by '**/pseudo**') are not included (except in the CDS counter).

In addition to the identifier, fasta headers include gene (**g=**), product (**p=**), organism (**o=**), and EC numbers (**ec=**), if these are present for a CDS. Individual EC numbers are separated by **semicolons**. The location/position (**l=** start..stop) of a CDS will always be included. If gene is used as fasta header identifier '**g=** gene' will only be included with option **-f**.

Fuzzy locations in the feature table of a sequence file are not taken into consideration for '**l=**'. If you set options **-u** and/or **-d** and the feature location overlaps a **circular** replicon boundary, positions are marked with '<' or '>' in the direction of the exceeded boundary. Features with overlapping locations in **linear** sequences (e.g. contigs) will be skipped and are **not** included in the output! A CDS feature is on the lagging strand if start > stop in the location. In the special case of overlapping circular sequence boundaries this is reversed.

Of course, the '**l=**' positions are separate for each sequence in a multi- sequence file. Thus, if you want continuous positions for the CDSs run these files first through `cat_seq.pl`.

Optionally, a file with locus tags can be given to extract only these CDS features with option **-l** (each locus tag in a new line).

## Usage

### Extract amino acid sequences

    perl cds_extractor.pl -i Ecoli_MG1655.gbk -p [-l locus_tags.txt -c MG1655 -f]

### Extract nucleotide sequences

    perl cds_extractor.pl -i Banthracis_Ames.embl -n [-l locus_tags.txt -u 100 -d 20 -c Ames -f]

### UNIX loop to extract sequences from all files in the current working directory

    for i in *.embl; do perl cds_extractor.pl -i $i -p [-l locus_tags.txt]; done

## Options

### Mandatory options

* -i, -input

Input RichSeq sequence file including CDS annotation (e.g. embl or genbank)

* -p, -protein

Extract **protein** sequence for each CDS feature, excludes option **-n**

* -n, -nucleotide

Extract **nucleotide** sequence for each CDS feature, excludes option **-p**

### Optional options

* -h, -help:   Help (perldoc POD)

* -u, -upstream

Include given number of flanking nucleotides upstream of each CDS feature, forces option **-n**

* -d, -downstream

Include given number of flanking nucleotides downstream of each CDS feature, forces option **-n**

* -c, -cds_prefix

Prefix for the internal CDS counter [default = 'CDS']

* -l, -locustag_list

List of locus tags to extract only those (each locus tag on a new line)

* -f, -full_header

If gene is used as identifier include additionally '**g=** gene' in fasta headers, so downstream analyses can recognize the gene tag (e.g. `prot_finder.pl`).

* -v, -version

Print version number to STDERR

## Output

* *\_cds\_aa.fasta

Multi-fasta file of CDS protein sequences

**or**

* *\_cds\_nuc.fasta

Multi-fasta file of CDS DNA sequences

* (no_annotation_err.txt)

Lists input files missing CDS annotation, script exited with **fatal error** i.e. no fasta output file

* (double_id_err.txt)

Lists input files with ambiguous fasta identifiers, script exited with **fatal error** i.e. no fasta output file

* (locus_tag_missing_err.txt)

Lists CDS features without locus tags

* (linear_seq_cds_overlap_err.txt)

Lists CDS features overlapping sequence border of a **linear** molecule, which are **not** included in the result multi-fasta file

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies (not in the core Perl modules)

* BioPerl (tested with version 1.006901)

## Author/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

* v0.7 (31.03.2014)
    - location (l=) and EC numbers (ec=) for CDS features are included in the fasta header
    - 'ec=', 'g=', 'p=', and 'o=' only included in fasta header if these tags are present for a CDS feature, or additionally for 'g=' with option **-f**
    - if, with options '-u' and/or '-d', the location of a CDS feature overlaps a sequence boundary, the positions are marked with '<' or '>' in 'l='
    - additionally, CDS features whose location overlaps the sequence boundary of a linear molecule will not be included in the output, but identifier written to an error file
    - new option **-c** to chose prefix for internal CDS counter
    - /product feature value will not be used as fasta identifier anymore, skip directly to internal CDS counter, if /locus_tag, /protein_id, or /gene is missing for a CDS (too many 'hypothetical proteins')
    - internal CDS counter counts all CDSs of multi-sequence files sequential (doesn't start new with each new sequence in the multi-sequence file)
    - 'control_double' subroutine also called if /gene is used as fasta identifier
    - fixed bug introduced in v0.6 to exit if no CDS primary features found, because a draft multi-sequence file might have unannotated small contigs
    - new error files: no_annotation_err.txt, double_id_err.txt, linear_seq_cds_overlap_err.txt (the first two come in handy if you run `cds_extractor.pl` in a UNIX loop with many files)
    - included 'use autodie'
    - included version switch
    - included pod2usage with Pod::Usage
    - reorganized code into more subroutines to remove useless double codings (which contained also some bugs) and to make the script more concise
    - minor changes to Perl syntax
* v0.6 (06.06.2013)
    - exit with error if no CDS primary features present in input file, as /translation feature only present in CDS features (some genbank files are only annotated with 'gene')
    - included Bio::SeqFeatureI's method *spliced-seq* for CDS with split nucleotide sequences (CDS position indicated by 'join')
    - minor changes how the optional list of locus tags is handled
* v0.5 (03.06.2013)
    - included a POD
    - options with Getopt::Long
    - option **-n** to alternatively extract nucleotide sequences for CDS features (optionally with upstream and downstream sequences)
    - option to include full fasta ID header for downstream `blast_prot_finder.pl` analysis
    - exit with error if the values for two (or more) /locus_tag or /protein_id tags are not unambiguous
    - print message to STDOUT if and which locus tags were not found in a given locus tag list (option **-l**)
* v0.4 (06.02.2013)
    - replace whitespaces of /product values with underscores
* v0.3 (06.09.2012)
    - internal CDS counter to use in fasta ID for CDS features without a /locus_tag, /protein_id, /gene, or /product tag
    - include also organism (and possible plasmid) information in fasta ID lines
    - give a warning to STDOUT if a CDS feature without a /locus_tag is found (but only for the first occurence)
    - additionally, *locus_tag_errors.txt* error file to list all CDSs without locus tags
    - catch errors with *eval* if a tag is missing
* v0.2 (04.09.2012)
    - if a CDS feature does not have a /locus_tag, then use the value for /protein_id, /gene, or /product (in this order) in the fasta ID lines of the result file
    - optional extract only CDSs with locus tags given in a file
* v0.1 (24.05.2012)
