cds_extractor
=============

A script to extract amino acid or nucleotide sequences from coding sequence (CDS) features in annotated genomes.

## Synopsis

    perl cds_extractor.pl -i seq_file.embl -p

## Description

This scripts extracts protein or DNA sequences of CDS features from a (multi)-embl or -genbankfile and writes them to a multi-fasta file. The fasta ID line includes either the locus tag (plus g=gene, p=product, o=organism; if existent), if that's not available protein id (plus g=gene, p=product, o=organism), gene (plus p=product, o=organism), product (plus o=organism), or an internal CDS counter (in this order). The organism info includes also possible plasmid names. Amino acid sequences are read from '/translation' tags. Pseudogenes (tagged by '/pseudo') are not included (except in the CDS counter).

CDSs without locus tags are written to an error file 'locus_tag_errors.txt' (if the file already exists, errors will be appended to it).

Optionally, a file with locus tags can be given to extract only these proteins with option **-l** (each locus tag in a new line).

## Usage

### Extract amino acid sequences

    perl cds_extractor.pl -i seq_file.gbk -p [-l locus_tags.txt]

### Extract nucleotide sequences

    perl cds_extractor.pl -i seq_file.embl -n [-l locus_tags.txt -u 100 -d 20]

### UNIX loop to extract sequences from all files in the current working directory

    for i in *.[embl|gbk]; do perl cds_extractor.pl -i $i -p [-l locus_tags.txt]; done

## Options for *cds_extractor.pl*

### Mandatory options

* -i, -input

Input RichSeq sequence file including annotation (embl or genbank)

* -p, -protein

Extract protein sequence for each CDS feature, excludes option **-n**

* -n, -nucleotide

Extract nucleotide sequence for each CDS feature, excludes option **-p**

### Optional options

* -h, -help:   Help (perldoc POD)

* -l, -locustag_list

List of locus tags to extract only those

* -u, -upstream

Include given number of flanking nucleotides upstream of each CDS feature, forces option **-n**

* -d, -downstream

Include given number of flanking nucleotides downstream of each CDS feature, forces option **-n**

* -f, -full_header

Include full ID header for downstream *blast_prot_finder.pl* analysis (all IDs include 'g=', 'p=', and 'o=' for parsing regardless of their existence)

## Output

* *\_cds\_aa.fasta

Multi-fasta file of CDS protein sequences

* *\_cds\_nucl.fasta 

Multi-fasta file of CDS DNA sequences

* (locus_tag_errors.txt)

Lists CDS features without locus tags

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Dependencies (not in the core Perl modules)

* BioPerl (tested with version 1.006901)

## Authors/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

* v0.6 (06.06.2013)
    - exit with error if no CDS primary features present in input file (some genbank files are only annotated with 'gene')
    - included Bio::SeqFeatureI's method *spliced-seq* for CDS with split nucleotide sequences (CDS position indicated by 'join')
    - minor changes how the optional list of locus tags is handled
* v0.5 (03.06.2013)
    - included a POD
    - options with Getopt::Long
    - option to alternatively extract nucleotide sequences for CDS features (optionally with upstream and downstream sequences)
    - option to include full fasta ID header for downstream '*blast_prot_finder*' analysis
    - exit with error if the values for two (or more) '/locus_tag' or '/protein_id' tags are not unambiguous
    - print message to STDOUT if and which locus tags were not found in a given locus tag list (option '-l')
* v0.4 (06.02.2013)
    - replace whitespaces of '/product' values with underscores
* v0.3 (06.09.2012)
    - internal CDS counter to use in fasta ID for CDS features without a '/locus_tag', '/protein_id', '/gene', or '/product' tag
    - include also organism (and possible plasmid) information in fasta ID lines
    - give a warning to STDOUT if a CDS feature without a '/locus_tag' is found (but only for the first occurence)
    - additionally, *locus_tag_errors.txt* info file to list all CDSs without locus tags
    - catch errors with *eval* if a tag is missing
* v0.2 (04.09.2012)
    - if a CDS feature does not have a '/locus_tag', then use the value for '/protein_id', '/gene', or '/product' (in this order) in the fasta ID lines of the result file
    - optional extract only CDSs with locus tags given in a file
* v0.1 (24.05.2012)