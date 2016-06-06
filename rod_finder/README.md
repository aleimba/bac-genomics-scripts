rod_finder
==========

A script to find regions of difference (RODs) between a query genome and reference genome(s).

## Synopsis

    blast_rod_finder_legacy.sh subject.fasta query.fasta query.[embl|gbk|fasta] 5000

or

    perl blast_rod_finder.pl -q query.embl -r blastn.out -m 2000

## Description

This script is intended to identify RODs between a nucleotide query and a nucleotide subject/reference sequence. In order to do so, a *blastn* (http://blast.ncbi.nlm.nih.gov/Blast.cgi) needs to be performed beforehand with the query and the subject sequences (see also *blast_rod_finder_legacy.sh* below). *blast_rod_finder.pl* is mainly designed to work with bacterial genomes, while a query genome can be blasted against several subject sequences to detect RODs over a number of references. Although the results are optimized towards a complete query genome, both the reference(s) as well as the query can be used in draft form. To create artificial genomes via concatenation use *cat_seq.pl* or the EMBOSS application union (http://emboss.sourceforge.net/).

The *blastn* report file, the query sequence file (preferably in RichSeq format, see below) and a minimum size for ROD detection have to be provided. Subsequently, RODs are summarized in a tab-separated summary file, a gff3 (usable e.g. in Artemis/DNAPlotter, http://www.sanger.ac.uk/resources/software/artemis/) and a BRIG (BLAST Ring Image Generator, http://brig.sourceforge.net/) output file. Nucleotide sequences of each ROD are written to a multi-fasta file.

The query sequence can be provided in RichSeq format (embl or genbank), but has to correspond to the fasta file used in querying the BLAST database (the accession numbers have to correspond to the fasta headers). Use *seq_format-converter.pl* to create a corresponding fasta file from embl|genbank files for *blastn* if needed. With RichSeq query files additional info is given in the result summary and the amino acid sequences of all non-pseudo CDSs, which are contained or overlap a ROD, are written to a result file. Furthermore, all detected RODs are saved in individual sequence files in the corresponding query sequence format.

Run *blastn* and the script *blast_rod_finder.pl* consecutively manually or use the bash shell wrapper script *blast_rod_finder_legacy.sh* (see usage below) to perform the pipeline with one command. The same folder has to contain the subject fasta file(s), the query fasta file, optionally the query RichSeq file and the script *blast_rod_finder.pl*! *blastn* is run **without** filtering of query sequences ('-F F') and an evalue cutoff of '2e-11' is set.

## Usage

### 1.) Manual consecutively

#### 1.1.) *blastn*

    formatdb -p F -i subject.fasta -n blast_db
    blastall -p blastn -d blast_db -i query.fasta -o blastn.out -e 2e-11 -F F

#### 1.2.) *blast_rod_finder.pl*

    perl blast_rod_finder.pl -q query.[embl|gbk|fasta] -r blastn.out -m 5000

### 2.) With one command: *blast_rod_finder_legacy.sh* pipeline

    blast_rod_finder_legacy.sh subject.fasta query.fasta query.[embl|gbk|fasta] 5000

## Options for *blast_rod_finder.pl*

### Mandatory options

* -m, -min

Minimum size of RODs that are reported

* -q, -query

Query sequence file [fasta, embl, or genbank format]

* -r, -report

*blastn* report/output file

### Optional options

* -h, -help:   Help (perldoc POD)

## Output

### a.) *blast_rod_finder_legacy.sh* or *blastn*

* *blastn* database files for subject sequence(s)

\*.nhr, \*.nin, \*.nsq

* *blastn* report

Text file named 'blastn.out'

### b.) *blast_rod_finder.pl*

* ./results

All output files are stored in this result folder

* rod_summary.txt

Summary of detected ROD regions (for embl/genbank queries includes annotation), tab-separated

* rod.gff

GFF3 file with ROD coordinates to use in Artemis/DNAPlotter etc.

* rod_BRIG.txt

ROD coordinates to use in BRIG (BLAST Ring Image Generator), tab-separated

* rod_seq.fasta

Nucleotide sequences of ROD regions (>ROD# size start..stop), multi-fasta

* rod_aa_fasta.txt

Only present if query is in RichSeq format. Amino acid sequences of all CDSs that are contained in or overlap a ROD region in multi-fasta format (>locus_tag gene product). RODs are seperated in the file via '\~\~' (\~\~ROD# size start..stop).

* ROD#.[embl|gbk]

Only present if query is in RichSeq format. Each identified ROD is written to an individual sequence file (in the same format as the query).

## Run environment

The Perl script runs under Windows and UNIX flavors, the bash-shell script of course only under UNIX.

## Dependencies (not in the core Perl modules)

* Legacy blast (tested version blastall 2.2.18)
* BioPerl (tested with version 1.006901)

## Authors/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

David Studholme (original code; D[dot]J[dot]Studholme[at]exeter[dot]ac[dot]uk; University of Exeter)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.4 (13.02.2013)
    - included a POD
    - options with Getopt::Long
    - results directory for output files
    - include accession number column for multi-sequence files in 'rod_summary.txt'
    - include locus_tags (or alternatively gene, product, note ...) in 'rod_summary.txt'
    - feature positions according to leading or lagging strand
    - indicate if a primary feature overlaps ROD boundaries
    - output each ROD in the query RichSeq format with BioPerl's Bio::SeqUtils
* v0.3 (23.11.2011)
    - status messages with autoflush
    - BRIG output file
    - extended primary tag output for RODs (in addition to CDS): tRNA, rRNA, tmRNA, ncRNA, misc_RNA, repeat_region, misc_binding, and mobile_element
* v0.1 (07.11.2011)
