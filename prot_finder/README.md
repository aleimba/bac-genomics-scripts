prot_finder
===========

`prot_finder.pl` is a script to search for homologous proteins in annotated bacterial genomes via *blastp*.

Included is the script `prot_binary_matrix.pl` to create a presence/abscence matrix (e.g. for iTOL) from the `prot_finder.pl` output.

## Synopsis

### *prot_finder*

    blast_prot_finder_legacy.sh subject_file-extension query_prot.fasta ident_cutoff cov_q_cutoff

or

    perl blast_prot_finder.pl -s subject_prot.fasta -r blastp-report.out -b

### *prot_binary_matrix*

    prot_binary_matrix.pl -i blast_hits.txt -c

## Description

This script is intended to search for homologous proteins in annotated bacterial genomes. Therefore, a previous *blastp* (http://blast.ncbi.nlm.nih.gov/Blast.cgi) needs to be run with known query protein sequences against a blast database of subject proteins (e.g. all proteins from several *E. coli* genomes).

For this purpose the script `cds_extractor.pl` (with option **-f**) can be used to create multi-fasta protein files of all non-pseudo CDS from genome sequence files to create the needed blast database. The *blastp* report file, the subject multi-fasta file, and optionally the query protein fasta file are then given to `blast_prot_finder.pl`. Significant blast hits are filtered and the result is written to a tab-separated file. The subject hits are also concatenated in a multi-fasta file for each query sequence (optionally including the query sequence).

Optionally, Clustal Omega (http://www.clustal.org/omega/) can be called (has to be in the **$PATH** or change variable `$clustal_call`) to create an alignment (fasta format) for each of the concatenated multi-fasta files. The alignments can then be used to calculate phylogenies. Use e.g. SeaView (http://pbil.univ-lyon1.fr/software/seaview.html) or `aln_format-converter.pl` to convert the alignment format to Phylip for RAxML (http://sco.h-its.org/exelixis/software.html). MEGA (http://www.megasoftware.net/) can work directly with the Clustal fasta format.

Run the script `cds_extractor.pl` (with option **-f**) and the *blastp* manually or use the bash shell wrapper script `blast_prot_finder_legacy.sh` (see usage below) to execute the whole pipeline (including `blast_prot_finder.pl` with optional options **-q**, **-b**, and **-a**) consecutively in one folder. For this purpose the same folder has to contain the annotated bacterial genome subject files (in RichSeq format, e.g. embl or genbank), the query protein fasta, and the scripts `cds_extractor.pl` and `blast_prot_finder.pl`! *blastp* is run **without** filtering of query sequences (**-F F**), an evalue cutoff of **1e-10**, local optimal Smith-Waterman alignments (**-s T**), and increasing the number of database sequences to show alignments to 500 (**-b 500**, to adjust it to the alignment summary for BioPerl).

Additionally, the result file 'blast_hits.txt' can be given to the script `prot_binary_matrix.pl` to create a presence/abscence binary matrix of the results, '\*binary_matrix.csv'. This comma-separated file can e.g. be loaded into iTOL (http://itol.embl.de/) to associate the data with a phylogenetic tree. However, the organism names have to have the same names as in the tree file, thus manual adaptation, e.g. in Excel, might be needed. **Careful**, organisms that didn't have a *blastp* hit, are obviously not present in 'blast_hits.txt', and hence can't be included by `prot_binary_matrix.pl`. Thus, add them manually to the result file(s).

## Usage

### *prot_finder*

#### 1.) Manual consecutively

##### 1.1.) `cds_extractor.pl` and concatenate subject proteins for *blastp* database

    for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i -p -f; done
    cat *_cds_aa.fasta > subject_prot.fasta
    rm -f *_cds_aa.fasta

##### 1.1.) Legacy *blastp*

    formatdb -p T -i subject_prot.fasta -n blast_finder
    blastall -p blastp -d blast_finder -i query_prot.fasta -o blastp.out -e 1e-10 -F F -s T -b 500

##### 1.2.) `blast_prot_finder.pl`

    perl blast_prot_finder.pl -q query_prot.fasta -s subject_prot.fasta -r blastp.out -i 50 -cov_q 50 -b -a

#### 2.) With one command: `blast_prot_finder_legacy.sh` wrapper script

    blast_prot_finder_legacy.sh subject_file-extension query_prot.fasta ident_cutoff cov_q_cutoff

### *prot_binary_matrix*

    perl prot_binary_matrix.pl -i blast_hits.txt [-c|-s] (-t)

## Options

### *prot_finder.pl*

#### Mandatory options

* -r, -report

*blastp* report/output file

* -s, -subject

Subject sequence file [fasta format]

#### Optional options

* -h, -help:   Help (perldoc POD)

* -q, -query

Query sequence file (to include each query protein sequence in the respective 'query-acc_hits.fasta' result file) [fasta format]

* -b, -best

Give only the best hit (i.e. highest identity) for each subject locus tag, if a locus tag has several hits with different queries

* -i, -ident

Query identity cutoff (not including gaps) [default 70]

* -cov_q, -cov_query

Query coverage cutoff [default 70]

* -cov_s, -cov_subject

Subject/hit coverage cutoff [default 0]

* -a, -align

Call Clustal Omega for alignment

### *prot_binary_matrix.pl*

#### Mandatory options

* -i, -input

Input 'blast_hits.txt' file from `blast_prot_finder.pl`

* -c, -collective

Collective presence/abscence file for all query proteins combined [default, excludes **-s**]

* -s, -separate

Separate presence/abscence files for each query protein [excludes **-c**]

#### Optional options

* -h, -help:   Print usage

* -t, -total

Count total occurences of query proteins not just presence/abscence binary

* -v, -version:   Print version number

## Output

### *prot_finder.pl*

#### a.) `cds_extractor.pl`

* blast_finder_aa.fasta

Concatenated subject proteins from all non-pseudo CDSs of all subject organisms

#### b.) `blast_prot_finder_legacy.sh` or *blastp*

* *blastp* database files for subject sequences

\*.phr, \*.pin, \*.psq

* *blastp* report

Text file named 'blastp.out'

#### c.) `blast_prot_finder.pl`

* blast_hits.txt

Contains a list of the filtered *blastp* hits, tab-separated

* query-acc_hits.fasta

Multi-fasta protein files of subject hits for each query protein (with acc#), optionally includes the respective query protein sequence

* \*.idx.dir and \*.idx.pag

Index files of the subject protein file for fast sequence retrieval (can be deleted if no further *blastp* analyses are needed with these subject sequences)

* (query-acc_aln.fasta)

Optional, Clustal Omega alignment of subject hits for each query [fasta format]

* (query-acc_tree.nwk)

Optional, Clustal Omega NJ-guide tree in Newick format

### *prot_binary_matrix.pl*

* \*binary_matrix.csv

Comma-separated presence/abscence matrix of the query proteins for the subject organisms. Option **-c** results in one collective file, option **-s** in separate files for each query protein.

## Run environment

The Perl scripts run under Windows and UNIX flavors, the bash-shell script of course only under UNIX.

## Dependencies (not in the core Perl modules)

* BioPerl (tested version 1.006901)
* Legacy blast (tested version blastall 2.2.18)
* Clustal Omega (tested version 1.1.0)

## Authors/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

### *prot_finder*

* v0.6 (10.06.2013)
    - included BioPerls 'frac\*' methods, which include hsp tiling, to correct bug in query coverage and identity calculations for **whole** hits
    - corrected bug to use **whole** hit e-value and not just first hsp e-value of a significant hit
    - option **-cov_subject** for subject/hit coverage cutoff (not only query coverage cutoff)
    - more info in POD
* v0.5 (14.05.2013)
    - optionally include query proteins in result subject multi-fasta files
    - new option **-best** to include only best blastp hit (highest identity) for each subject locus_tag, if a subject protein has hits to several query proteins
    - print queries with no subject blast hits to STDOUT
* v0.4 (24.01.2013)
    - corrected bug in hash structure to store blastp hits, query acc (keys) and array reference of subject locus_tags (values), as several queries can have the same subject locus_tag as hit
    - include 'subject protein_function' in 'blast_hits.txt' output
* v0.3 (12.09.2012)
    - **original** script name 'blast_prot_finder'
    - blastp hits are stored in a hash with subject locus_tags (keys) and query accessions (values)


### *prot_binary_matrix*

* v0.5 (05.03.2014)
    - option **-total** to count all occurences of query proteins within one organism (paralogs), instead of just binary presence/abscence
    - enforce mandatory options
    - changed script name to 'prot_binary_matrix.pl'
    - version switch
    - included 'use autodie'
    - options with Getopt::Long
    - changed usage to HERE document
* v0.4 (29.04.2013)
    - adapted to new *prot_finders* 'blast_hits.txt' layout, which includes an additional column for 'subject protein_function'
    - changed script name to 'blast_binary_matrix.pl'
* v0.3 (15.02.2013)
    - status of created result files in STDOUT for option **-s**
* v0.2 (21.12.2012)
    - prot_finder output file 'blast_hits.txt' doesn't have to be ordered by query protein accessions anymore
* v0.1 (25.10.2012)
    - **original** script name 'blastp_iTOL_binary.pl'
