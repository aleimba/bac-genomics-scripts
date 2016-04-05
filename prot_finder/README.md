prot_finder
===========

`prot_finder.pl` is a script to search for query protein homologs in annotated bacterial genomes with **BLASTP**. A companion bash shell script pipeline is available, `prot_finder_pipe.sh`.

Included is also the script `prot_binary_matrix.pl` to create a presence/absence matrix (e.g. for [**iTOL**](http://itol.embl.de/)) from the `prot_finder.pl` output.

* [Synopsis](#synopsis)
  * [prot_finder synopsis](#prot_finder-synopsis)
  * [prot_binary_matrix synopsis](#prot_binary_matrix-synopsis)
* [Description](#description)
* [Usage](#usage)
  * [prot_finder usage](#prot_finder-usage)
    * [Manual consecutively](#manual-consecutively)
      * [cds_extractor for subject proteins](#cds_extractor-for-subject-proteins)
      * [Legacy BLASTP](#legacy-blastp)
      * [BLASTP plus](#blastp-plus)
      * [prot_finder](#prot_finder-1)
    * [prot_finder_pipe bash script pipeline](#prot_finder_pipe-bash-script-pipeline)
  * [prot_binary_matrix usage](#prot_binary_matrix-usage)
* [Options](#options)
  * [prot_finder.pl options](#prot_finderpl-options)
    * [Mandatory prot_finder.pl options](#mandatory-prot_finderpl-options)
    * [Optional prot_finder.pl options](#optional-prot_finderpl-options)
  * [prot_finder_pipe.sh options](#prot_finder_pipesh-options)
    * [Mandatory prot_finder_pipe.sh options](#mandatory-prot_finder_pipesh-options)
    * [Optional prot_finder_pipe.sh options](#optional-prot_finder_pipesh-options)
  * [prot_binary_matrix.pl options](#prot_binary_matrixpl-options)
* [Output](#output)
  * [cds_extractor.pl output](#cds_extractorpl-output)
  * [prot_finder.pl output](#prot_finderpl-output)
  * [prot_finder_pipe.sh output](#prot_finder_pipesh-output)
  * [prot_binary_matrix.pl output](#prot_binary_matrixpl-output)
* [Dependencies](#dependencies)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Changelog](#changelog)
  * [prot_finder changelog](#prot_finder-changelog)
  * [prot_binary_matrix changelog](#prot_binary_matrix-changelog)

## Synopsis

### prot_finder synopsis

    ./prot_finder_pipe.sh -q query.faa (-s subject.faa|-f (embl|gbk)) > blast_hits.tsv

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa > blast_hits.tsv

### prot_binary_matrix synopsis

    perl prot_binary_matrix.pl blast_hits.tsv > binary_matrix.tsv

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl > binary_matrix.tsv

## Description

The script `prot_finder.pl` is intended to search for homologous proteins in annotated bacterial genomes. For this purpose, a previous [**BLASTP**](http://blast.ncbi.nlm.nih.gov/Blast.cgi), either [legacy or plus](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), needs to be run with query protein sequences against a **BLASTP** database of subject proteins (e.g. all proteins from several *Escherichia coli* genomes).

The script [`cds_extractor.pl`](/cds_extractor) (with options **-p -f**) can be used to create multi-FASTA protein files of all non-pseudo CDS from RichSeq genome files to create the needed subject **BLASTP** database. Present locus tags will be used as FASTA IDs, but see [`cds_extractor.pl`](/cds_extractor) for a description of the format. Query protein sequences for the **BLASTP** need a **unique** FASTA ID.

The **BLASTP** report file (option **-r**), the subject protein multi-FASTA file (option **-s**), and optionally the query protein (multi-)FASTA (option **-q**) file are then given to `prot_finder.pl`. Significant **BLASTP** subject hits are filtered according to the given cutoffs (options **-i**, **-cov_q**, and **-cov_s**) and the result is printed as an informative tab-separated result table to *STDOUT*. To apply global identity/coverage cutoffs to subject hits high-scoring pairs (HSPs) are tiled (see http://www.bioperl.org/wiki/HOWTO:Tiling and http://search.cpan.org/dist/BioPerl/Bio/Search/Hit/GenericHit.pm). Additionally, the subject protein sequences with significant query hits are written to result multi-FASTA files, named according to the respective query FASTA IDs (optionally including the query sequence with option **-q**).

Optionally, [**Clustal Omega**](http://www.clustal.org/omega/) can be called (option **-a** with optional **-p**) to create multiple alignments (FASTA format) for each of the resulting multi-FASTA files. These alignments can be used to calculate phylogenies e.g. with **RAxML** (http://sco.h-its.org/exelixis/software.html) or **MEGA** (http://www.megasoftware.net/).

Run the script [`cds_extractor.pl`](/cds_extractor) (with options **-p -f**) and the **BLASTP** manually or use the bash shell wrapper script `prot_finder_pipe.sh` (see below ['prot_finder_pipe bash script pipeline'](#prot_finder_pipe-bash-script-pipeline)) to execute the whole pipeline including `prot_finder.pl` (with optional option **-q**). For additional options of the pipeline shell script see below ['prot_finder_pipe.sh options'](#prot_finder_pipesh-options). Be aware that some options in `prot_finder_pipe.sh` corresponding to options in `prot_finder.pl` have different names (**-c** instead of **-cov_q**, **-k** instead of **-cov_s**, and **-o** instead of **-p**; also **-f** has a different meaning). If [`cds_extractor.pl`](/cds_extractor) is used in the pipeline (option **-f** of the shell script) the working folder has to contain the annotated bacterial genome subject files (in RichSeq format, e.g. EMBL or GENBANK format). Also, the Perl scripts [`cds_extractor.pl`](/cds_extractor) (only for `prot_finder_pipe.sh` option **-f**) and `prot_finder.pl`have to be either contained in the current working directory or installed in the global *PATH*. **BLASTP** (legacy and/or plus) and **Clustal Omega** binaries have to be installed in global *PATH*, or for **Clustal Omega** you can give the path to the binary with option **-o**. In the pipeline **BLASTP** is run with **disabled** query filtering, locally optimal Smith-Waterman alignments, and increasing the number of database sequences to show alignments to 500 for [**BioPerl**](http://www.bioperl.org) parsing (legacy: **-F F -s T -b 500**, plus: **-seg no -use_sw_tback -num_alignments 500**). The pipeline script ends with the *STDERR* message 'Pipeline finished!', if this is not the case have a look at the log files in the result directory for errors.

At last, the resulting tab-separated table with significant **BLASTP** hits (from `prot_finder.pl` or `prot_finder_pipe.sh`) can be given to the script `prot_binary_matrix.pl`, either as *STDIN* or as a file, to create a presence/absence matrix of the results. See below ['prot_binary_matrix.pl options'](#prot_binary_matrixpl-options) for the `prot_binary_matrix.pl` options. By default a tab-delimited binary presence/absence matrix for query hits per subject organism will be printed to *STDOUT*. Use option **-t** to count all query hits per subject organism, not just the binary presence/absence.

The presence/absence matri(x|ces) can e.g. be loaded into the Interactive Tree Of Life website ([**iTOL**](http://itol.embl.de/)) to associate the data with a phylogenetic tree. [**iTOL**](http://itol.embl.de/) likes individual comma-separated input files, thus use `prot_binary_matrix.pl` options **-s -c** for this purpose. However, the organism names have to have identical names to the leaves of the phylogenetic tree, thus manual adaptation, e.g. in a spreadsheet software (like [**LibreOffice Calc**](https://www.libreoffice.org/discover/calc/)), might be needed. **Careful**, subject organisms without a significant **BLASTP** hit won't be included in the `prot_finder.pl` result table and hence can't be included by `prot_binary_matrix.pl`. If needed add them manually to the result matri(x|ces).

## Usage

### prot_finder usage

#### Manual consecutively

##### cds_extractor for subject proteins

    for file in *.(gbk|embl); do perl cds_extractor.pl -i "$file" -p -f; done
    cat *.faa > subject.faa
    rm !(subject).faa

##### Legacy **BLASTP**

    formatdb -p T -i subject.faa -n prot_finder_db
    blastall -p blastp -d prot_finder_db -i query.faa -o prot_finder.blastp -e 1e-10 -F F -s T -b 500

**or**

##### **BLASTP** plus

    makeblastdb -dbtype prot -in subject.faa -out prot_finder_db
    blastp -db prot_finder_db -query query.faa -out prot_finder.blastp -evalue 1e-10 -seg no -use_sw_tback -num_alignments 500

##### prot_finder

    perl prot_finder.pl -r prot_finder.blastp -s subject.faa -cov_s 80 > blast_hits.tsv
    perl prot_finder.pl -r prot_finder.blastp -s subject.faa -d result_dir -f -q query.faa -i 50 -cov_q 50 -b -a -p ~/bin/clustalo -t 6 > result_dir/blast_hits.tsv

#### prot_finder_pipe bash script pipeline

    ./prot_finder_pipe.sh -q query.faa -s subject.faa > blast_hits.tsv

**or**

    ./prot_finder_pipe.sh -q query.faa -f embl -d result_dir -p legacy -e 0 -t 12 -i 50 -c 50 -k 30 -b -a -o ~/bin/clustalo -m > result_dir/blast_hits.tsv

### prot_binary_matrix usage

    perl prot_binary_matrix.pl -s -d result_dir -t blast_hits.tsv
    perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl -l -c > binary_matrix.csv
    mkdir result_dir && ./prot_finder_pipe.sh -q query.faa -s subject.faa -d result_dir -m | tee result_dir/blast_hits.tsv | perl prot_binary_matrix.pl > binary_matrix.tsv

## Options

### `prot_finder.pl` options

#### Mandatory `prot_finder.pl` options

- **-r**=_str_, **-report**=_str_

    Path to **BLASTP** report/output

- **-s**=_str_, **-subject**=_str_

    Path to subject multi-FASTA protein sequence file (\*.faa) created with [`cds_extractor.pl`](/cds_extractor) (and its options **-p -f**), which was used to create the **BLASTP** database

#### Optional `prot_finder.pl` options

- **-h**, **-help**

    Help (perldoc POD)

- **-d**=_str_, **-dir_result**=_str_

    Path to result folder [default = query identity and coverage cutoffs, './results_i#_cq#']

- **-f**, **-force_dir**

    Force output to an existing result folder, otherwise ask user to remove content of existing folder. Careful, files from a previous analysis might not be overwritten if different to current analysis.

- **-q**=_str_, **-query**=_str_

    Path to query (multi-)FASTA protein sequence file (\*.faa) with **unique** FASTA IDs, which was used as query in the **BLASTP**. Will include each query protein sequence in the respective multi-FASTA 'query-ID_hits.faa' result file.

- **-b**, **-best_hit**

    Give only the best hit (i.e. highest identity) for each subject sequence if a subject has several hits with different queries

- **-i**=_int_, **-ident_cutoff**=_int_

    Query identity cutoff for significant hits (not including gaps), has to be an integer number >= 0 and <= 100 [default = 70]

- **-cov_q**=_int_, **-cov_query_cutoff**=_int_

    Query coverage cutoff, has to be an integer number >= 0 and <= 100 [default = 70]

- **-cov_s**=_int_, **-cov_subject_cutoff**=_int_

    Subject/hit coverage cutoff, has to be an integer >= 0 and <= 100 [default = 0]

- **-a**, **-align_clustalo**

    Call [**Clustal Omega**](http://www.clustal.org/omega/) for multiple alignment of each 'query-ID_hits.faa' result file

- **-p**=_str_, **-path_clustalo**=_str_

    Path to executable **Clustal Omega** binary if not present in global *PATH* variable; requires option **-a**

- **-t**=_int_, **-threads_clustalo**=_int_

    Number of threads for **Clustal Omega** to use; requires option **-a** [default = all processors on system]

- **-v**, **-version**

    Print version number to *STDERR*

### `prot_finder_pipe.sh` options

#### Mandatory `prot_finder_pipe.sh` options

- **-q**=_str_

    Path to query protein (multi-)FASTA file (\*.faa) with **unique** FASTA IDs

- **-f**=_str_

    File extension for files in the **current** working directory to use for [`cds_extractor.pl`](/cds_extractor) (e.g. 'embl' or 'gbk'); excludes shell script option **-s**

**or**

- **-s**=_str_

    Path to subject protein multi-FASTA file (\*.faa) already created with [`cds_extractor.pl`](/cds_extractor) (and its options **-p -f**), will not run [`cds_extractor.pl`](/cds_extractor); excludes shell script option **-f**

#### Optional `prot_finder_pipe.sh` options

- **-h**

    Print usage

- **-d**=_str_

    Path to result folder [default = query identity and coverage cutoffs,'./results_i#_cq#']

- **-p**=(legacy|plus)

    **BLASTP** suite to use [default = plus]

- **-e**=_real_

    E-value for **BLASTP** [default = 1e-10]

- **-t**=_int_

    Number of threads to be used for **BLASTP** and **Clustal Omega** [default = all processors on system]

- **-i**=_int_

    Query identity cutoff for significant hits [default = 70]

- **-c**=_int_

    Query coverage cutoff (corresponds to [`prot_finder.pl` option](#optional-prot_finderpl-options) **-cov_q**) [default = 70]

- **-k**=_int_

    Subject coverage cutoff (corresponds to [`prot_finder.pl` option](#optional-prot_finderpl-options) **-cov_s**) [default = 0]

- **-b**

    Give only the best hit (i.e. highest identity) for each subject sequence

- **-a**

    Multiple alignment of each multi-FASTA result file with [**Clustal Omega**](http://www.clustal.org/omega/)

- **-o**=_str_

    Path to executable **Clustal Omega** binary if not in global *PATH*; requires shell script option **-a** (corresponds to [`prot_finder.pl` option](#optional-prot_finderpl-options) **-p**)

- **-m**

    Clean up all non-essential files, see below ['`prot_finder_pipe.sh` output'](#prot_finder_pipesh-output)

### `prot_binary_matrix.pl` options

- **-h**, **-help**

    Help (perldoc POD)

- **-s**, **-separate**

    Separate presence/absence files for each query protein printed to the result directory [default without **-s** = *STDOUT* matrix for all query proteins combined]

- **-d**=_str_, **-dir_result**=_str_

    Path to result folder, requires option **-s** [default = './binary_matrix_results']

- **-t**, **-total**

    Count total occurrences of query proteins, not just presence/absence binary

- **-c**, **-csv**

    Output matri(x|ces) in comma-separated format (\*.csv) instead of tab-delimited format (\*.tsv)

- **-l**, **-locus_tag**

    Use the locus_tag **prefixes** in the subject_ID column of the `prot_finder.pl` output (instead of the subject_organism column) as organism IDs to associate query hits to organisms. The subject_ID column will include locus_tags if they're annotated for a genome (see the [`cds_extractor.pl`](/cds_extractor) format description). Useful if the [`cds_extractor.pl`](/cds_extractor) output doesn't include strain names for 'o=' in the FASTA IDs, because the prefix of a locus_tag should be unique for a genome (see http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation).

- **-v**, **-version**

    Print version number to *STDERR*

## Output

### `cds_extractor.pl` output

- \*.faa

    Multi-FASTA file(s) of subject CDS protein sequences; will be removed with [`prot_finder_pipe.sh` option](#optional-prot_finder_pipesh-options) **-m**

For optional error files from [`cds_extractor.pl`](/cds_extractor) see its documentation.

### `prot_finder.pl` output

- *STDOUT*

    The resulting tab-delimited output table with the significant subject **BLASTP** hits is printed to *STDOUT*. Redirect (e.g. to a file in the result directory, options **-d -f**) or pipe into another tool as needed (e.g. `prot_binary_matrix.pl`).

- ./results_i#_cq#

    All output files are stored in a result folder

- ./results_i#_cq#/query-ID_hits.faa

    Multi-FASTA protein files of significant subject hits for each query protein (named after the respective query FASTA ID), optionally includes the respective query protein sequence (with option **-q**)

- subject.faa.idx

    Index file of the subject protein file for fast sequence retrieval (can be deleted if no further **BLASTPs** are needed with these subject sequences); will be removed with [`prot_finder_pipe.sh` option](#optional-prot_finder_pipesh-options) **-m**

- (./results_i#_cq#/queries_no_blastp-hits.txt)

    Lists all query sequence IDs without significant subject hits; with option **-b** includes also queries with significant hits but *without* a best blast hit for a subject

- (./results_i#_cq#/clustal_omega.log)

    Optional log file of verbose **Clustal Omega** *STDOUT*/*STDERR* messages; will be removed with [`prot_finder_pipe.sh` option](#optional-prot_finder_pipesh-options) **-m**

- (./results_i#_cq#/query-ID_aln.fasta)

    Optional **Clustal Omega** multiple alignment of each 'query-ID_hits.faa' result file in FASTA alignment format

- (./results_i#_cq#/query-ID_tree.nwk)

    Optional, **Clustal Omega** NJ-guide tree in Newick format

### `prot_finder_pipe.sh` output

In addition to the [`cds_extractor.pl` output](#cds_extractorpl-output) and the [`prot_finder.pl` output](#prot_finderpl-output) the pipeline also creates the following non-essential output files, which will be removed with [`prot_finder_pipe.sh` option](#optional-prot_finder_pipesh-options) **-m**:

- ./results_i#_cq#/cds_extractor.log

    Log file of [`cds_extractor.pl`](/cds_extractor) *STDOUT*/*STDERR* messages (with `prot_finder_pipe.sh` option **-f**)

- ./results_i#_cq#/prot_finder.faa

    Concatenated [`cds_extractor.pl`](/cds_extractor) output files to create the subject **BLASTP** database

- prot_finder_db.phr, prot_finder_db.pin, prot_finder_db.psq

    **BLASTP** database files from the concatenated subject sequences ('prot_finder.faa')

- formatdb.log, error.log **or** ./results_i#_cq#/makeblastdb.log

    Legacy **BLASTP** or **BLASTP+** log files

- ./results_i#_cq#/prot_finder.blastp

    **BLASTP** report/output

- ./results_i#_cq#/prot_finder.log

    Log file of `prot_finder.pl` *STDOUT*/*STDERR* messages

### `prot_binary_matrix.pl` output

- *STDOUT*

    The resulting presence/absence matrix is printed to *STDOUT* without option **-s**. Redirect or pipe into another tool as needed.

- (./binary_matrix_results)

    Separate query presence/absence files are stored in a result folder with option **-s**

- (./binary_matrix_results/query-ID_binary_matrix.(tsv|csv))

    Separate query presence/absence files with option **-s**

## Dependencies

- [**BioPerl**](http://www.bioperl.org) (tested version 1.006923)
- [`cds_extractor`](/cds_extractor) (tested version 0.7.1)
- [**BLASTP**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - legacy **BLASTP** (tested version 2.2.26)
    - **BLASTP+** (tested version 2.2.28+)
- [**Clustal Omega**](http://www.clustal.org/omega/) (tested version 1.2.1)

## Run environment

The scripts run under UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

### prot_finder changelog

* v0.7.1 (05.04.2016)
    * bug fix: significant but non-best blast hits with option **-b** now listed in *queries_no_blastp-hits.txt*
* v0.7 (23.11.2015)
    * changed script name from `blast_prot_finder.pl` to `prot_finder.pl`
    * fixed bug introduced in v0.6 with a `seek`, because option **-query** didn't pull query sequences anymore
    * included version switch
    * included 'use autodie;' pragma
    * included `pod2usage` with Pod::Usage and `pod2usage`-die for Getopt::Long call
    * tab-delimited output table with filtered **BLASTP** hits now printed to *STDOUT* instead of file and can be used as *STDIN* for `prot_binary_matrix.pl`
    * result files now created in result directory to unclutter working dir (new options **-dir_result** and **-force_dir**) and replaced subroutine 'file_exist' with 'empty_dir' (also no *STDOUT* message which files were created anymore)
    * new output file *queries_no_blastp-hits.txt* to list all queries without **BLASTP** hits
    * new output file *clustal_omega.log* to log **Clustal Omega** verbose output
    * additional new options, **-path_clustalo** and **-threads_clustalo**
    * major code changes/restructuring and additions:
        * POD syntax changes and additions for new code
        * Perl syntax changes and removing code redundancies
        * including script run status messages
        * changes to subroutine 'split_fasta_header' to work more robustly with `split` instead of a regex and adapted to `cds_extractor.pl` v0.7+ output format
    * included several additional option and input checks
    * replaced simple `blast_prot_finder_legacy.sh` with more elaborate bash script pipeline `prot_finder_pipe.sh`
    * changed some output column names and moved the subject_ID column before the subject_gene column
* v0.6 (10.06.2013)
    * included BioPerls 'frac\*' methods, which include hsp tiling, to correct bug in query coverage and identity calculations for **whole** hits
    * corrected bug to use **whole** hit e-value and not just first hsp e-value of a significant hit
    * option **-cov_subject** for subject/hit coverage cutoff (not only query coverage cutoff)
    * more info in POD
* v0.5 (14.05.2013)
    * optionally include query proteins in result subject multi-FASTA files
    * new option **-best_hit** to include only best **BLASTP** hit (highest identity) for each subject locus_tag, if a subject protein has hits to several query proteins
    * print queries with no subject **BLASTP** hits to *STDOUT*
* v0.4 (24.01.2013)
    * corrected bug in hash structure to store **BLASTP** hits, query acc/IDs (keys) and array reference of subject locus_tags (values), as several queries can have the same subject locus_tag as hit
    * include 'subject protein_function' in 'blast_hits.txt' output
* v0.3 (12.09.2012)
    * **original** script name `blast_prot_finder.pl`
    * **BLASTP** hits are stored in a hash with subject locus_tags (keys) and query accessions/IDs (values)

### prot_binary_matrix changelog

* v0.6 (23.11.2015)
    * adapted to `prot_finder.pl` v0.7 output
    * included a POD, `pod2usage` with Pod::Usage and `pod2usage`-die for Getopt::Long call
    * Perl syntax changes and some simpler loop structures
    * Removed option **-input** and instead accept `prot_finder.pl` output as *STDIN* or file as argument
    * Removed option **-c|-collective** (actually repurposed, see below), option **-separate** to indicate separate output is enough
    * **Without** option **-separate** result matrix is now printed to *STDOUT*
    * **With** option **-separate** result matrices now printed to result directory (name optionally given with new option **-dir_result**)
    * New option **-locus_tag** to use the locus_tag prefix in column subject_ID of the `prot_finder.pl` output as ID, which is in most cases the locus_tag, instead of the subject_organism; controls if the locus_tags follow NCBI standards
    * Output now default tab-separated format (*tsv*), optionally with new option **-csv** comma-separated (*csv*)
    * Removed subroutine 'result' and statement which result files have been created
* v0.5 (05.03.2014)
    * option **-total** to count all occurences of query proteins within one organism (paralogs), instead of just binary presence/absence
    * enforce mandatory options
    * changed script name to 'prot_binary_matrix.pl'
    * version switch
    * included 'use autodie'
    * options with Getopt::Long
    * changed usage to HERE document
* v0.4 (29.04.2013)
    * adapted to new *prot_finders* 'blast_hits.txt' layout, which includes an additional column for 'subject protein_function'
    * changed script name to 'blast_binary_matrix.pl'
* v0.3 (15.02.2013)
    * status of created result files in *STDOUT* for option **-s**
* v0.2 (21.12.2012)
    * prot_finder output file 'blast_hits.txt' doesn't have to be ordered by query protein accessions/IDs anymore
* v0.1 (25.10.2012)
    * **original** script name 'blastp_iTOL_binary.pl'
