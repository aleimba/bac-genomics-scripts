prot_finder
===========

`prot_finder.pl` is a script to search for query protein homologs in annotated bacterial genomes with **BLASTP**. A companion bash shell script pipeline is available, `prot_finder_pipe.sh`.

Included is also the script `prot_binary_matrix.pl` to create a binary presence/absence matrix (e.g. for [**iTOL**](http://itol.embl.de/)) from the `prot_finder.pl` output. Additionally, two downstream scripts are provided to wrangle these binary presence/absence matrix: `transpose_matrix.pl` to transpose a delimited TEXT matrix and `binary_group_stats.pl` to get overall presence/absence statistics for groups of columns in a delimited binary TEXT matrix (in the style of [`po2group_stats.pl`](/po2group_stats)).

* [Synopsis](#synopsis)
  * [prot_finder synopsis](#prot_finder-synopsis)
  * [prot_binary_matrix synopsis](#prot_binary_matrix-synopsis)
  * [transpose_matrix synopsis](#transpose_matrix-synopsis)
  * [binary_group_stats synopsis](#binary_group_stats-synopsis)
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
  * [transpose_matrix usage](#transpose_matrix-usage)
  * [binary_group_stats usage](#binary_group_stats-usage)
* [Options](#options)
  * [prot_finder.pl options](#prot_finderpl-options)
    * [Mandatory prot_finder.pl options](#mandatory-prot_finderpl-options)
    * [Optional prot_finder.pl options](#optional-prot_finderpl-options)
  * [prot_finder_pipe.sh options](#prot_finder_pipesh-options)
    * [Mandatory prot_finder_pipe.sh options](#mandatory-prot_finder_pipesh-options)
    * [Optional prot_finder_pipe.sh options](#optional-prot_finder_pipesh-options)
  * [prot_binary_matrix.pl options](#prot_binary_matrixpl-options)
  * [transpose_matrix.pl options](#transpose_matrixpl-options)
  * [binary_group_stats.pl options](#binary_group_statspl-options)
    * [Mandatory binary_group_stats.pl options](#mandatory-binary_group_statspl-options)
    * [Optional binary_group_stats.pl options](#optional-binary_group_statspl-options)
* [Output](#output)
  * [cds_extractor.pl output](#cds_extractorpl-output)
  * [prot_finder.pl output](#prot_finderpl-output)
  * [prot_finder_pipe.sh output](#prot_finder_pipesh-output)
  * [prot_binary_matrix.pl output](#prot_binary_matrixpl-output)
  * [transpose_matrix.pl output](#transpose_matrixpl-output)
  * [binary_group_stats.pl output](#binary_group_statspl-output)
* [Dependencies](#dependencies)
  * [`prot_finder.pl`/`prot_finder_pipe.sh` dependencies](#prot_finderplprot_finder_pipesh-dependencies)
  * [`binary_group_stats.pl` dependencies](#binary_group_statspl-dependencies)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Acknowledgements](#acknowledgements)
* [Changelog](#changelog)
  * [prot_finder changelog](#prot_finder-changelog)
  * [prot_binary_matrix changelog](#prot_binary_matrix-changelog)
  * [transpose_matrix changelog](#transpose_matrix-changelog)
  * [binary_group_stats changelog](#binary_group_stats-changelog)

## Synopsis

### prot_finder synopsis

    ./prot_finder_pipe.sh -q query.faa (-s subject.faa|-f (embl|gbk)) > blast_hits.tsv

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa > blast_hits.tsv

### prot_binary_matrix synopsis

    perl prot_binary_matrix.pl blast_hits.tsv > binary_matrix.tsv

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl > binary_matrix.tsv

### transpose_matrix synopsis

    perl transpose_matrix.pl input_matrix.tsv > input_matrix_transposed.tsv

**or**

    perl prot_binary_matrix.pl blast_hits.tsv | perl transpose_matrix.pl > binary_matrix_transposed.tsv

### binary_group_stats synopsis

    perl binary_group_stats.pl -i binary_matrix.tsv -g group_file.tsv -p > overall_stats.tsv

## Description

The script `prot_finder.pl` is intended to search for homologous proteins in annotated bacterial genomes. For this purpose, a previous [**BLASTP**](http://blast.ncbi.nlm.nih.gov/Blast.cgi), either [legacy or plus](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), needs to be run with query protein sequences against a **BLASTP** database of subject proteins (e.g. all proteins from several *Escherichia coli* genomes).

The script [`cds_extractor.pl`](/cds_extractor) (with options **-p -f**) can be used to create multi-FASTA protein files of all non-pseudo CDS from RichSeq genome files to create the needed subject **BLASTP** database. Present locus tags will be used as FASTA IDs, but see [`cds_extractor.pl`](/cds_extractor) for a description of the format. Query protein sequences for the **BLASTP** need a **unique** FASTA ID.

The **BLASTP** report file (option **-r**), the subject protein multi-FASTA file (option **-s**), and optionally the query protein (multi-)FASTA (option **-q**) file are then given to `prot_finder.pl`. Significant **BLASTP** subject hits are filtered according to the given cutoffs (options **-i**, **-cov_q**, and **-cov_s**) and the result is printed as an informative tab-separated result table to *STDOUT*. To apply global identity/coverage cutoffs to subject hits high-scoring pairs (HSPs) are tiled (see http://www.bioperl.org/wiki/HOWTO:Tiling and http://search.cpan.org/dist/BioPerl/Bio/Search/Hit/GenericHit.pm). Additionally, the subject protein sequences with significant query hits are written to result multi-FASTA files, named according to the respective query FASTA IDs (optionally including the query sequence with option **-q**).

Optionally, [**Clustal Omega**](http://www.clustal.org/omega/) can be called (option **-a** with optional **-p**) to create multiple alignments (FASTA format) for each of the resulting multi-FASTA files. These alignments can be used to calculate phylogenies e.g. with **RAxML** (http://sco.h-its.org/exelixis/software.html) or **MEGA** (http://www.megasoftware.net/).

Run the script [`cds_extractor.pl`](/cds_extractor) (with options **-p -f**) and the **BLASTP** manually or use the bash shell wrapper script `prot_finder_pipe.sh` (see below ['prot_finder_pipe bash script pipeline'](#prot_finder_pipe-bash-script-pipeline)) to execute the whole pipeline including `prot_finder.pl` (with optional option **-q**). For additional options of the pipeline shell script see below ['prot_finder_pipe.sh options'](#prot_finder_pipesh-options). Be aware that some options in `prot_finder_pipe.sh` corresponding to options in `prot_finder.pl` have different names (**-c** instead of **-cov_q**, **-k** instead of **-cov_s**, and **-o** instead of **-p**; also **-f** has a different meaning). If [`cds_extractor.pl`](/cds_extractor) is used in the pipeline (option **-f** of the shell script) the working folder has to contain the annotated bacterial genome subject files (in RichSeq format, e.g. EMBL or GENBANK format). Also, the Perl scripts [`cds_extractor.pl`](/cds_extractor) (only for `prot_finder_pipe.sh` option **-f**) and `prot_finder.pl`have to be either contained in the current working directory or installed in the global *PATH*. **BLASTP** (legacy and/or plus) and **Clustal Omega** binaries have to be installed in global *PATH*, or for **Clustal Omega** you can give the path to the binary with option **-o**. In the pipeline **BLASTP** is run with **disabled** query filtering, locally optimal Smith-Waterman alignments, and increasing the number of database sequences to show alignments to 500 for [**BioPerl**](http://www.bioperl.org) parsing (legacy: **-F F -s T -b 500**, plus: **-seg no -use_sw_tback -num_alignments 500**). The pipeline script ends with the *STDERR* message 'Pipeline finished!', if this is not the case have a look at the log files in the result directory for errors.

The resulting tab-separated table with significant **BLASTP** hits (from `prot_finder.pl` or `prot_finder_pipe.sh`) can be given to the script `prot_binary_matrix.pl`, either as *STDIN* or as a file, to create a presence/absence matrix of the results. See below ['prot_binary_matrix.pl options'](#prot_binary_matrixpl-options) for the `prot_binary_matrix.pl` options. By default a tab-delimited binary presence/absence matrix for query hits per subject organism will be printed to *STDOUT*. Use option **-t** to count all query hits per subject organism, not just the binary presence/absence. This presence/absence matrix can be given to the script `transpose_matrix.pl`, either as *STDIN* or as a file, to transpose the matrix, i.e. rows will become columns and columns rows. Actually, `transpose_matrix.pl` can be used to transpose any delimited TEXT matrix (see below ['transpose_matrix.pl options'](#transpose_matrixpl-options)).

The presence/absence matri(x|ces) can e.g. be loaded into the Interactive Tree Of Life website ([**iTOL**](http://itol.embl.de/)) to associate the data with a phylogenetic tree. [**iTOL**](http://itol.embl.de/) likes individual comma-separated input files, thus use `prot_binary_matrix.pl` options **-s -c** for this purpose. However, the organism names have to have identical names to the leaves of the phylogenetic tree, thus manual adaptation, e.g. in a spreadsheet software (like [**LibreOffice Calc**](https://www.libreoffice.org/discover/calc/)), might be needed. **Careful**, subject organisms without a significant **BLASTP** hit won't be included in the `prot_finder.pl` result table and hence can't be included by `prot_binary_matrix.pl`. If needed add them manually to the result matri(x|ces).

At last, script `binary_group_stats.pl` can be used to categorize columns of the binary presence/absence matrix from `prot_binary_matrix.pl` according to group affiliations. `binary_group_stats.pl` is based upon [`po2group_stats.pl`](/po2group_stats), which does the same thing for genomes in an ortholog/paralog output matrix from a [Proteinortho5](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) calculation. Actually, `binary_group_stats.pl` can work with any delimited TEXT **binary** matrix (option **-i**). But, all fields of the binary matrix need to be filled with either a **0** indicating absence or a **1** indicating presence, i.e. all rows need to have the same number of columns. Use option **-d** to set the delimiter of the input matrix, default is set to tab-delimited/separated matrices.

Also, column headers in the first row and row headers in the first column are **mandatory** for the input binary matrix. Only alphanumeric (a-z, A-Z, 0-9), underscore (_), dash (-), and period (.) characters are allowed for the **column headers** and **group names** in the group file (option **-g**) to avoid downstream problems with the operating/file system. As a consequence, also no whitespaces are allowed in these! Additionally, **column headers**, **row headers**, and **group names** need to be **unique**.

The group affiliations of the columns are intended to get overall presence/absence statistics for groups of columns and not simply single columns of the matrix. Percentage inclusion (option **-cut_i**) and exclusion (option **-cut_e**) cutoffs can be set to define how strict the presence/absence of column groups within a row are defined. Of course groups can also hold only single column headers to get single column statistics. Group affiliations are defined in a mandatory **tab-delimited** group input file (option **-g**), including the column headers from the input binary matrix, with **minimal two** and **maximal four** groups.

See the *README.md* of [`po2group_stats.pl`](/po2group_stats) for an explanation of the logic behind the categorization and the resulting group binary matrix and venn diagram of `binary_group_stats.pl` (and of course its [options](#binary_group_statspl-options) and [output](#binary_group_statspl-output)).

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

**or**

    perl prot_finder.pl -r prot_finder.blastp -s subject.faa -d result_dir -f -q query.faa -i 50 -cov_q 50 -b -a -p ~/bin/clustalo -t 6 > result_dir/blast_hits.tsv

#### prot_finder_pipe bash script pipeline

    ./prot_finder_pipe.sh -q query.faa -s subject.faa > blast_hits.tsv

**or**

    ./prot_finder_pipe.sh -q query.faa -f embl -d result_dir -p legacy -e 0 -t 12 -i 50 -c 50 -k 30 -b -a -o ~/bin/clustalo -m > result_dir/blast_hits.tsv

### prot_binary_matrix usage

    perl prot_binary_matrix.pl -s -d result_dir -t blast_hits.tsv

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl -l -c > binary_matrix.csv

**or**

    mkdir result_dir && ./prot_finder_pipe.sh -q query.faa -s subject.faa -d result_dir -m | tee result_dir/blast_hits.tsv | perl prot_binary_matrix.pl > binary_matrix.tsv

### transpose_matrix usage

    perl transpose_matrix.pl -d ' ' -e NA input_matrix_space-delimit.txt > input_matrix_space-delimit_transposed.txt

**or**

    perl prot_finder.pl -r report.blastp -s subject.faa | perl prot_binary_matrix.pl -l -c | perl transpose_matrix.pl -d , > binary_matrix_transposed.csv

**or**

    mkdir result_dir && ./prot_finder_pipe.sh -q query.faa -s subject.faa -d result_dir -m | tee result_dir/blast_hits.tsv | perl prot_binary_matrix.pl | tee result_dir/binary_matrix.tsv | perl transpose_matrix.pl > result_dir/binary_matrix_transposed.tsv

Transpose all matrices in a folder:

    for matrix in *.tsv; do perl transpose_matrix.pl "$matrix" > "${matrix%.*}_transposed.tsv"; done

### binary_group_stats usage

    perl binary_group_stats.pl -i binary_matrix_transposed.csv -g group_file.tsv -d , -r result_dir -cut_i 0.7 -cut_e 0.2 -b -p -co -s -u -a > overall_stats.tsv

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

### `transpose_matrix.pl` options

- **-h**, **-help**

    Help (perldoc POD)

- **-d**=_str_, **-delimiter**=_str_

    Set delimiter of input and output matrix (e.g. comma ',', single space ' ' etc.) [default = tab-delimited/separated]

- **-e**=_str_, **-empty**=_str_

    Fill empty cells of the input matrix with a value in the transposed matrix (e.g. 'NA', '0' etc.)

- **-v**, **-version**

    Print version number to *STDERR*

### `binary_group_stats.pl` options

#### Mandatory `binary_group_stats.pl` options

- **-i**=*str*, **-input**=*str*

    Input delimited TEXT binary matrix (e.g. *.tsv, *.csv, or *.txt), see option **-d**

- **-g**=*str*, **-groups_file**=*str*

    Tab-delimited file with group affiliation for the columns from the input binary matrix with **minimal two** and **maximal four** groups (easiest to create in a spreadsheet software and save in tab-separated format). **All** column headers from the input binary matrix need to be included. Column headers and group names can only include alphanumeric (a-z, A-Z, 0-9), underscore (_), dash (-), and period (.) characters (no whitespaces allowed either). Example format with two column headers in group A, three in group B and D, and one in group C:

    group\_A&emsp;group\_B&emsp;group\_C&emsp;group\_D<br>
    column\_header1&emsp;column\_header9&emsp;column\_header3&emsp;column\_header8<br>
    column\_header7&emsp;column\_header6&emsp;&emsp;column\_header5<br>
    &emsp;column\_header4&emsp;&emsp;column\_header2

#### Optional `binary_group_stats.pl` options

- **-h**, **-help**

    Help (perldoc POD)

- **-d**=*str*, **-delimiter**=*str*

    Set delimiter of input binary matrix (e.g. comma ',', single space ' ' etc.) [default = tab-delimited/separated]

- **-r**=*str*, **-result\_dir**=*str*

    Path to result folder \[default = inclusion and exclusion percentage cutoff, './results\_i#\_e#'\]

- **-cut\_i**=*float*, **-cut\_inclusion**=*float*

    Percentage inclusion cutoff for column presence counts in a group per row, has to be > 0 and <= 1. Cutoff will be rounded according to the column header number in each group and has to be > the rounded exclusion cutoff in this group. \[default = 0.9\]

- **-cut\_e**=*float*, **-cut\_exclusion**=*float*

    Percentage exclusion cutoff, has to be >= 0 and < 1. Rounded cutoff has to be < rounded inclusion cutoff. \[default = 0.1\]

- **-b**, **-binary\_group\_matrix**

    Print a group binary matrix with the presence/absence column group results according to the cutoffs (excluding 'unspecific' category rows)

- **-p**, **-plot\_venn**

    Plot venn diagram from the group binary matrix (except 'unspecific' and 'underrepresented' categories)  with function `venn` from **R** package **gplots**, requires option **-b**

- **-co**, **-core_strict**

    Include 'strict core' category in output for rows where **all** columns have a '1'

- **-s**, **-singletons**

    Include singleton/column-specific rows for each column header in the output, activates also overall column header presence ('1') counts in final stats matrix for columns with singletons

- **-u**, **-unspecific**

    Include 'unspecific' category in output

- **-a**, **-all\_column\_presence\_overall**

    Report overall presence counts for all column headers (appended to the final stats matrix), also those without singletons; will include all overall column header presence counts without option **-s**

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

### `transpose_matrix.pl` output

- *STDOUT*

    The transposed matrix is printed to *STDOUT*. Redirect or pipe into another tool as needed.

### `binary_group_stats.pl` output

- *STDOUT*

    The tab-delimited final stats matrix is printed to *STDOUT*. Redirect or pipe into another tool as needed.

- ./results_i#_e#

    All output files are stored in a results folder

- ./results_i#_e#/[\*_specific|\*_absent|cutoff_core|underrepresented]_rows.txt

    Files including the row headers for rows in non-optional categories

- (./results_i#_e#/[\*_singletons|strict_core|unspecific]_rows.txt)

    Optional category output files with the respective row headers

- (./results_i#_e#/binary_matrix.tsv)

    Tab-delimited binary matrix of group presence/absence results according to cutoffs (excluding 'unspecific' rows)

- (./results_i#_e#/venn_diagram.pdf)

    Venn diagram for non-optional categories (except 'unspecific' and 'underrepresented' categories)

## Dependencies

### `prot_finder.pl`/`prot_finder_pipe.sh` dependencies

- [**BioPerl**](http://www.bioperl.org) (tested version 1.006923)
- [`cds_extractor`](/cds_extractor) (tested version 0.7.1)
- [**BLASTP**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - legacy **BLASTP** (tested version 2.2.26)
    - **BLASTP+** (tested version 2.2.28+)
- [**Clustal Omega**](http://www.clustal.org/omega/) (tested version 1.2.1)

### `binary_group_stats.pl` dependencies

- **Statistical computing language [R](http://www.r-project.org/)**

    `Rscript` is needed to plot the venn diagram with option **-p**, tested with version 3.2.2

- **gplots (https://cran.r-project.org/web/packages/gplots/index.html)**

    Package needed for **R** to plot the venn diagram, includes function `venn`. Tested with **gplots** version 2.17.0.

## Run environment

The scripts run under UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Acknowledgements

The Perl implementation for transposing a matrix on Stack Overflow
was very useful for `transpose_matrix.pl`: https://stackoverflow.com/questions/1729824/transpose-a-file-in-bash

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

### transpose_matrix changelog

* v0.1 (12.04.2016)

### binary_group_stats changelog

* v0.1 (06.06.2016)
