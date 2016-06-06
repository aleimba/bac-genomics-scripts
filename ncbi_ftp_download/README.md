ncbi_ftp_download
=================

**This pipeline is NOT working at the moment, as NCBI reorganized the structure of their [FTP server for genomes](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/). As an alternative way to fetch bacterial genomes from NCBI I recommend [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download) from @kbiln, or [`Bio-RetrieveAssemblies`](https://github.com/andrewjpage/Bio-RetrieveAssemblies) from @andrewjpage from the Wellcome Trust Sanger Institute.**

Scripts to batch download all bacterial genomes of a genus/species from NCBI's FTP site (RefSeq and GenBank) for easy access.

## Synopsis

    ncbi_ftp_download.sh Genus_species

## Description

These scripts are intended to download all bacterial genomes for a particular genus or species from NCBI's FTP site (http://www.ncbi.nlm.nih.gov/Ftp/ and ftp://ftp.ncbi.nlm.nih.gov/) and copy them to result folders for easy access.

`ncbi_ftp_download.sh` is a bash shell wrapper script that employs UNIX's `wget` to download microbial genomes in genbank (\*.gbk) and fasta (\*.fna) format from the GenBank and RefSeq databases (NCBI Reference Sequence Database, http://www.ncbi.nlm.nih.gov/refseq/) on NCBI's FTP server, which can be accessed anonymously. As first argument it takes the bacterial genus or species name you want to download (it uses that name with a glob inside the script, e.g. Escherichia_coli will be used as Escherichia_coli\*), see examples below in [usage](#usage). Have a look on the NCBI FTP server to get the correct name (either with your browser or e.g. with FileZilla, http://filezilla-project.org/). If you want to download genomes for several distinct species just run the script with different arguments repeatedly.

The `wget` parameters are specified to keep the FTP server folder structure and mirror it locally downstream from the current working directory (folder 'ftp.ncbi.nlm.nih.gov' will be the top folder of the new folder structure). If you update an already existing folder structure, `wget` will only download and replace files if they are in a newer version on NCBI's FTP server. **But** be aware that NCBI shuffles files around (including new ones, deleting old ones etc.), thus it might be useful to remove 'ftp.ncbi.nlm.nih.gov' and download everything new.

After the download with `wget`, `ncbi_ftp_download.sh` will run the Perl script `ncbi_ftp_concat_unpack.pl`. This script unpacks (draft genomes are stored as tarballs, \*.tgz) and concatenates all complete and draft genomes, which are present in the folder 'ftp.ncbi.nlm.nih.gov' in the current working directory. The script traverses the downloaded NCBI ftp-folder structure and thus has to be called from the top level (containing the folder 'ftp.ncbi.nlm.nih.gov'). `ncbi_ftp_download.sh` runs `ncbi_ftp_concat_unpack.pl` with both **genbank** and **refseq** options, as well as option **y** to overwrite the old result folders (see below [options](#options)). Both scripts have to be in the same directory (or in the path) to run `ncbi_ftp_download.sh`.

For **complete** genomes **plasmids** are concatenated to the **chromosomes** to create multi-genbank/-fasta files (script `split_multi-seq_file.pl` can be used to split the multi-sequence file to single-sequence files).

In **draft** genomes, **scaffold** and/or **contig** files, designated by 'draft_scaf' or 'draft_con', are controlled for annotation (i.e. if gene primary feature tags exist); usually only one of those contains annotations. The one with annotation is then used to create multi-genbank files. Multi-fasta files are created for the corresponding genbank file or, if no annotation exists, for the file which contains more sequence information (either contigs or scaffolds). In the case, that the sequence information is equal, scaffold files are preferred. If sequence size discrepancies between a genbank and its corresponding fasta file are found, error file 'seq_errors.txt' will be created and indicate the villains.

As a suggestion, pick the genomes you're looking for **first** out of './refseq' and the rest out of './genbank'. RefSeq genomes have a higher annotation quality, while GenBank includes more genomes.

Depending on the amount of data to download, the whole process can take quite a while. Also have a mind for space requirements, e.g. all *E. coli*/*Shigella* genomes (March 2014) have a final total space requirement of ~58 GB ('ftp.ncbi.nlm.nih.gov' = ~18 GB; ./genbank = ~25 GB; ./refseq = ~16 GB)!

If you're new to the NCBI FTP site you should read an excellent overview for microbial RefSeq genomes on NCBI's FTP site on Torsten Seemann's blog: http://thegenomefactory.blogspot.de/2012/07/navigating-microbial-genomes-on-ncbi.html.

You can also access an introductory talk for the microbial NCBI FTP resources at figshare (http://figshare.com/articles/Introduction_to_NCBI_s_FTP_server_for_bacterial_genomes/972893). It might be a good idea to read the blog post and have a look in the PDF to have a general idea what's going on, but of course you can just run the scripts and work with the genome files.

## Usage

### 1.) Manual consecutively

#### 1.1.) `wget`

Download RefSeq complete genomes (in fasta and genbank format):

    wget -cNrv -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Genus_species*" -P .

Download RefSeq draft genomes as tarballs:

    wget -cNrv -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT/Genus_species*" -P .

The same procedure has to be followed for GenBank files, here complete genomes:

    wget -cNrv -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Genus_species*" -P .

And finally download GenBank draft genomes:

    wget -cNrv -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/Genus_species*" -P .

#### 1.2.) `ncbi_ftp_concat_unpack.pl`

    perl ncbi_ftp_concat_unpack.pl refseq y
    perl ncbi_ftp_concat_unpack.pl genbank y

### 2.) With one command: `ncbi_ftp_download.sh` wrapper script

Some examples how you can use the shell script, e.g. download all *E. coli* genomes from NCBI's ftp server:

    ncbi_ftp_download.sh Escherichia_coli

Download all *B. cereus* genomes:

    ncbi_ftp_download.sh Bacillus_cereus

Download all *Paenibacillus* genomes:

    ncbi_ftp_download.sh Paenibacillus

## Options

### *ncbi_ftp_concat_unpack.pl*

* genbank (as first argument)

Copy GenBank genomes (from './ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria\*') as (multi-)sequence files in the result folder './genbank'.

* refseq (as first argument)

Copy RefSeq genomes (from './ftp.ncbi.nlm.nih.gov/genomes/Bacteria\*') as (multi-)sequence files in the result folder './refseq'.

* y (as second argument)

Will delete previous result folders and create new ones (otherwise, the script will ask user if to proceed with overwriting)

## Output

### `ncbi_ftp_download.sh`

* './ftp.ncbi.nlm.nih.gov/'

Mirrors NCBI's FTP server structure and downloads the wanted bacterial genome files in this folder with subfolders

### `ncbi_ftp_concat_unpack.pl`

* './genbank'

Result folder for all **GenBank** genomes

* './refseq'

Result folder for all **RefSeq** genomes

* (seq_errors.txt)

Lists \*.gbk and corresponding \*.fasta files with sequence size discrepancies.

## Run environment

Both the Perl script and the bash-shell script run only under UNIX flavors.

## Dependencies (not in the core Perl modules)

* no extra dependencies

## Authors/contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

### *ncbi_ftp_concat_unpack.pl*

* v0.2.1 (13.07.2015)
    - Adapted all scripts to the new NCBI FTP server address: 'ftp://ftp.ncbi.nlm.nih.gov/'
* v0.2 (21.02.2013)
    - 'seq_errors.txt' error file if sequence size discrepancies between genbank and corresponding fasta file found
    - die with error if 'genbank|refseq' not given as first argument
    - print status message which genome is being processed and what file is kept for draft genomes (e.g. scaffold or contig etc.)
    - bug fixes to test for file existence before running code
    - changed usage to HERE document
* v0.1 (15.09.2012)
