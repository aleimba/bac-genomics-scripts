bac-genomics-scripts
====================

A collection of scripts intended for bacterial genomics (some might also be useful for eukaryotes).

* [Summary](#summary)
* [Introduction](#introduction)
* [Installation recommendations](#installation-recommendations)
* [Dependencies](#dependencies)
* [UNIX loops](#unix-loops)
* [Windows - UNIX linebreak problems](#windows---unix-linebreak-problems)
* [License](#license)

## Summary

* Basic stats for bases and reads in FASTQ files: [`calc_fastq-stats`](/calc_fastq-stats)
* Concatenate multi-sequence files (RichSeq EMBL or GENBANK format, or FASTA format) to a single artificial file: [`cat_seq`](/cat_seq)
* COG ([cluster of orthologous groups](http://www.ncbi.nlm.nih.gov/COG/)) classification of proteins: [`cdd2cog`](/cdd2cog)
* Extraction of protein/nucleotide sequences from CDSs: [`cds_extractor`](/cds_extractor)
* MLST (multilocus sequence typing) assignment and allele extraction for *E. coli* ([Achtman scheme](http://mlst.warwick.ac.uk/mlst/)): [`ecoli_mlst`](/ecoli_mlst)
* Create a feature table for all annotated primary features in RichSeq (EMBL or GENBANK format) files: [`genomes_feature_table`](/genomes_feature_table)
* Batch downloading of sequences from NCBI's FTP server: [`ncbi_ftp_download`](/ncbi_ftp_download) and `ncbi_e-utilities`
* Order sequence entries in FASTA/FASTQ files according to an ID list: [`order_fastx`](/order_fastx)
* Create an ortholog/paralog annotation comparison matrix from [*Proteinortho5*](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) output: [`po2anno`](/po2anno)
* Protein search with *blastp* plus concise hit summary and optional alignment: [`prot_finder`](/prot_finder)
* Rename FASTA ID lines and optionally numerate them: [`rename_fasta_id`](/rename_fasta_id)
* Regions of difference (ROD) detection in genomes with *blastn*: [`rod_finder`](/rod_finder)
* NGS paired-end library insert size estimation from BAM/SAM: [`sam_insert-size`](/sam_insert-size)
* Randomly subsample FASTA, FASTQ, or TEXT files with *reservoir sampling*: [`sample_fastx-txt`](/sample_fastx-txt)
* Convert a sequence file to another format with [BioPerl](http://www.bioperl.org): [`seq_format-converter`](/seq_format-converter)
* Manual curation of annotation in NCBI's TBL format (e.g. from [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) automatic annotation) in a spreadsheet software: [`tbl2tab`](/tbl2tab)
* And an assortment of smaller scripts for tasks like (not yet uploaded to GitHub): alignment format converters, dnadiff, GC% calculation etc.

## Introduction

All the scripts here are written in [**Perl**](https://www.perl.org/) (some include bash shell wrappers).

Each script is hosted in its own folder, so that a separate 'README.md' can be included for more information. However, all of the Perl scripts include additionally a usage/help text or a comprehensive POD (Plain Old Documentation) by calling the script either without arguments/options or option **-h**.

The scripts are only tested under UNIX, some won't run in a Windows environment (because of included UNIX commands). If you are on Windows an alternative might be [Cygwin](http://cygwin.com/).

## Installation recommendations

To download the repository, use either the 'Download zip' button on the right hand side or clone the repository with `git`:

    git clone https://github.com/aleimba/bac-genomics-scripts.git

If there is an update to this GitHub repository (see above [commits](https://github.com/aleimba/bac-genomics-scripts/commits/master) and [releases](https://github.com/aleimba/bac-genomics-scripts/releases)), you can refresh your **local** repository by using the following command **inside** the local folder:

    git pull

To install the scripts, copy them to your '*/bin*' folder and make them executable

    $ chmod u+x script.pl

the scripts are then in your '*$PATH*' and can be run everywhere on your system. Of course you can just call them directly by prefexing '*perl*' to the command or a './' for bash wrappers:

    $ perl /path/to/script/script.pl <options>

or

    $ ./script.sh <arguments>

**Single** scripts can be downloaded as well. For this purpose click on the folder you're interested in and then on the link of the script. There click on the **Raw** button and save this page to a file (without **Raw** you'll get an unusable html file). This is also true for other files (e.g. PDFs etc.).

## Dependencies

All scripts are tested with Perl v5.18.2.

Most of the Perl scripts include modules from [BioPerl](http://www.bioperl.org) as stated in their respective 'README.md', which as a consequence has to be installed on your system. For BioPerl installation instructions see the website (**How Do I...?...install BioPerl?**).

Some scripts need additional Perl modules, which will be stated in the associated 'README.md' or POD. If they're not installed yet on your system get them from [CPAN](http://www.cpan.org/) (installation instructions can be found on the website, see e.g. **Getting Started** or **FAQ**).

## UNIX loops

A very handy tip, if you want to run a script on all files in the current working directory you can use a **loop** in UNIX, e.g.:

    $ for i in *.fasta; do perl script.pl -i $i; done

## Windows - UNIX linebreak problems

At last, some of the scripts don't like Windows formatted line breaks, you might consider running these input files through a nifty UNIX utility called [dos2unix](http://dos2unix.sourceforge.net/):

    $ dos2unix input

## License

All scripts are licensed under GPLv3 which is contained in the file *LICENSE*.

For help, suggestions, bugs etc. use the GitHub issues or write an email to aleimba [at] gmx [dot] de.
