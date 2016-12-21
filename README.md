[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.215824.svg)](http://dx.doi.org/10.5281/zenodo.215824)

bac-genomics-scripts
====================

A collection of scripts intended for **bacterial genomics** (some might also be useful for eukaryotes) from **high-throughput sequencing** (aka next-generation sequencing).

* [Summary](#summary)
* [Introduction](#introduction)
* [Installation recommendations](#installation-recommendations)
* [Dependencies](#dependencies)
* [UNIX loops](#unix-loops)
* [Windows - UNIX linebreak problems](#windows---unix-linebreak-problems)
* [Citation](#citation)
* [License](#license)
* [Author - contact](#author---contact)

## Summary

* Basic stats for bases and reads in FASTQ files: [`calc_fastq-stats`](/calc_fastq-stats)
* Concatenate multi-sequence files (RichSeq EMBL or GENBANK format, or FASTA format) to a single artificial file: [`cat_seq`](/cat_seq)
* COG ([cluster of orthologous groups](http://www.ncbi.nlm.nih.gov/COG/)) classification of proteins: [`cdd2cog`](/cdd2cog)
* Extraction of protein/nucleotide sequences from CDSs: [`cds_extractor`](/cds_extractor)
* MLST (multilocus sequence typing) assignment and allele extraction for *Escherichia coli* ([Achtman scheme](http://mlst.warwick.ac.uk/mlst/)): [`ecoli_mlst`](/ecoli_mlst)
* Create a feature table for all annotated primary features in RichSeq (EMBL or GENBANK format) files: [`genomes_feature_table`](/genomes_feature_table)
* **Deprecated!** Batch downloading of sequences from NCBI's FTP server: [`ncbi_ftp_download`](/ncbi_ftp_download)
* Order sequence entries in FASTA/FASTQ files according to an ID list: [`order_fastx`](/order_fastx)
* Create an ortholog/paralog annotation comparison matrix from [*Proteinortho5*](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) output: [`po2anno`](/po2anno)
* Calculate stats and plot venn diagrams for genome groups according to orthologs/paralogs from [*Proteinortho5*](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) output, i.e. overall presence/absence statistics for groups of genomes and not simply single genomes: [`po2group_stats`](/po2group_stats)
* Strain panel query protein search with **BLASTP** plus concise hit summary, optional alignment, and presence/absence matrix. Also included, scripts to transpose the matrix and calculate overall presence/absence statistics for groups of columns in the matrix: [`prot_finder`](/prot_finder)
* Rename FASTA ID lines and optionally numerate them: [`rename_fasta_id`](/rename_fasta_id)
* Reverse complement (multi-)sequence files (RichSeq EMBL or GENBANK format, or FASTA format): [`revcom_seq`](/revcom_seq)
* Regions of difference (ROD) detection in genomes with **BLASTN**: [`rod_finder`](/rod_finder)
* NGS paired-end library insert size estimation from BAM/SAM: [`sam_insert-size`](/sam_insert-size)
* Randomly subsample FASTA, FASTQ, or TEXT files with [*reservoir sampling*](https://en.wikipedia.org/wiki/Reservoir_sampling): [`sample_fastx-txt`](/sample_fastx-txt)
* Convert a sequence file to another format with [BioPerl](http://www.bioperl.org): [`seq_format-converter`](/seq_format-converter)
* Manual curation of annotation in NCBI's TBL format (e.g. from [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) automatic annotation) in a spreadsheet software: [`tbl2tab`](/tbl2tab)
* Truncate sequence files (RichSeq EMBL or GENBANK format, or FASTA format) according to given coordinates: [`trunc_seq`](/trunc_seq)
* And an assortment of smaller scripts for tasks like (not yet uploaded to GitHub): alignment format converters, dnadiff, GC% calculation etc.

## Introduction

All the scripts here are written in [**Perl**](https://www.perl.org/) (some include bash shell wrappers).

Each script is hosted in its own folder, so that a separate *README.md* can be included for more information. However, all of the Perl scripts include additionally a usage/help text or a comprehensive [POD](http://perldoc.perl.org/perlpod.html) (Plain Old Documentation) by calling the script either without arguments/options or option **-h|-help**.

The scripts are only tested under UNIX, some won't run in a Windows environment (because of included UNIX commands). If you are on Windows an alternative might be [Cygwin](http://cygwin.com/).

## Installation recommendations

To download the repository, use either the '[Download ZIP](https://github.com/aleimba/bac-genomics-scripts/archive/master.zip)' link after clicking the green 'Clone or download' button at the top or clone the repository with `git`:

    git clone https://github.com/aleimba/bac-genomics-scripts.git

If there is an update to this GitHub repository (see above [commits](https://github.com/aleimba/bac-genomics-scripts/commits/master) and [releases](https://github.com/aleimba/bac-genomics-scripts/releases)), you can refresh your **local** repository by using the following command **inside** the local folder:

    git pull

To install the scripts, copy them e.g. to a home */bin* folder in your *PATH* and make them executable

    $ find . \( -name '*.pl' -o -name '*.sh' -o -name '*.fas' -o -name '*.txt' \) -exec cp {} ~/bin \;
    $ chmod u+x ~/bin/*.pl

the scripts can then be run everywhere on your system. Of course you can just call them directly by prefexing `perl` to the command or a './' for bash wrappers:

    $ perl /path/to/script/script.pl <options>

or

    $ ./script.sh <arguments>

**Single** scripts can be downloaded as well. For this purpose click on the folder you're interested in and then on the link of the script. There click on the **Raw** button and save this page to a file (without **Raw** you'll get an unusable html file). This is also true for other files (e.g. PDFs etc.).

## Dependencies

All scripts are tested with Perl v5.22.1.

Most of the Perl scripts include modules from [BioPerl](http://www.bioperl.org) as stated in their respective *README.md* or POD, which as a consequence has to be installed on your system. For BioPerl installation instructions see the website ([**Installation**](http://bioperl.org/INSTALL.html)).

Some scripts need additional Perl modules, which will be stated in the associated *README.md* or POD. If they're not installed yet on your system get them from [CPAN](http://www.cpan.org/) (installation instructions can be found on the website, see e.g. [**Getting Started...Installing Perl Modules**](http://www.cpan.org/modules/INSTALL.html) or [**FAQ**](http://www.cpan.org/misc/cpan-faq.html#How_install_Perl_modules)).

Furthermore, some scripts call upon statistical computing language [**R**](http://www.r-project.org/) and dependent packages for plotting purposes (again see the respective *README.md* or POD).

## UNIX loops

A very handy tip, if you want to run a script on all files in the current working directory you can use a **loop** in UNIX, e.g.:

    $ for file in *.fasta; do perl script.pl "$file"; done

## Windows - UNIX linebreak problems

At last, some of the scripts don't like Windows formatted line breaks, you might consider running these input files through a nifty UNIX utility called [dos2unix](http://dos2unix.sourceforge.net/):

    $ dos2unix input

## Citation
For now cite the latest major release (tag: [***bovine_ecoli_mastitis***](https://github.com/aleimba/bac-genomics-scripts/releases)) hosted on [Zenodo](https://zenodo.org/):

**Leimbach A**. 2016. bac-genomics-scripts: Bovine *E. coli* mastitis comparative genomics edition. Zenodo. <http://dx.doi.org/10.5281/zenodo.215824>.

Also, all scripts have a version number (see option **-v**), which might be included in a materials and methods section.

## License

All scripts are licensed under GPLv3 which is contained in the file [*LICENSE*](./LICENSE).

## Author - contact
For help, suggestions, bugs etc. use the GitHub issues or write an email to aleimba [at] gmx [dot] de.

Andreas Leimbach (Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)
