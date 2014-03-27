bac-genomics-scripts
====================

A collection of scripts intended for bacterial genomics (some might be useful also for eukaryotes ;-)).

Here is a short summary of the tools (if not here yet will be uploaded soon):

* Protein search with *blastp* plus concise hit summary and optional alignment: `prot_finder`
* Regions of difference (rod) detection: `rod_finder`
* Extraction of protein/nucleotide sequences from CDSs: `cds_extractor`
* COG classification of proteins: `cdd2cog`
* MLST assignment for *E. coli* (Achtman scheme): `ecoli_mlst`
* Batch downloading of sequences from NCBI's FTP server: `ncbi_ftp_download` and `ncbi_e-utilities`
* NGS paired-end library insert size estimation from BAM/SAM: `sam_insert-size`
* Count all annotated primary features from RichSeq (e.g. embl or genbank) files: `get_genome_features`
* And an assortment of smaller scripts for tasks like: concatenation of RichSeq sequence files, sequence/alignment format converters, dnadiff, GC% calculation etc.

Anyways you get the picture ...

## Introduction

All the scripts here are written in **Perl** (some include bash shell wrappers).

Each script is hosted in its own folder, so that a separate 'README.md' can be included for more information. However, all of the Perl scripts include additionally a usage/help text or a comprehensive POD by calling the script either without arguments/options or option **-h**.

The scripts are only tested under UNIX, some won't run in a Windows environment (because of included UNIX commands). If you are on Windows an alternative might be **Cygwin**: http://cygwin.com/

## Installation recommendations

To download the repository, use either the 'Download zip' button on the right hand side or clone the repository with `git`:

    git clone https://github.com/aleimba/bac-genomics-scripts.git

To install the scripts, copy them to your '*/bin*' folder and make them executable

    $ chmod u+x script.pl

the scripts are then in your '*$PATH*' and can be run everywhere on your system. Of course you can just call them directly by prefexing '*perl*' to the command or a './' for bash wrappers:

    $ perl /path/to/script/script.pl <options>

or

    $ ./script.sh <arguments>

## Dependencies

Most of the Perl scripts include modules from **BioPerl**, which as a consequence has to be installed on your system. For BioPerl installation instructions see: http://www.bioperl.org

Some scripts need additional Perl modules, which will be stated in the associated 'README.md' or POD. Install them from CPAN (installation instructions can be found on the website as well, see FAQ): http://www.cpan.org/

## UNIX loops

A very handy tip, if you want to run a script on all files in the current working directory you can use a **loop** in UNIX, e.g.:

    $ for i in *.fasta; do perl script.pl -i $i; done

## Windows UNIX linebreak problems

At last, some of the scripts don't like Windows formatted line breaks, you might consider running these input files through a nifty UNIX utility called **dos2unix**:

    $ dos2unix input

## License

All scripts are licensed under GPLv3 which is contained in the file *LICENSE*.

For help, suggestions, bugs etc. use the GitHub issues or write an email to aleimba [at] gmx [dot] de.
