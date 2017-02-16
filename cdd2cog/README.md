cdd2cog
=======

`cdd2cog.pl` is a script to assign COG categories to query protein sequences.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [RPS-BLAST+](#rps-blast)
  * [cdd2cog](#cdd2cog)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Acknowledgements](#acknowledgements)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl cdd2cog.pl -r rps-blast.out -c cddid.tbl -f fun.txt -w whog

## Description
For troubleshooting and a working example please see issue [#1](https://github.com/aleimba/bac-genomics-scripts/issues/1).

The script assigns COG ([cluster of orthologous
groups](http://www.ncbi.nlm.nih.gov/COG/)) categories to proteins.
For this purpose, the query proteins need to be blasted with
RPS-BLAST+ ([Reverse Position-Specific BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download))
against NCBI's Conserved Domain Database
([CDD](http://www.ncbi.nlm.nih.gov/cdd)). Use
[`cds_extractor.pl`](/cds_extractor) beforehand to extract multi-fasta protein
files from GENBANK or EMBL files.

Both tab-delimited RPS-BLAST+ outformats, **-outfmt 6** and **-outfmt
7**, can be processed by `cdd2cog.pl`. By default, RPS-BLAST+ hits
for each query protein are filtered for the best hit (lowest
e-value). Use option **-a|all\_hits** to assign COGs to all BLAST hits
and e.g. do a downstream filtering in a spreadsheet application.
Results are written to tab-delimited files in the './results'
folder, overall assignment statistics are printed to *STDOUT*.

Several files are needed from NCBI's FTP server to run the RPS-BLAST+ and `cdd2cog.pl`:

1. **CDD** (ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/)

    More information about the files in the CDD FTP archive can be found in the respective 'README' file.

  1. 'cddid.tbl.gz'

    The file needs to be unpacked:

    `gunzip cddid.tbl.gz`

    Contains summary information about the CD models in a tab-delimited format. The columns are: PSSM-Id, CD accession (e.g. COG#), CD short name, CD description, and PSSM (position-specific scoring matrices) length.

  2. './little_endian/Cog_LE.tar.gz'

    Unpack and untar via:

    `tar xvfz Cog_LE.tar.gz`

    Preformatted RPS-BLAST+ database of the CDD COG distribution for Intel CPUs and Unix/Windows architectures.

2. **COG** (ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/)

    Read 'readme' for more information about the respective files in the COG FTP archive.

  1. 'fun.txt'

    One-letter functional classification used in the COG database.

  2. 'whog'

    Name, description, and corresponding functional classification of each COG.

## Usage

### RPS-BLAST+

    rpsblast -query protein.fasta -db Cog -out rps-blast.out -evalue 1e-2 -outfmt 6
    rpsblast -query protein.fasta -db Cog -out rps-blast.out -evalue 1e-2 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'

### cdd2cog

    perl cdd2cog.pl -r rps-blast.out -c cddid.tbl -f fun.txt -w whog -a

## Options

### Mandatory options

- -r, -rps\_report

    Path to RPS-BLAST+ report/output, outfmt 6 or 7

- -c, -cddid

    Path to CDD's 'cddid.tbl' file

- -f, -fun

    Path to COG's 'fun.txt' file

- -w, -whog

    Path to COG's 'whog' file

### Optional options

- -h, -help

    Help (perldoc POD)

- -a, -all\_hits

    Don't filter RPS-BLAST+ output for the best hit, rather assign COGs to all hits

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    Overall assignment statistics

- ./results

    All tab-delimited output files are stored in this result folder

- rps-blast_cog.txt

    COG assignments concatenated to the RPS-BLAST+ results for filtering

- protein-id_cog.txt

    Slimmed down 'rps-blast_cog.txt' only including query id (first BLAST report column), COGs, and functional categories

- cog_stats.txt

    Assignment counts for each used COG

- func_stats.txt

    Assignment counts for single-letter functional categories

## Run environment

The Perl script runs under UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Acknowledgements

I got the idea for using NCBI's CDD PSSMs for COG assignment from JGI's [IMG/ER annotation system](http://img.jgi.doe.gov/), which employes the same technique.

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2 (2017-02-16)
    * Adapted to new NCBI FASTA header format for CDD RPS-BLAST+ output
* v0.1 (2013-08-01)
