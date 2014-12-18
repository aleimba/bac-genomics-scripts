po2anno
=======

`po2anno.pl` is a script to create an annotation comparison matrix from Proteinortho5 output.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [cds_extractor](#cds_extractor)
  * [Proteinortho5](#proteinortho5)
  * [po2anno](#po2anno)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Changelog](#changelog)

## Synopsis

    perl po2anno.pl -i matrix.proteinortho -g genome_fasta_dir/ -l -a > annotation_comparison.tab

## Description

Supplement an ortholog/paralog output matrix from a
[*Proteinortho5*](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)
calculation with annotation information. The resulting tab-separated
annotation comparison matrix (ACM) is mainly intended for the
transfer of high quality annotations from reference genomes to
homologs (orthologs and co-orthologs/paralogs) in a query genome
(e.g. in conjunction with [`tbl2tab.pl`](/tbl2tab)). But of course
it can also be used to have a quick glance at the annotation of
genes present only in a couple of input genomes in comparison to the
others.

Annotation is retrieved from multi-FASTA files created with
[`cds_extractor.pl`](/cds_extractor). See
[`cds_extractor.pl`](/cds_extractor) for a description of the
format. These files are used as input for the PO analysis.

*Proteinortho5* (PO) has to be run with option **-singles** to include
also genes without orthologs, so-called singletons/ORFans, for each
genome in the PO matrix (see the
[PO manual](http://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html)).
Additionally, option **-selfblast** is recommended to enhance paralog
detection by PO.

Each orthologous group (OG) is listed in a row of the resulting ACM,
the first column holds the OG numbers from the PO input matrix (i.e.
line number minus one). The following columns specify the
orthologous CDS for each input genome. For each CDS the ID,
optionally the length in bp (option **-l**), gene, EC number(s), and
product are shown depending on their presence in the CDS's
annotation. The ID is in most cases the locus tag (see
[`cds_extractor.pl`](/cds_extractor)). If several EC numbers exist
for a single CDS they're separated by ';'. If an OG includes
paralogs, i.e. co-orthologs from a single genome, these will be
printed in the following row(s) **without** a new OG number in the
first column. The order of paralogous CDSs within an OG is
arbitrarily.

The OGs are sorted numerically via the query ID (see option **-q**).
If option **-a** is set, the non-query OGs are appended to the output
after the query OGs, sorted numerically via OG number.

## Usage

### cds_extractor

    for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i [-p|-n]; done
    rename 's/_cds_[aa|nuc].fasta/.[faa|fna]/' *_cds_[aa|nuc].fasta

### Proteinortho5

    proteinortho5.pl -graph [-synteny] -cpus=# -selfblast -singles -identity=50 -cov=50 -blastParameters='-use_sw_tback' *.[faa|fna]

### po2anno

    perl po2anno.pl -i matrix.[proteinortho|poff] -g genome_fasta_dir/ -q query.[faa|fna] -l -a > annotation_comparison.tab

## Options

### Mandatory options

- -i, -input

    Proteinortho (PO) result matrix (\*.proteinortho or \*.poff), or piped *STDIN* (-)

- -g, -genome_dir

    Path to the directory including the genome multi-FASTA PO input files, created with `cds_extractor.pl`

### Optional options

- -h, -help

    Help (perldoc POD)

- -q, -query

    Query genome (has to be identical to the string in the PO matrix) [default = first one in alphabetical order]

- -l, -length

    Include length of each CDS in bp

- -a, -all

    Append non-query orthologous groups (OGs) to the output

- -v, -version

    Print version number to *STDERR*

## Output

- *STDOUT*

    The resulting ACM is printed to *STDOUT*. Redirect or pipe into another tool as needed (e.g. cut, grep, head, or tail).

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Changelog

* v0.1 (18.12.2014)
