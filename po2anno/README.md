po2anno
=======

`po2anno.pl` is a script to create an annotation comparison matrix from [Proteinortho5](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) output.

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
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl po2anno.pl -i matrix.proteinortho -d genome_fasta_dir/ -l -a > annotation_comparison.tsv

## Description

Supplement an ortholog/paralog output matrix from a
[**Proteinortho5**](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)
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
format. These files are used as input for the PO analysis and option
**-d** for `po2anno.pl`.

**Proteinortho5** (PO) has to be run with option **-singles** to include
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

### [`cds_extractor`](/cds_extractor)

    for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i [-p|-n]; done

### [**Proteinortho5**](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)

    proteinortho5.pl -graph [-synteny] -cpus=# -selfblast -singles -identity=50 -cov=50 -blastParameters='-use_sw_tback [-seg no|-dust no]' *.[faa|ffn]

### po2anno

    perl po2anno.pl -i matrix.[proteinortho|poff] -d genome_fasta_dir/ -q query.[faa|ffn] -l -a > annotation_comparison.tsv

## Options

### Mandatory options

- **-i**=_str_, **-input**=_str_

    Proteinortho (PO) result matrix (\*.proteinortho or \*.poff), or piped *STDIN* (-)

- **-d**=_str_, **-dir\_genome**=_str_

    Path to the directory including the genome multi-FASTA PO input files (\*.faa or \*.ffn), created with [`cds_extractor.pl`](/cds_extractor)

### Optional options

- **-h**, **-help**

    Help (perldoc POD)

- **-q**=_str_, **-query**=_str_

    Query genome (has to be identical to the string in the PO matrix) [default = first one in alphabetical order]

- **-l**, **-length**

    Include length of each CDS in bp

- **-a**, **-all**

    Append non-query orthologous groups (OGs) to the output

- **-v**, **-version**

    Print version number to *STDERR*

## Output

- *STDOUT*

    The resulting tab-delimited ACM is printed to *STDOUT*. Redirect or pipe into another tool as needed (e.g. `cut`, `grep`, `head`, or `tail`).

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2.2 (23.10.2015)
    * minor syntax changes to `po2anno.pl` and README
    * changed option **-g|-genome_dir** to **-d|-dir_genome** for consistency with [`po2group_stats.pl`](/po2group_stats)
* v0.2.1 (07.09.2015)
    * get rid of underscores in product annotation strings (from [`cds_extractor.pl`](/cds_extractor))
    * debugged hard-coded relative path for `$genome_file_path`
* v0.2 (15.01.2015)
    * give number of query-specific OGs and total query singletons/ORFans in final stat output
    * changed final stat output to an easier readable format
    * fixed bug: %Query_ID_Seen included also non-query IDs, which luckily had no consequences
* v0.1 (18.12.2014)
