po2group_stats
==============

`po2group_stats.pl` is a script to categorize orthologs from [Proteinortho5](http://www.bioinf.uni-leipzig.de/Software/proteinortho/) output according to genome groups. In the [prot_finder](/prot_finder) workflow is a script, `binary_group_stats.pl`, which does the same thing for column groups in a delimited TEXT binary matrix.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
  * [cds_extractor](#cds_extractor)
  * [Proteinortho5](#proteinortho5)
  * [po2group_stats](#po2group_stats)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Dependencies](#dependencies)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl po2group_stats.pl -i matrix.proteinortho -d genome_fasta_dir/ -g group_file.tsv -p > overall_stats.tsv

## Description

Categorize the genomes in an ortholog/paralog output matrix (option **-i**) from a
[**Proteinortho5**](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)
calculation according to group affiliations. The group
affiliations of the genomes are intended to get overall
presence/absence statistics for groups of genomes and not simply
single genomes (e.g. comparing 'marine', 'earth', 'commensal',
'pathogenic' etc. genome groups). Percentage inclusion (option
**-cut\_i**) and exclusion (option **-cut\_e**) cutoffs can be set to
define how strict the presence/absence of genome groups within an
orthologous group (OG) are defined. Of course groups can also hold
only single genomes to get single genome statistics. Group
affiliations are defined in a mandatory **tab-delimited** group input
file (option **-g**) with **minimal two** and **maximal four** groups.

Only alphanumeric (a-z, A-Z, 0-9), underscore (\_), dash (-), and
period (.) characters are allowed for the **group names** in the
group file to avoid downstream problems with the operating/file
system. As a consequence, also no whitespaces are allowed in these!
Additionally, **group names**, **genome filenames** (should be
enforced by the file system), and **FASTA IDs** considering **all**
genome files (mostly locus tags; should be enforced by Proteinortho5)
need to be **unique**.

**Proteinortho5** (PO) has to be run with option **-singles** to
include also genes without orthologs, so-called singletons/ORFans,
for each genome in the PO matrix (see the
[PO manual](http://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html)).
Additionally, option **-selfblast** is recommended to enhance
paralog detection by PO.

To explain the logic behind the categorization, the following
annotation for example groups will be used. A '1' exemplifies a
group genome count in a respective OG >= the rounded inclusion
cutoff, a '0' a group genome count <= the rounded exclusion cutoff.
The presence and absence of OGs for the group affiliations are
structured in different categories depending on the number of
groups. For **two groups** (e.g. A and B) there are five categories:
'A specific' (A:B = 1:0), 'B specific' (0:1), 'cutoff core' (1:1),
'underrepresented' (0:0), and 'unspecific'. Unspecific OGs have a
genome count for at least **one** group outside the cutoffs
(exclusion cutoff < genome count < inclusion cutoff) and
thus cannot be categorized. These 'unspecific' OGs will only be
printed to a final annotation result file with option **-u**. Overall
stats for all categories are printed to *STDOUT* in a final
tab-delimited output matrix.

**Three groups** (A, B, and C) have the following nine categories: 'A
specific' (A:B:C = 1:0:0), 'B specific' (0:1:0), 'C specific'
(0:0:1), 'A absent' (0:1:1), 'B absent' (1:0:1), 'C absent' (1:1:0),
'cutoff core' (1:1:1), 'underrepresented' (0:0:0), and 'unspecific'.

**Four groups** (A, B, C, and D) are classified in 17 categories: 'A
specific' (A:B:C:D = 1:0:0:0), 'B specific' (0:1:0:0), 'C specific'
(0:0:1:0), 'D specific' (0:0:0:1), 'A-B specific' (1:1:0:0), 'A-C
specific' (1:0:1:0), 'A-D specific' (1:0:0:1), 'B-C specific'
(0:1:1:0), 'B-D specific' (0:1:0:1), 'C-D specific' (0:0:1:1), 'A
absent' (0:1:1:1), 'B absent' (1:0:1:1), 'C absent' (1:1:0:1), 'D
absent' (1:1:1:0), 'cutoff core' (1:1:1:1), 'underrepresented'
(0:0:0:0), and 'unspecific'.

The resulting group presence/absence (according to the cutoffs) can
also be printed to a binary matrix (option **-b**) in the result
directory (option **-r**), excluding the 'unspecific' category. Since
the categories are the logics underlying venn diagrams, you can also
plot the results in a venn diagram using the binary matrix (option
**-p**). The 'underrepresented' category is exempt from the venn
diagram, because it is outside of venn diagram logics.

Here are venn diagrams illustrating the logic categories (see folder ['pics'](./pics)):

<p align="center">
  <img alt="venn diagram logics" src="https://github.com/aleimba/bac-genomics-scripts/raw/master/po2group_stats/pics/venn_diagram_logics.png">
</p>

There are two optional categories (which are only considered for the
final print outs and in the final stats matrix, not for the binary
matrix and the venn diagram): 'strict core' (option **-co**) for
OGs where **all** genomes have an ortholog, independent of the
cutoffs. Of course all the 'strict core' OGs are also included in
the 'cutoff\_core' category ('strict core' is identical to 'cutoff
core' with **-cut\_i** 1 and **-cut\_e** 0). Option **-s** activates the
detection of 'singleton/ORFan' OGs present in only **one** genome.
Depending on the cutoffs and number of genomes in the groups,
category 'underrepresented' includes most of these singletons.

Additionally, annotation is retrieved from multi-FASTA files created
with [`cds_extractor.pl`](/cds_extractor). See
[`cds_extractor.pl`](/cds_extractor) for a description of the
format. These files are used as input for the PO analysis and with
option **-d** for `po2group_stats.pl`. The annotations are printed
in category output files in the result directory.

Annotations are only pulled from one representative genome for each
category present in the current OG. With option **-co** you can set a
specific genome for the representative annotation for category
'strict core'. For all other categories the representative genome is
chosen according to the order of the genomes in the group files,
depending on the presence in each OG. Thus, the best annotated
genome should be in the first group at the topmost position
(especially for 'cutoff core'), as well as the best annotated ones
at the top in all other groups.

In the result files, each orthologous group (OG) is listed in a row
of the resulting category files, the first column holds the OG
numbers from the PO input matrix (i.e. line number minus one). The
following columns specify the ID for each CDS, gene, EC number(s),
product, and organism are shown depending on their presence in the
CDS's annotation. The ID is in most cases the locus tag (see
[`cds_extractor.pl`](/cds_extractor)). If several EC numbers exist
for a single CDS they are separated by a ';'. If the representative
genome within an OG includes paralogs (co-orthologs) these will be
printed in the following row(s) **without** a new OG number in the
first column.

The number of OGs in the category annotation result files are the
same as listed in the venn diagram and the final stats matrix.
However, since only annotation from one representative annotation is
used the CDS number will be different to the final stats. The final
stats include **all** the CDS in this category in **all** genomes
present in the OG in groups >= the inclusion cutoff (i.e. for
'strict core' the CDS for all genomes in this OG are counted). Two
categories are different, for 'unspecific' all unspecific groups are
included, for 'underrepresented' all groups <= the exclusion
cutoffs. This is also the reason, the 'pangenome' CDS count is
greater than the 'included in categories' CDS count in the final
stats matrix, as genomes in excluded groups are exempt from the CDS
counts for most categories. 'Included in categories' is the OG/CDS
sum of all non-optional categories ('\*specific', '\*absent', 'cutoff
core', 'underrepresented', and 'unspecific'), since the optional
categories are included in non-optionals. An exception to the
difference in CDS counts are the 'singletons' category where OG and
CDS counts are identical in the result files and in the overall
final output matrix (as there is only one genome), as well as in
group-'specific' categories for groups including only one genome.

At last, if you want the respective representative sequences for a
category you can first filter the locus tags from the result file
with Unix command-line tools:

    grep -v "^#" result_file.tsv | cut -f 2 > locus_tags.txt

And then feed the locus tag list to
[`cds_extractor.pl`](/cds_extractor) with option **-l**.

As a final note, in the [prot_finder](/prot_finder) workflow is a
script, `binary_group_stats.pl`, based upon `po2group_stats.pl`,
which can calculate overall presence/absence statistics for column
groups in a delimited TEXT binary matrix (as with genomes here).

## Usage

### [`cds_extractor`](/cds_extractor)

    for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i [-p|-n]; done

### [**Proteinortho5**](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)

    proteinortho5.pl -graph [-synteny] -cpus=# -selfblast -singles -identity=50 -cov=50 -blastParameters='-use_sw_tback [-seg no|-dust no]' *.[faa|ffn]

### po2group_stats

    perl po2group_stats.pl -i matrix.[proteinortho|poff] -d genome_fasta_dir/ -g group_file.tsv -r result_dir -cut_i 0.7 -cut_e 0.2 -b -p -co genome4.[faa|ffn] -s -u -a > overall_stats.tsv

## Options

### Mandatory options

- **-i**=_str_, **-input**=_str_

    Proteinortho (PO) result matrix (\*.proteinortho or \*.poff)

- **-d**=_str_, **-dir\_genome**=_str_

    Path to the directory including the genome multi-FASTA PO input files (\*.faa or \*.ffn), created with [`cds_extractor.pl`](/cds_extractor)

- **-g**=_str_, **-groups\_file**=_str_

    Tab-delimited file with group affiliation for the genomes with **minimal two** and **maximal four** groups (easiest to create in a spreadsheet software and save in tab-separated format). **All** genomes from the PO matrix need to be included. Group names can only include alphanumeric (a-z, A-Z, 0-9), underscore (\_), dash (-), and period (.) characters (no whitespaces allowed either). Example format with two genomes in group A, three genomes in group B and D, and one genome in group C:

    group\_A&emsp;group\_B&emsp;group\_C&emsp;group\_D<br>
    genome1.faa&emsp;genome2.faa&emsp;genome3.faa&emsp;genome4.faa<br>
    genome5.faa&emsp;genome6.faa&emsp;&emsp;genome7.faa<br>
    &emsp;genome8.faa&emsp;&emsp;genome9.faa

### Optional options

- **-h**, **-help**

    Help (perldoc POD)

- **-r**=_str_, **-result\_dir**=_str_

    Path to result folder \[default = inclusion and exclusion percentage cutoff, './results\_i#\_e#'\]

- **-cut\_i**=_float_, **-cut\_inclusion**=_float_

    Percentage inclusion cutoff for genomes in a group per OG, has to be > 0 and <= 1. Cutoff will be rounded according to the genome number in each group and has to be > the rounded exclusion cutoff in this group. \[default = 0.9\]

- **-cut\_e**=_float_, **-cut\_exclusion**=_float_

    Percentage exclusion cutoff, has to be >= 0 and < 1. Rounded cutoff has to be < rounded inclusion cutoff. \[default = 0.1\]

- **-b**, **-binary\_matrix**

    Print a binary matrix with the presence/absence genome group results according to the cutoffs (excluding 'unspecific' category OGs)

- **-p**, **-plot\_venn**

    Plot venn diagram from the binary matrix (except 'unspecific' and 'underrepresented' categories) with function `venn` from **R** package **gplots**, requires option **-b**

- **-co**=(_str_), **-core_strict**=(_str_)

    Include 'strict core' category in output. Optionally, give a genome name from the PO matrix to use for the representative output annotation. \[default = topmost genome in first group\]

- **-s**, **-singletons**

    Include singletons/ORFans for each genome in the output, activates also overall genome OG/CDS stats in final stats matrix for genomes with singletons

- **-u**, **-unspecific**

    Include 'unspecific' category representative annotation file in result directory

- **-a**, **-all\_genomes\_overall**

    Report overall stats for all genomes (appended to the final stats matrix), also those without singletons; will include all overall genome stats without option **-s**

- **-v**, **-version**

    Print version number to *STDERR*

## Output

- *STDOUT*

    The tab-delimited final stats matrix is printed to *STDOUT*. Redirect or pipe into another tool as needed.

- ./results_i#_e#

    All output files are stored in a results folder

- ./results_i#_e#/[\*_specific|\*_absent|cutoff_core|underrepresented]_OGs.tsv

    Tab-delimited files with OG annotation from a representative genome for non-optional categories

- (./results_i#_e#/[\*_singletons|strict_core|unspecific]_OGs.tsv)

    Optional category tab-delimited output files with representative annotation

- (./results_i#_e#/binary_matrix.tsv)

    Tab-delimited binary matrix of group presence/absence results according to cutoffs (excluding 'unspecific' category)

- (./results_i#_e#/venn_diagram.pdf)

    Venn diagram for non-optional categories (except 'unspecific' and 'underrepresented' categories)

## Dependencies

- **Statistical computing language [R](http://www.r-project.org/)**

    `Rscript` is needed to plot the venn diagram with option **-p**, tested with version 3.2.2

- **gplots (https://cran.r-project.org/web/packages/gplots/index.html)**

    Package needed for **R** to plot the venn diagram with function `venn`. Tested with **gplots** version 2.17.0.

## Run environment

The Perl script runs under UNIX and Windows flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.1.3 (06.06.2016)
    * included check for file system conformity for group names
    * some minor syntax changes and additions to error messages, basically adapting to [`binary_group_stats.pl`](/prot_finder)
* v0.1.2 (19.11.2015)
    * added `pod2usage`-die for Getopts::Long call
    * minor POD/README change
* v0.1.1 (30.10.2015)
    * fixed bug for representative annotation in output files, the representative genome was not chosen according to genome order in the groups file
* v0.1 (23.10.2015)
