tbl2tab
=======

`tbl2tab.pl` is a script to convert tbl to tab-separated format and back.

* [Synopsis](#synopsis)
* [Description](#description)
* [Usage](#usage)
* [Options](#options)
  * [Mandatory options](#mandatory-options)
  * [Optional options](#optional-options)
* [Output](#output)
* [Run environment](#run-environment)
* [Author - contact](#author---contact)
* [Citation, installation, and license](#citation-installation-and-license)
* [Changelog](#changelog)

## Synopsis

    perl tbl2tab.pl -m tbl2tab -i feature_table.tbl -s -l locus_prefix

**or**

    perl tbl2tab.pl -m tab2tbl -i feature_table.tab -g -l locus_prefix -p "gnl|dbname|"

## Description

NCBI's feature table (**tbl**) format is needed for the submission of genomic data to GenBank with the NCBI tools [Sequin](http://www.ncbi.nlm.nih.gov/Sequin/) or [tbl2asn](http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2). tbl files can be created with automatic annotation systems like [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml). `tbl2tab.pl` can convert a tbl file to a tab-separated format (tab) and back to the tbl format. The tab-delimited format is useful to manipulate the data more comfortably in a spreadsheet software (e.g. LibreOffice or MS Excel). For a conversion back to tbl format save the file in the spreadsheet software as a tab-delimited text file. The script is intended for microbial genomes, but might also be useful for eukaryotes.

Regular expressions are applied in mode '**tbl2tab**' to correct gene names and words in '/product' values  to lowercase initials (with the exception of 'Rossman' and 'Willebrand'). The resulting tab file can then be used to check for possible errors.

The first four header columns of the **tab** format are mandatory, 'seq_id' for the SeqID, and for each primary tag/feature (e.g. CDS, RNAs, repeat_region etc.), 'start', 'stop', and 'primary_tag'. These mandatory columns have to be filled in every row in the tab file. All the following columns will be included as tags/qualifiers (e.g. '/locus_tag', '/product', '/EC_number', '/note' etc.) in the conversion to the tbl file if a value is present.

There are three special cases:

**First**, '/pseudo' will be included as a tag if *any* value (the script uses 'T' for true) is present in the **tab** format. If a primary tag is indicated as pseudo both the primary tag and the accessory 'gene' primary tag (for CDS/RNA features with option **-g**) will include a '/pseudo' qualifier in the resulting **tbl** file. *Pseudo-genes* are indicated by 'pseudo' in the 'primary_tag' column, thus the 'pseudo' column is ignored in these cases.

**Second**, tag '/gene_desc' is reserved for the 'product' values of pseudo-genes, thus a 'gene_desc' column in a tab file will be ignored in the conversion to tbl.

**Third**, column 'protein_id' in a tab file will also be ignored in the conversion. '/protein_id' values are created from option **-p** and the locus_tag for each CDS primary feature.

Furthermore, with option **-s** G2L-style spreadsheet formulas ([Goettingen Genomics Laboratory](http://appmibio.uni-goettingen.de/)) can be included with additional columns, 'spreadsheet_locus_tag', 'position', 'distance', 'gene_number', and 'contig_order'. These columns will not be included in a conversion to the tbl format. Thus, if you want to include e.g. the locus_tags from the formula in column 'spreadsheet_locus_tag' in the resulting tbl file copy the *values* to the column 'locus_tag'!

To illustrate the process two example files are included in the repository, 'example.tbl' and 'example2.tab', which are interconvertible (see "[USAGE](#usage)" below).

**Warning**, be aware of possible errors introduced by automatic format conversions using a spreadsheet software like MS Excel, see e.g. Zeeberg *et al.* 2004 (http://www.ncbi.nlm.nih.gov/pubmed/15214961).

For more information regarding the feature table and the submission process see NCBI's [prokaryotic annotation guide](http://www.ncbi.nlm.nih.gov/genbank/genomesubmit) and the [bacterial genome submission guide](http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation).

## Usage

### Conversion from tbl to tab format

    perl tbl2tab.pl -m tbl2tab -i example.tbl -s -l EPE

### Conversion from tab to tbl format

    perl tbl2tab.pl -m tab2tbl -i example2.tab -g -l EPE

## Options

### Mandatory options

* -m, -mode

Conversion mode, either 'tbl2tab' or 'tab2tbl' [default = 'tbl2tab']

* -i, -input

Input tbl or tab file to be converted to the other format

### Optional options

* -h, -help

Help (perldoc POD)

* -v, -version

Print version number to *STDERR*

#### Mode *tbl2tab*

* -l, -locus_prefix

Only in combination with option **-s** and there mandatory to include the locus_tag prefix in the formula for column 'spreadsheet_locus_tag'

* -c, -concat

Concatenate values of identical tags within one primary tag with '~' (e.g. several '/EC_number' or '/inference' tags)

* -e, -empty

String used for primary features without value for a tag [default = '']

* -s, -spreadsheet

Include formulas for spreadsheet editing

* -f, -formula_lang

Syntax language of the spreadsheet formulas, either 'English' or 'German'. If you're still encountering problems with the formulas set the decimal and thousands separator manually in the options of the spreadsheet software (instead of using the operating system separators). [default = 'e']

#### Mode *tab2tbl*

* -l, -locus_prefix

Prefix to the SeqID if not present already in the SeqID

* -g, -gene

Include accessory 'gene' primary tags (with '/gene', '/locus_tag' and possibly '/pseudo' tags) for 'CDS/RNA' primary tags; NCBI standard

* -t, -tags_full

Only in combination with option **-g**, include '/gene' and '/locus_tag' tags additionally in primary tag, not only in accessory 'gene' primary tag

* -p, -protein_id_prefix

Prefix for '/protein_id' tags; don't forget the double quotes for the string, otherwise the shell will intepret as pipe [default = 'gnl|goetting|']

## Output

* *.tab|tbl

Result file in the opposite format

* (hypo_putative_genes.txt)

Created in mode **tab2tbl**, indicates if CDSs are annotated as
'hypothetical/putative/predicted protein' but still have a gene name

## Run environment

The Perl script runs under Windows and UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.2 (29.10.2014)
    * fixed bug: message which file was created was mixed up
    * *hypo_putative_genes.txt* includes now also 'predicted protein' annotations
    * additions and syntax changes to POD and README.md
* v0.1 (24.06.2014)
