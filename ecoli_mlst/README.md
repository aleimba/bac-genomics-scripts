ecoli_mlst
==========

`ecoli_mlst` is a script to determine MLST sequence types for *E. coli* genomes and extract allele sequences.

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

# Synopsis

    perl ecoli_mlst.pl -a fas -g fasta

# Description

The script searches for multilocus sequence type (MLST) alleles in *E. coli* genomes according to
Mark Achtman's scheme with seven house-keeping genes (*adk*, *fumC*, *gyrB*,
*icd*, *mdh*, *purA*, and *recA*) [Wirth et al., 2006]. *NUCmer* from the
[*MUMmer package*](http://mummer.sourceforge.net/) is used to compare the given allele
sequences to bacterial genomes via nucleotide alignments.

Download the allele files (adk.fas ...) and the sequence type file
('publicSTs.txt') from this website:
    http://mlst.ucc.ie/mlst/dbs/Ecoli

To run `ecoli_mlst.pl` include all *E. coli* genome files (file
extension e.g. 'fasta'), all allele sequence files (file extension
'fas') and 'publicSTs.txt' in the current working directory. The
allele profiles are parsed from the created \*.coord files and written
to a result file, plus additional information from the file
'publicSTs.txt'. Also, the corresponding allele sequences (obtained
from the allele input files) are concatenated for each *E. coli* genome
into a result multi-fasta file. Option **-c** can be used to initiate
an alignment for this multi-fasta file with [*ClustalW*](http://www.clustal.org/clustal2/) (standard
alignment parameters; has to be in the `$PATH` or change variable
`$clustal_call`). The alignment fasta output file can be used
directly for [*RAxML*](http://sco.h-its.org/exelixis/web/software/raxml/index.html). CAREFUL the Phylip alignment format from
*ClustalW* allows only 10 characters per strain ID.

`ecoli_mlst.pl` works with complete and draft genomes. However, several genomes cannot be included in a single input file!

Obviously, only for those genomes whose allele sequences have been
deposited in Achtman's allele database results can be obtained. If an
allele is not found in a genome it is marked by a '?' in the result
profile file and a place holder 'XXX' in the result fasta file. For
these cases a manual *NUCmer* or *BLASTN* might be useful to fill the
gaps and [`run_sub_seq.pl`](/run_sub_seq) to get the corresponding 'new' allele
sequences.

Non-NCBI fasta headers for the genome files have to have a
unique ID directly following the '>' (e.g. 'Sakai', '55989' ...).

# Usage

    perl ecoli_mlst.pl -a fas -g fasta -c

# Options

## Mandatory options

- -a, -alleles

    File extension of the MLST allele fasta files, e.g. 'fas' (<=> **-g**).

- -g, -genomes

    File extension of the *E. coli* genome fasta files, e.g. 'fasta' (<=> **-a**).

## Optional options

- -h, -help

    Help (perldoc POD)

- -c, -clustalw

    Call [*ClustalW*](http://www.clustal.org/clustal2/) for alignment

# Output

- ecoli_mlst_profile.txt

    Tab-separated allele profiles for the *E. coli* genomes, plus additional info from 'publicSTs.txt'

- ecoli_mlst_seq.fasta

    Multi-fasta file of all concatenated allele sequences for each genome

- *.coord

    Text files that contain the coordinates of the *NUCmer* hits for each genome and allele

- (errors.txt)

    Error file, summarizing number of not found alleles or unclear *NUCmer* hits

- (ecoli_mlst_seq_aln.fasta)

    Optional, [*ClustalW*](http://www.clustal.org/clustal2/) alignment in Phylip format

- (ecoli_mlst_seq_aln.dnd)

    Optional, *ClustalW* alignment guide tree

## Run environment

The Perl script runs only under UNIX flavors.

## Author - contact

Andreas Leimbach (aleimba[at]gmx[dot]de; Microbial Genome Plasticity, Institute of Hygiene, University of Muenster)

## Citation, installation, and license

For [citation](https://github.com/aleimba/bac-genomics-scripts#citation), [installation](https://github.com/aleimba/bac-genomics-scripts#installation-recommendations), and [license](https://github.com/aleimba/bac-genomics-scripts#license) information please see the repository main [*README.md*](https://github.com/aleimba/bac-genomics-scripts/blob/master/README.md).

## Changelog

* v0.3 (30.01.2013)
    - additional info in POD
    - check if result files already exist and ask user what to do
    - changed script name from `ecoli_mlst_alleles.pl` to `ecoli_mlst.pl`
* v0.2 (20.10.2012)
    - included a POD
    - options with Getopt::Long
    - don't consider input *E. coli* genome query files, which are too big (set cutoff at 9 MB for a fasta *E. coli* file)
    - draft *E. coli* genomes can now be used as input query files
    - additional info in 'publicSTs.txt' now associated to found ST types in output
    - give text to STDOUT which files were created
    - new option **-c** to align the resulting allele sequences via *ClustalW*
* v0.1 (25.10.2011)
