#!/bin/bash
# Download/update RefSeq complete genomes
echo "#### Updating RefSeq complete $1 genomes"
wget -cNrv -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/$1*" -P .
# Download/update RefSeq draft genomes
echo "#### Updating RefSeq draft $1 genomes"
wget -cNrv -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT/$1*" -P .
# Download/update GenBank complete genomes
echo "#### Updating GenBank complete $1 genomes"
wget -cNrv -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/$1*" -P .
# Download/update GenBank draft genomes
echo "#### Updating GenBank draft $1 genomes"
wget -cNrv -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/$1*" -P .
# Run script 'ncbi_concat_unpack.pl' to fill the result folders './refseq' and './genbank'
echo "#### Copying files to result folder './refseq'"
perl ncbi_ftp_concat_unpack.pl refseq y
echo "#### Copying files to result folder './genbank'"
perl ncbi_ftp_concat_unpack.pl genbank y
