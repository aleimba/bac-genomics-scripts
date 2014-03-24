#!/bin/bash
# Download/update RefSeq complete genomes
echo "#### Updating RefSeq complete $1 genomes"
wget -c -N -v -r -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/$1*" -P .
# Download/update RefSeq draft genomes
echo "#### Updating RefSeq draft $1 genomes"
wget -c -N -v -r -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/$1*" -P .
# Download/update Genbank complete genomes
echo "#### Updating Genbank complete $1 genomes"
wget -c -N -v -r -t 45 -A *.gbk,*.fna "ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/$1*" -P .
# Download/update Genbank draft genomes
echo "#### Updating Genbank draft $1 genomes"
wget -c -N -v -r -t 45 -A *.gbk.tgz,*.fna.tgz "ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria_DRAFT/$1*" -P .
# Run script 'ncbi_concat_unpack.pl' to fill the result folders './refseq' and './genbank'
echo "#### Copying files to result folder './refseq'"
perl ncbi_ftp_concat_unpack.pl refseq y
echo "#### Copying files to result folder './genbank'"
perl ncbi_ftp_concat_unpack.pl genbank y
