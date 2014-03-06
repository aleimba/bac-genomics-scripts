#!/bin/bash
echo "###Running cds_extractor.pl on all '$1'-files in the current working directory"
for i in *.$1; do perl cds_extractor.pl -i $i -p -f; done
cat *_cds_aa.fasta > blast_finder_aa.fasta
rm -f *_cds_aa.fasta
echo "###Running legacy BLASTP with query '$2'"
formatdb -p T -i blast_finder_aa.fasta -n blast_finder
blastall -p blastp -d blast_finder -i $2 -o blastp.out -e 1e-10 -F F -s T -b 500
echo "###Running blast_prot_finder.pl with an identity cutoff of '$3' and a query coverage cutoff of '$4' plus Clustal Omega alignment"
perl blast_prot_finder.pl -q $2 -s blast_finder_aa.fasta -r blastp.out -i $3 -cov_q $4 -a -b
