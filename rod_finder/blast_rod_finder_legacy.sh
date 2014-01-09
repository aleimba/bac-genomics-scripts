#!/bin/bash
echo "###Running legacy BLASTN"
formatdb -p F -i $1 -n ROD
blastall -p blastn -d ROD -i $2 -o blastn.out -e 2e-11 -F F
echo "###Running blast_rod_finder.pl"
perl blast_rod_finder.pl -q $3 -r blastn.out -m $4