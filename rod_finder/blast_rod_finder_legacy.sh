#!/bin/bash
echo "###Running legacy BLASTN with subject '$1' and query '$2'"
formatdb -p F -i $1 -n ROD
blastall -p blastn -d ROD -i $2 -o blastn.out -e 2e-11 -F F
echo "###Running blast_rod_finder.pl with query '$3' and minimum ROD size '$4'"
perl blast_rod_finder.pl -q $3 -r blastn.out -m $4