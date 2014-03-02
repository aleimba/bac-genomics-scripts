#!/usr/bin/perl
use warnings;
use strict;

my $usage = "\n".
   "\t#############################################################\n".
   "\t# $0 blast_hits.txt (c|s)                #\n".
   "\t# Use the result file from blast_prot_finder.pl,            #\n".
   "\t# 'blast_hits.txt', to create a comma-separated binary      #\n".
   "\t# input file for iTOL.                                      #\n".
   "\t# However, the organism names have to have the same names   #\n".
   "\t# as in the tree file, thus manual copying might be needed. #\n".
   "\t# Option 'c' results in a collective file for all query     #\n".
   "\t# proteins, 's' in seperate files for each query protein.   #\n".
   "\t# 'c' is the default if none is given.                      #\n".
   "\t#                                                           #\n".
   "\t# version 0.1                              Andreas Leimbach #\n".
   "\t# 25.10.2012                          aleimba[at]gmx[dot]de #\n".
   "\t#############################################################\n\n";


### Print usage if -h|--h|--help is given as argument or input file not given
my ($blast_hits, $option) = @ARGV;
if (!defined $blast_hits) {
    die $usage;
} elsif ($blast_hits =~ /-h/) {
    die $usage;
} elsif (!defined $option) {
    $option = 'c';
}


### Parse the input 'blast_hits.txt' file
my @queries; # store all query proteins
my $query = ''; # 'blast_hits.txt' ordered by query proteins; each query only once in @queries
my %hits; # hash to associate organism with query hit
open (BLAST, "<$blast_hits") or die "Can't open file \'$blast_hits\': $!\n";
my $header = <BLAST>; # skip header of 'blast_hits.txt'
while (<BLAST>) {
    my @line = split (/\t/, $_);
    $hits{$line[0]}{$line[3]} = 1; # $line[0] = organism, [3] = query protein
    if ($query ne $line[3]) { # push each query only once in @queries
    push (@queries, $line[3]);
    }
    $query = $line[3];
}
close BLAST;


### Print binary data to a collective or seperate (for each query; as needed by iTOL) file(s)
if ($option eq 'c') { # collective binary file
    open (BINARY, ">binary_iTOL.txt") or die "Can't create file 'binary_iTOL.txt': $!\n";
    ### Print header
    print BINARY "Organism";
    foreach my $query (sort @queries) {
    print BINARY ",$query";
    }
    print BINARY "\n";
    ### Print data in file
    foreach my $organism (sort keys %hits) {
    print BINARY "$organism";
    foreach my $query (sort @queries) {
        if (defined $hits{$organism}->{$query}) {
        print BINARY ",$hits{$organism}->{$query}";
        } else {
        print BINARY ",0";
        }
    }
    print BINARY "\n";
    }
    close BINARY;
} elsif ($option eq 's') { # separate binary files
    foreach my $query (sort @queries) {
    open (BINARY, ">$query\_binary\_iTOL.txt") or die "Can't create file '$query\_binary\_iTOL.txt': $!\n";
    foreach my $organism (sort keys %hits) {
        print BINARY "$organism";
        if (defined $hits{$organism}->{$query}) {
        print BINARY ",$hits{$organism}->{$query}\n";
        } else {
        print BINARY ",0\n";
        }
    }
    close BINARY;
    }
}


### Print to STDOUT which file has been created
if (-e 'binary_iTOL.txt') {
    print "\nResult file 'binary_iTOL.txt' has been created!\n\n";
} else {
    die "\nResults could not be written to 'binary_iTOL.txt', quitting program!\n\n";
}


exit;
