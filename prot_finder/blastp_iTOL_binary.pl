#!/usr/bin/perl
use warnings;
use strict;

my $usage = "\n".
   "\t#####################################################################\n".
   "\t# $0 blast_hits.txt [c|s]                        #\n".
   "\t#                                                                   #\n".
   "\t# Use the result file from 'blast_prot_finder.pl', 'blast_hits.txt',#\n".
   "\t# to create a comma-separated binary input file for iTOL            #\n".
   "\t# (Interactive Tree Of Life, http://itol.embl.de/). However, the    #\n".
   "\t# organism names have to have the same names as in the tree file,   #\n".
   "\t# thus manual adaptation, e.g. in Excel, might be needed. CAREFUL,  #\n".
   "\t# organisms that didn't have a blastp hit, are obviously not        #\n".
   "\t# present in 'blast_hits.txt', and hence can't be included by       #\n".
   "\t# 'blastp_iTOL_binary.pl'. Thus, add them manually to the result    #\n".
   "\t# file(s).                                                          #\n".
   "\t# Option 'c' results in a collective file for all query proteins,   #\n".
   "\t# 's' in separate files for each query protein. Although iTOL       #\n".
   "\t# needs a separate file for each binary dataset, 'c' is the default #\n".
   "\t# if none is given, as this makes adapting in Excel easier.         #\n".
   "\t#                                                                   #\n".
   "\t# version 0.3 updated: 15.02.2013                  Andreas Leimbach #\n".
   "\t# 25.10.2012                                  aleimba[at]gmx[dot]de #\n".
   "\t#####################################################################\n\n";


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
my $query = ''; # each query only once in @queries
my %hits; # hash-in-hash to associate organism with query hit
open (BLAST, "<$blast_hits") or die "Can't open file \'$blast_hits\': $!\n";
my $header = <BLAST>; # skip header of 'blast_hits.txt'
while (<BLAST>) {
    my @line = split (/\t/, $_);
    $hits{$line[0]}{$line[3]} = 1; # $line[0] = organism, [3] = query protein
    $query = $line[3];
    foreach (@queries) { # control if the query already exists
        if ($_ eq $query) {
            $query = '';
        }
    }
    if ($query ne '') { # push each query only once in @queries
        push (@queries, $query);
    }
}
close BLAST;


### Print binary data to a collective or separate (for each query; as needed by iTOL) file(s)
if ($option eq 'c') { # collective binary file
    open (BINARY, ">binary_iTOL.csv") or die "Can't create file 'binary_iTOL.txt': $!\n";
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
        open (BINARY, ">$query\_binary\_iTOL.csv") or die "Can't create file '$query\_binary\_iTOL.txt': $!\n";
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
if ($option eq 'c') {
    result('binary_iTOL.csv'); # subroutine
} elsif ($option eq 's') {
    foreach my $query (sort @queries) {
        result("$query\_binary\_iTOL.csv"); # subroutine
    }
}

exit;


###############
# Subroutines #
###############

### Subroutine to test if result file(s) exists and print result files to STDOUT
sub result {
    my $file = shift;
    if (-e $file) {
        print "\nResult file \"$file\" has been created!\n";
    } else {
        die "\nResults could not be written to $file, quitting program!\n";
    }
    return 1;
}
