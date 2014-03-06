#!/usr/bin/perl
use warnings;
use strict;
use autodie;
use Getopt::Long;

my $usage = << "USAGE";

  #####################################################################
  # $0 -i blast_hits.txt [-c|-s] (-t)              #
  #                                                                   #
  # Use the result file from 'prot_finder.pl', 'blast_hits.txt', to   #
  # create a comma-separated presence/abscence file,                  #
  # '*binary_matrix.csv'. This 'binary matrix' can be used e.g. for   #
  # iTOL (Interactive Tree Of Life, http://itol.embl.de/).            #
  # However, the organism names have to have the same names as in the #
  # tree file, thus manual adaptation, e.g. in Excel, might be needed.#
  # CAREFUL, organisms that didn't have a blastp hit, are obviously   #
  # not present in 'blast_hits.txt', and hence can't be included by   #
  # '$0'. Thus, add them manually to the result    #
  # file(s).                                                          #
  #                                                                   #
  # Mandatory options:                                                #
  # -i, -input       input 'blast_hits.txt' file                      #
  # -c, -collective  collective presence/abscence file for all query  #
  #                  proteins combined [default, excludes '-s']       #
  # -s, -separate    separate presence/abscence files for each query  #
  #                  protein [excludes '-c']                          #
  # Optional options:                                                 #
  # -h, -help        print usage                                      #
  # -t, -total       count total occurences of query proteins not     #
  #                  just presence/abscence binary                    #
  # -v, -version     print version number                             #
  #                                                                   #
  # version 0.5 updated: 05.03.2014                  Andreas Leimbach #
  # 25.10.2012                                  aleimba[at]gmx[dot]de #
  #####################################################################

USAGE
;

### Get options with Getopt::Long
my $blast_hits; # result file 'blast_hits.txt' from 'prot_finder.pl' as input
my $opt_collective; # a single presence/absence file for all queries
my $opt_separate; # separate presence/absence files for each query
my $opt_total; # count total occurences of query proteins not just presence/abscence binary
my $version = 0.5;
my ($opt_version, $opt_help);
GetOptions ('input=s' => \$blast_hits,
            'collective' => \$opt_collective,
            'separate' => \$opt_separate,
            'total' => \$opt_total,
            'version' => \$opt_version,
            'help|?' => \$opt_help);



### Print usage
if ($opt_help) {
    die $usage;
} elsif ($opt_version) {
    die "$0 $version\n";
} elsif (!$blast_hits) {
    die $usage, "### Fatal error: Option or argument for \'-i\' is missing!\n\n";
}



### Enforce mandatory option
if (!$opt_collective && !$opt_separate) {
    warn "### Warning: None of the mandatory options (\'-c\' or \'-s\') given. Forcing default option \'-c\'!\n";
    $opt_collective = 1;
} elsif ($opt_collective && $opt_separate) {
    die "\n### Fatal error: Both mandatory options (\'-c\' and \'-s\') given! Choose only one of the output options!\n";
}



### Parse the input 'blast_hits.txt' file; associate organism with query hit and store all query proteins in a separate array
my @queries; # store all query proteins
my $query = ''; # each query only once in @queries
my %hits; # hash-in-hash to associate organism with query hit
open (BLAST, "<$blast_hits");
my $header = <BLAST>; # skip header of 'blast_hits.txt'
while (<BLAST>) {
    my @line = split (/\t/, $_); # $line[0] = subject_organism, $line[4] = query_acc
    if ($opt_total) { # count total occurences of query proteins
        $hits{$line[0]}{$line[4]}++;
    } elsif (!$opt_total) { # count only presence/abscence of query proteins, i.e. only binary output (0 or 1)
        $hits{$line[0]}{$line[4]} = 1; # $line[0] = organism, [4] = query acc
    }
    $query = $line[4];
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
if ($opt_collective) { # collective binary file
    open (BINARY, ">binary_matrix.csv");
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
} elsif ($opt_separate) { # separate binary files
    foreach my $query (sort @queries) {
        open (BINARY, ">$query\_binary\_matrix.csv");
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
if ($opt_collective) {
    result('binary_matrix.csv'); # subroutine
} elsif ($opt_separate) {
    foreach my $query (sort @queries) {
        result("$query\_binary\_matrix.csv"); # subroutine
    }
}

exit;



###############
# Subroutines #
###############

### Subroutine to test if result file(s) exists and print result file(s) to STDOUT
sub result {
    my $file = shift;
    if (-e $file) {
        print "\nResult file \"$file\" has been created!\n";
    } else {
        die "\nResults could not be written to \"$file\", quitting program!\n";
    }
    return 1;
}
