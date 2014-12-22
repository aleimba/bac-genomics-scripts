#!/usr/bin/perl

use strict;
use warnings;

my $Usage = "\n".
   "\t###################################################################\n".
   "\t#  $0 fasta fas                                #\n". #$0 = program name
   "\t#  Searches for MLST alleles in E. coli genomes according to Mark #\n".
   "\t#  Achtman's scheme [Wirth et al., 2006] and gives out allele     #\n".
   "\t#  profile and allele sequences (from the allele input files).    #\n".
   "\t#  Uses nucmer from the MUMmer package: mummer.sourceforge.net    #\n".
   "\t#                                                                 #\n".
   "\t#  fasta = file extension of E. coli genome-fasta files.          #\n".
   "\t#  fas = file extension of allele fasta-files, e.g. 'fas'.        #\n".
   "\t#                                                                 #\n".
   "\t#  version 0.1                                        A Leimbach  #\n".
   "\t#  25.10.2011                              aleimba[at]gmx[dot]de  #\n".
   "\t###################################################################\n\n";

#print usage if -h|--h|--help is given as argument
my ($ext_seq, $ext_allele) = @ARGV;
if (!defined($ext_seq) || !defined($ext_allele)) {
    die $Usage;
} elsif ($ext_seq =~ m/-h/) {
    die $Usage;
}

#The MLST alleles from the Achtman scheme
my @alleles = ('adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA');

#read directory and look for corresponding files
my $dirname = ".";
my @genome_files;
my %allele_files;
opendir(DIR, $dirname) or die "Can't opendir $dirname: $!\n";
while (defined(my $file = readdir(DIR))) { #go through each file in the directory
    if ($file =~ m/^(\S+)\.$ext_seq$/ && (-s $file < 9000000)) { #don't use files > 10 MB, E. coli fasta file shouldn't be to big (normally around 5 MB), otherwise probably wrong fasta file in folder
        push (@genome_files, $file);
    } elsif ($file eq 'Ecoli_MLST_seq.fasta') { #don't use the later created allele multi fasta file 'E_coli_MLST_multi.fasta' (see below)!
    next;
    }
    foreach my $allele (@alleles) {
    if ($file =~ m/^$allele\S*\.$ext_allele$/i) {
        $allele_files{$allele} = $file;
    }
    }
}
closedir(DIR);
die "No E. coli genome fasta-files were found!\n" unless scalar @genome_files; #empty array in scalar context returns 0
die "No allele fasta-files were found!\n" unless scalar %allele_files; #empty hash in scalar context returns 0

#concatenate E. coli genome sequence files to one multi-fasta file for nucmer run
my @delete; #Files to be deleted after program run
my $genome_file = 'Ecoli_multi.fasta'; #result filename of the multi-fasta file
open(GENOMES, ">$genome_file") or die "Failed to create file '$genome_file': $!\n";
foreach (sort @genome_files) {
    open(FILE, "<$_") or die "Failed to open E. coli genome file $_: $!\n";
    while (<FILE>) {
    chomp;
    print GENOMES "$_\n";
    }
    close FILE;
}
close GENOMES;
if (-e $genome_file) {
    push (@delete, $genome_file);
} else {
    die "The multi-fasta E. coli genome file \'$genome_file\' could not be created: $!\n";
}

#run nucmer with allele sequences as references and E. coli multi-fasta file as query
my @created_files; #Used to print out files that are created
my %coord_files;
foreach (keys %allele_files) {
    system ("nucmer --prefix $_-all_Ecoli $allele_files{$_} $genome_file"); #seperate args not possible, probably because nucmer not a system command?!
    system ("show-coords -lTH $_-all_Ecoli.delta > $_-all_Ecoli.coords");
    if (-e "$_-all_Ecoli.coords") {
    $coord_files{$_} = "$_-all_Ecoli.coords";
    push (@delete, "$_-all_Ecoli.delta");
    push (@created_files, "\t\t\t$_-all_Ecoli.coords\n");
    } else {
    die "$_-all_Ecoli.coords could not be created: $!\n";
    }
}

#Parse E. coli accession#s and genome names from multi-fasta file with NCBI headers
my $genome_number = 0; #to control if lines in *.coords files correspond to genome number (in comparison to $lines, see sub parse_coords)
my %genome_desc; #hash to associate accession# with genome descriptions
open(GENOMES, "<$genome_file") or die "Failed to open file '$genome_file': $!\n";
while (<GENOMES>) {
    chomp;
    if ( m/^>/ ) { #Shorten NCBI fasta headers
        $genome_number++;
    $_ =~ s/>//;
    $_ =~ s/\|//g;
    $_ =~ s/gi[0-9]*(emb|gb|dbj|ref)//; #ref for refseq sequences
    $_ =~ s/\'//g;
    $_ =~ s/( DNA, complete genome| chromosome, complete genome|, complete genome| chromosome, complete sequence| complete genome)//;
        $_ =~ s/, strain (\S+)$//;
        $_ =~ s/\s/_/g;
    if ( m/(\w*\d*\.\d)_(\S+)$/ ) { #split shortened headers in accession-nr. ($1) and genome description ($2)
        $genome_desc{$1} = $2;
    }
    }
}
close GENOMES;

#Parse *.coords nucmer result files with the MLST alleles and the corresponding accession-nrs.
my @summary_err; #Array to print out error summary in error file
my @detail_err; #Array to print out the more detailed errors at the end of error file
my $error = 0; #switches to '1' if error is detected!
my %hash; #used as anonymous hash in sub parse_coords to store MLST alleles and accesssion-nrs.
my $allele_hashref; #%hash returned form subroutine as reference
foreach (sort keys %coord_files) { ###sort not needed, but MLST alleles are given ordered alphabetically into subroutine
    $allele_hashref = &parse_coord($coord_files{$_}, $_);
}

#Print the corresponding allele sequences (from the allele sequence files) and sequence type for each E. coli genome
my $ecoli_allele_seq = 'Ecoli_MLST_seq.fasta'; #result file with allele sequences for each E. coli genome in multi-fasta format
open(SEQ, ">$ecoli_allele_seq") or die "Failed to create file '$ecoli_allele_seq': $!\n";
my $ecoli_allele_profile = 'Ecoli_MLST_profile.txt'; #result file for the sequence type profile of each E. coli genome
open(PROFILE, ">$ecoli_allele_profile") or die "Failed to create file '$ecoli_allele_profile': $!\n";
print PROFILE ("Strain");
foreach (sort @alleles) {
    print PROFILE ("\t$_");
}
print PROFILE ("\n");
foreach my $acc (sort {lc $genome_desc{$a} cmp lc $genome_desc{$b}} keys %genome_desc) { #call hash %genome_desc by keys, but sort by values of the hash case-insensitively (all in lowercase, because Perl sorts lowercase after uppercase!)
    print SEQ (">$genome_desc{$acc}\n");
    print PROFILE ("$genome_desc{$acc}");
    foreach my $key (sort keys %$allele_hashref) {
    my $match = 0; #switches to 1 if a corresponding allele was found in the allele sequence files###but should be a given, since nucmer found it!
    if (defined $allele_hashref->{$key}->{$acc}) {
        open(ALLELE, "<$allele_files{$key}") or die "Failed to open allele file $allele_files{$key}: $!\n";
        while (my $line = <ALLELE>) { #Search for the corresponding ID line in the multi-fasta allele seq-file
        chomp $line;
        if ($line =~ m/^>$allele_hashref->{$key}->{$acc}$/i) {
            $match = 1; #corresponding allele was found###
            my $nextline = <ALLELE>; #don't need the ID line but the following sequence lines
            while ($nextline !~ m/^>/) { #Print sequence in file until another '>' is found
            chomp $nextline;
            print SEQ ("$nextline");
            $nextline = <ALLELE>; #jump to the next line
            }
        }
        }
        close ALLELE;
        if ($match == 0) { ###should be a given, since nucmer found it!
        push (@summary_err, "$allele_hashref->{$key}->{$acc} allele sequence could not be found!\n");
        unshift (@detail_err, "No $allele_hashref->{$key}->{$acc} allele sequence could be found for $genome_desc{$acc}!\n");
        }
        $allele_hashref->{$key}->{$acc} =~ s/^(\D+)(\d*)$/$2/; #get rid of allele description (e.g. ADK6 -> 6)
        print PROFILE ("\t$allele_hashref->{$key}->{$acc}");
    } else {
        print SEQ ("XXX");
        print PROFILE ("\t?");
        unshift (@detail_err, "$genome_desc{$acc} didn't give a hit with $key, marked by \'XXX\' in allele sequences and \'?\' in allele profile!\n");
        next;
    }
    }
    print SEQ ("\n\n"); #a blank line in front of each ID line
    print PROFILE ("\n");
}
close SEQ;
close PROFILE;
if (-e $ecoli_allele_seq && $ecoli_allele_profile) {
    push (@created_files, "\n\tResult files:\n\t\t\t$ecoli_allele_seq\t->Allele sequences!\n", "\t\t\t$ecoli_allele_profile\t->Allele profile!\n");
} else {
    die "The result files $ecoli_allele_seq and $ecoli_allele_profile could not be created: $!\n";
}

#Print errors in error file
if ($error == 1) {
    my $err_file = 'errors.txt'; #Write detected errors in the alignment to 'error.txt' file, see sub parse_coords
    open(ERR, ">$err_file") or die "Failed to create file '$err_file': $!\n";
    print ERR ("Error summary:\n");
    print ERR (@summary_err);
    print ERR ("\nDetailed error output:\n");
    print ERR (@detail_err);
    push (@created_files, "\t\t\t$err_file\t\t->Error file!\n");
    close ERR;
}

#Delete unneeded files
foreach my $goners (@delete) {
    unlink $goners or warn "Could not unlink $goners: $!";
}

#Print, which files have been created
print "\nThe following files were created:\n";
print "\tCoordinates of MLST alleles in each genome:\n";
print @created_files;


exit;


#############
#Subroutines#
#############

#Subroutine to parse *.coords nucmer result files with the MLST alleles and the corresponding accession#s
sub parse_coord {
    my ($file, $allele) = @_;
    my $lines = 0; #lines of respective *.coords file, if hit is missing or insensitive hits are present in comparison to $genome_number
    my $discarded = 0; #number of discarded hits/lines in coord file
    open (FILE, "$file") or die "Failed to open file '$file': $!\n"; #open allele-.coord file
    while (<FILE>) {
        chomp;
        $lines++;
        if ( m/\d\t\d+\t\d+\t\d+\t(\d+)\t(\d+)\t(\d+\.\d+)\t(\d+)\t\d+\t(\D+\d+)\t\S*\|(\w*\d+\.\d)\|$/ ) {
        #$1 = length of ref alignment, $2 = length of query alignment, $3 = identity percentage,
        #$4 = length of ref seq, $5 = allele-nr., $6 = accession-nr.
            if ($1 != $4) { #Write error to @detail_err, will be printed in error.txt-file later
        push (@detail_err, "$genome_desc{$6}, $5: Reference (allele) alignment doesn't have a correct length, therefore allele is not included in output!\n");
        $discarded++;
        next; #alignment hit is useless and will be skipped
            } elsif ($2 != $4) {
        push (@detail_err, "$genome_desc{$6}, $5: Query (genome) alignment doesn't have a correct length, therefore allele is not included in output!\n");
        $discarded++;
        next;
            } elsif ($3 != '100.00') {
        push (@detail_err, "$genome_desc{$6}, $5: Identity is not 100\%, therefore is not included in output!\n");
        $discarded++;
        next;
            }
        $hash{$allele}{$6} = $5; #$5 format,e.g. ADK15
#            print "hash $allele alleles $6: $hash{$allele}->{$6}\n"; #Print out to test the hash
    }
    }
    close FILE;
    #Error identifications for later print out in error.txt (see above)
    if ($discarded == 0) {
    if ($lines < $genome_number) {
        $error = 1;
        push (@summary_err, $genome_number - $lines, " $allele allele(s) missing!\n");
    }
    } elsif ($discarded > 0) {
    $error = 1;
    if (($lines - $discarded) == $genome_number) {
        push (@summary_err, "Total number of specific assigned $allele alleles is correct, but $discarded non-specific hit(s) discarded!\n");
    } elsif (($lines - $discarded) < $genome_number) {
        push (@summary_err, $genome_number - ($lines - $discarded), " $allele allele(s) missing and $discarded non-specific hit(s)!\n");
    }
    }
    return \%hash;
}
