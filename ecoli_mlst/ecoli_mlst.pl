#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

ecoli_mlst.pl                                       25-10-2011

=head1 SYNOPSIS

C<perl ecoli_mlst.pl -a fas -g fasta>

=head1 DESCRIPTION

The script searches for multilocus sequence type (MLST) alleles in
I<E. coli> genomes according to Mark Achtman's scheme with seven
house-keeping genes (I<adk>, I<fumC>, I<gyrB>, I<icd>, I<mdh>,
I<purA>, and I<recA>) [Wirth et al., 2006]. I<NUCmer> from the
L<I<MUMmer package>|http://mummer.sourceforge.net/> is used to
compare the given allele sequences to bacterial genomes via
nucleotide alignments.

Download the allele files (adk.fas ...) and the sequence type file
('publicSTs.txt') from this website:
                  http://mlst.warwick.ac.uk/mlst/dbs/Ecoli

To run C<ecoli_mlst.pl> include all I<E. coli> genome files (file
extension e.g. 'fasta'), all allele sequence files (file extension
'fas') and 'publicSTs.txt' in the current working directory. The
allele profiles are parsed from the created *.coord files and written
to a result file, plus additional information from the file
'publicSTs.txt'. Also, the corresponding allele sequences (obtained
from the allele input files) are concatenated for each I<E. coli>
genome into a result multi-fasta file. Option B<-c> can be used to
initiate an alignment for this multi-fasta file with
L<I<ClustalW>|http://www.clustal.org/clustal2/> (standard alignment
parameters; has to be in the C<$PATH> or change variable
C<$clustal_call>). The alignment fasta output file can be used
directly for L<I<RAxML>|http://sco.h-its.org/exelixis/web/software/raxml/index.html>. CAREFUL the Phylip alignment format from
I<ClustalW> allows only 10 characters per strain ID.

C<ecoli_mlst.pl> works with complete and draft genomes. However,
several genomes cannot be included in a single input file!

Obviously, only for those genomes whose allele sequences have been
deposited in Achtman's allele database results can be obtained. If an
allele is not found in a genome it is marked by a '?' in the result
profile file and a place holder 'XXX' in the result fasta file. For
these cases a manual I<NUCmer> or I<BLASTN> might be useful to fill the
gaps and L<C<run_sub_seq.pl>|/run_sub_seq> to get the corresponding
'new' allele sequences.

Non-NCBI fasta headers for the genome files have to have a unique ID
directly following the '>' (e.g. 'Sakai', '55989' ...).

=head1 OPTIONS

=head2 Mandatory options

=over 22

=item B<-a>=I<str>, B<-alleles>=I<str>

File extension of the MLST allele fasta files, e.g. 'fas' (<=> B<-g>).

=item B<-g>=I<str>, B<-genomes>=I<str>

File extension of the I<E. coli> genome fasta files, e.g. 'fasta' (<=> B<-a>).

=back

=head2 Optional options

=over

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-c>, B<-clustalw>

Call L<I<ClustalW>|http://www.clustal.org/clustal2/> for alignment

=back

=head1 OUTPUT

=over 17

=item F<ecoli_mlst_profile.txt>

Tab-separated allele profiles for the I<E. coli> genomes, plus additional info from 'publicSTs.txt'

=item F<ecoli_mlst_seq.fasta>

Multi-fasta file of all concatenated allele sequences for each genome

=item F<*.coord>

Text files that contain the coordinates of the I<NUCmer> hits for each genome and allele

=item (F<errors.txt>)

Error file, summarizing number of not found alleles or unclear I<NUCmer> hits

=item (F<ecoli_mlst_seq_aln.fasta>)

Optional, L<I<ClustalW>|http://www.clustal.org/clustal2/> alignment in Phylip format

=item (F<ecoli_mlst_seq_aln.dnd>)

Optional, I<ClustalW> alignment guide tree

=back

=head1 EXAMPLE

=over

=item C<perl ecoli_mlst.pl -a fas -g fasta -c>

=back

=head1 VERSION

0.3                                                update: 30-01-2013

=head1 AUTHOR

Andreas Leimbach                                aleimba[at]gmx[dot]de

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 (GPLv3) of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see L<http://www.gnu.org/licenses/>.

=cut


########
# MAIN #
########

use strict;
use warnings;
use Getopt::Long; # module to get options from the command line


### Get the options with Getopt::Long, works also abbreviated and with two "--": -g, --g, -genomes ...
my $usage = "perldoc $0";
die system($usage) unless @ARGV;
my $allele_ext = ''; # file extension of allele fasta files (<=> $genome_ext)
my $genome_ext = ''; # file extension of E. coli genome fasta files (<=> allele_ext)
my $clustalw = ''; # optionally, call ClustalW for alignment
my $help = ''; # run perldoc on the POD
GetOptions ('alleles=s' => \$allele_ext, 'genomes=s' => \$genome_ext, 'clustalw' => \$clustalw, 'help' => \$help);


### Run perldoc on POD if no arguments or help
if (!$genome_ext || !$allele_ext) {
    die system($usage);
} elsif ($help) {
    die system($usage);
}


### Check if result files already exist, overwrite or exit script after user question
my $ecoli_allele_seq = 'ecoli_mlst_seq.fasta'; # multi-fasta result file with concatenated allele seqs for each E. coli genome (sequences are from the MLST allele files)
my $clustal_aln = 'ecoli_mlst_seq_aln.fasta'; # fasta alignment file from optional clustalW alignment
if (-e $ecoli_allele_seq) {
    print "A previous analysis exists, overwrite the old result files [y|n]? ";
    my $stdin = <STDIN>;
    chomp $stdin;
    if ($stdin =~ /n/i) {
         die "Script abborted!\n\n";
    } elsif ($stdin =~ /y/i) {
         if (-e $clustal_aln) { # get rid of the optional clustalW alignment fasta-file before program run
             unlink $clustal_aln;
         }
    }
}


### The MLST alleles from the Achtman scheme
my @alleles = ('adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA');


### Read directory and look for corresponding files
my $dirname = ".";
my @genome_files; # array to save all the fasta genome files
my %allele_files; # hash to save all the multi-fasta allele files
opendir(DIR, $dirname) or die "Can't opendir $dirname: $!\n";
while (defined(my $file = readdir(DIR))) { # go through each file in the directory
    if ($file eq 'ecoli_mlst_seq.fasta') { # don't use the result file from a previous analysis, the allele multi fasta file 'ecoli_mlst_multi.fasta' (see below)!
    next;
    } elsif ($file =~ /^(\S+)\.$genome_ext$/ && (-s $file < 9000000)) { # don't use files > 10 MB, E. coli fasta file shouldn't be too big (normally around 5 MB), otherwise probably wrong fasta file in folder
        push (@genome_files, $file);
    }
    foreach my $allele (@alleles) {
    if ($file =~ /^$allele\S*\.$allele_ext$/i) {
        $allele_files{$allele} = $file;
    }
    }
}
closedir(DIR);
die "No E. coli genome fasta-files were found!\n" unless scalar @genome_files; # empty array in scalar context returns 0
die "No allele fasta-files were found!\n" unless scalar %allele_files; # empty hash in scalar context returns 0


### Look for multi-fasta genome files and concatenate them to a single-fasta entry (e.g. draft genomes). Subsequently concatenate all E. coli genome sequence files to one multi-fasta/-genome file for the subsequent nucmer run. Additionally parse and associate the E. coli accession#s/unique IDs and descriptions (strain names ...) from the genomes
my $genome_file = 'ecoli_multi.fasta'; # multi-genome/-fasta file for the nucmer run, used as a temp file, will be deleted at the end of the script
open(GENOMES, ">$genome_file") or die "Failed to create file '$genome_file': $!\n";
my $genome_number = 0; # used to control if lines in *.coords files correspond to the overall genome count
my %genome_desc; # hash to associate accession#/unique ID with genome descriptions
foreach my $file (@genome_files) {
    $genome_number++;
    open(MULTI, "<$file") or die "Failed to open E. coli genome file '$file': $!\n";
    my @multi = <MULTI>; # read the whole genome file (potential multi-fasta file)
    my @IDs = grep(/^>/, @multi); # get ID lines in fasta file
    if (scalar @IDs > 1) { # multi-fasta files
    my $new_id; # new ID line for the multi-fasta file
    foreach (@IDs) { # discard plasmid ID lines in complete genomes for new file (name with chromosome ID line)
        if (!/plasmid/) {
        chomp;
        $new_id = $_;
        last; # jump out of the loop if a non-plasmid ID found
        }
    }
    if (!defined($new_id)) {
        die "The file '$file' only contains plasmids, program exits!\n";
    }
    my ($acc, $desc) = acc_desc($new_id); # subroutine to fill hash %genome_desc with acc#/unique ID and genome description
    print GENOMES ">$acc $desc; artificial genome\n"; # print the now shortened ID-line for the new single-fasta entry in the result multi-genome file
    foreach (@multi) { # print the rest of the new single-fasta file
        if (/^>/) { # skip ID lines of the old multi-fasta file
        next;
        }
        print GENOMES;
    }
    } elsif (scalar @IDs == 1) { # non-multi-fasta files
    my ($acc, $desc) = acc_desc($IDs[0]);
    print GENOMES ">$acc $desc\n"; # print the new shortened ID-line
    foreach (@multi) {
        if (/^>/) { # skip ID line, since new one already printed
        next;
        }
        print GENOMES;
    }
    } else { # wrong fasta files
    die "File '$file' does not include a fasta ID line, program exits!\n";
    }
    print GENOMES "\n";
    close MULTI;
}
close GENOMES;
my @delete; # files to be deleted after program run
if (-e $genome_file) {
    push (@delete, $genome_file);
} else {
    die "The multi-fasta E. coli genome file \'$genome_file\' could not be created for the NUCmer run: $!\n";
}


### Parse file 'publicSTs.txt' to get additional info for each sequence type
my $st_file = 'publicSTs.txt';
open(ST, "<$st_file")  or die "Failed to open sequence type file '$st_file': $!\n";
my %seq_type; # hash to save ST info to each allele profile
my $header = <ST>; # get rid of file header
while (<ST>) {
    chomp;
    my @st = split (/\t/, $_);
    my $profile = "$st[3] $st[4] $st[5] $st[6] $st[7] $st[8] $st[9]"; # allele profile
    my $info = "$st[0]\t$st[1]\t$st[2]\t$st[10]\t$st[11]"; # associated info to allele profile
    # $st[0] = ST, 1 = ST_COMPLEX, 2 = ANCESTRAL_GROUP, 3-9 = adk-recA alleles, 10 = SOURCE, 11 = REFSTRAIN
    $seq_type{$profile} = $info;
}
close ST;


### Run nucmer with allele sequences as REFERENCES and the created 'ecoli_multi.fasta' file as QUERY
my @created_files; # used to print out files that are created at the end of the script
my %coord_files;
foreach (keys %allele_files) {
    system ("nucmer --prefix $_-all_ecoli $allele_files{$_} $genome_file"); # system call; seperate args not possible, probably because nucmer not a system command?!
    system ("show-coords -lT $_-all_ecoli.delta > $_-all_ecoli.coords"); # '-l' include seq length in output, '-T' tab-separated output
    if (-e "$_-all_ecoli.coords") {
    $coord_files{$_} = "$_-all_ecoli.coords";
    push (@delete, "$_-all_ecoli.delta");
    push (@created_files, "\t\t\t$_-all_ecoli.coords\n");
    } else {
    die "Coord file '$_-all_ecoli.coords' could not be created: $!\n";
    }
}


### Parse *.coords nucmer result files for MLST alleles and corresponding accession#s
my @summary_err; # array to store error summary for error file
my @detail_err; # array to store the more detailed errors for error file
my $error = 0; # switches to '1' if an error is detected
my %strain_allele; # declare hash, that's subsequently used as 'hash in hash' in sub parse_coords to associate MLST alleles and accession#s
foreach (sort keys %coord_files) {
    parse_coord($coord_files{$_}, $_); # subroutine to fill %strain_allele
}


### Print the corresponding allele sequences and allele profile (plus additional info) for each E. coli genome
open(SEQ, ">$ecoli_allele_seq") or die "Failed to create file '$ecoli_allele_seq': $!\n";
my $ecoli_allele_profile = 'ecoli_mlst_profile.txt'; # tab-separated result file for the allele profile of each E. coli genome plus additional info from file 'publicSTs.txt' from the Achtman MLST DB
open(PROFILE, ">$ecoli_allele_profile") or die "Failed to create file '$ecoli_allele_profile': $!\n";
print PROFILE "Strain"; # header for profile file
foreach (sort @alleles) {
    print PROFILE "\t$_";
}
print PROFILE "\tST\tST_COMPLEX\tANCESTRAL_GROUP\tSOURCE\tREFSTRAIN\n"; # info from file 'publicSTs.txt'
foreach my $acc (sort {lc $genome_desc{$a} cmp lc $genome_desc{$b}} keys %genome_desc) { # call hash %genome_desc by keys (acc#s), but sort by values (genome desc.s) of the hash case-insensitively (all in lowercase, because Perl sorts lowercase after uppercase!)
    print SEQ ">$genome_desc{$acc} $acc\n"; # print fasta ID line for 'ecoli_mlst_seq.fasta'
    print PROFILE "$genome_desc{$acc}"; # print genome desc for strain in first column for 'ecoli_mlst_profile.txt'
    my $profile = ''; # allele profile to extract info from %seq_type
    foreach my $allele (sort keys %strain_allele) {
    if (defined $strain_allele{$allele}->{$acc}) {
        open(ALLELE, "<$allele_files{$allele}") or die "Failed to open allele file $allele_files{$allele}: $!\n";
        while (my $line = <ALLELE>) { # search for corresponding allele in multi-fasta allele file
        chomp $line;
        if ($line =~ /^>$strain_allele{$allele}->{$acc}$/i) {
            $line = <ALLELE>; # don't need the ID line but the following seq lines
            while ($line !~ /^>/) { # print sequence in result file until another '>' is found
            chomp $line;
            print SEQ "$line";
            $line = <ALLELE>;
            }
        }
        }
        close ALLELE;
        $strain_allele{$allele}->{$acc} =~ s/^(\D+)(\d+)$/$2/; # only keep allele number
        print PROFILE "\t$strain_allele{$allele}->{$acc}";
        $profile .= "$strain_allele{$allele}->{$acc} "; # concat allele numbers in $profile to extract info from %seq_type
    } else {
        print SEQ "XXX"; # place holder if allele of genome not determined
        print PROFILE "\t?"; # place holder as well
        unshift (@detail_err, "$genome_desc{$acc} didn't give a hit with $allele, marked by \'XXX\' in allele sequences and \'?\' in allele profile!\n");
        next;
    }
    }
    print SEQ "\n\n"; # blank line in front of each ID line
    chop $profile; # get rid of the last space, so it can be used as the key in %seq_type
    if (defined $seq_type{$profile}) { # print ST and additional info from file 'publicST.txt' in 'ecoli_mlst_profile.txt'
    print PROFILE "\t$seq_type{$profile}\n";
    } else {
    print PROFILE "\t?\t?\t?\t?\t?\n";
    }
}
close SEQ;
close PROFILE;
if (-e $ecoli_allele_seq && $ecoli_allele_profile) { # push new created files in array for print out
    push (@created_files, "\n\tResult files:\n\t\t\t$ecoli_allele_seq\t-> Allele sequences!\n", "\t\t\t$ecoli_allele_profile\t-> Allele profile plus info from 'publicSTs.txt'!\n");
} else {
    die "The result files $ecoli_allele_seq and $ecoli_allele_profile could not be created: $!\n";
}


### Print errors in error file,
if ($error == 1) { # error(s) occurred
    my $err_file = 'errors.txt';
    open(ERR, ">$err_file") or die "Failed to create file '$err_file': $!\n";
    print ERR "Error summary:\n";
    print ERR @summary_err; # filled in sub 'parse_coords'
    print ERR "\nDetailed error output:\n";
    print ERR @detail_err;
    push (@created_files, "\t\t\t$err_file\t\t-> Error file!\n");
    close ERR;
}


### Delete unneeded files, temp file 'ecoli_multi.fasta' (see above) and the *.delta files from the NUCmer run
foreach my $goners (@delete) {
    unlink $goners or warn "Could not remove file '$goners': $!";
}


### Print newly created files
print "\n###########################################################################\n";
print "The following files were created:\n";
print "\tCoordinates of MLST alleles in each genome:\n";
print @created_files;
print "\n###########################################################################\n";


### Align with ClustalW if option '-c' is given
if ($clustalw) {
    print "\nStarting ClustalW alignment with file $ecoli_allele_seq ...";
    my $out = $ecoli_allele_seq;
    $out =~ s/\.fasta$//;
    my $clustal_call = "clustalw -infile=$ecoli_allele_seq -outfile=$clustal_aln -align -type=DNA -output=phylip -tree -newtree=$out\_aln.dnd -outputtree=phylip";
    system ($clustal_call);
}

exit;


###############
# Subroutines #
###############

### Subroutine to associate the acc# to the genome description in hash %genome_desc (see above)
sub acc_desc {
    my $ID = shift;
    chomp $ID;
    if ($ID =~ /^>gi\|\d+\|(emb|gb|dbj|ref)\|(\w+\d+\.\d)\|\s(.+)$/) { # NCBI fasta header
    # e.g.: >gi|387605479|ref|NC_017626.1| Escherichia coli 042, complete genome
    # $1 = DB (embl, genbank, ddbj, refseq), $2 = acc#, $3 = genome description
    my $desc = shorten_desc($3); # subroutine to shorten the genome description
    $genome_desc{$2} = $desc; # associate accession# with genome description
    return ($2, $desc);
    } elsif ($ID =~ /^>(\S+)\s*(\S*)\s*(.*)$/) { # headers after EMBOSS's union of multi-fasta files, and other headers with a unique ID after '>'
    # e.g.: >NC_011748.1 NC_011748.1 Escherichia coli 55989 chromosome, complete genome
    # $1 = acc#, $2 = genome desc for drafts|duplicated acc# for complete genomes, $3 = genome desc for complete genomes
    my $desc = $2 . ' ' . $3;
    if ($1 eq $2) { # complete genomes have acc# double after EMBOSS's union (see above)
        $desc = $3;
    }
    $desc = shorten_desc($desc);
    $genome_desc{$1} = $desc;
    return ($1, $desc);
    }
    return 0;
}


### Subroutine to parse *.coord nucmer result files for MLST alleles and corresponding accession#s
sub parse_coord {
    my ($coord_file, $allele) = @_;
    my $lines = 0; # lines of respective *.coords file (without header); control if hit is missing or insensitive hits are present in comparison to $genome_number for the error file
    my $discarded = 0; # number of discarded hits/lines in coord file for error file
    open (COORD, "$coord_file") or die "Failed to open file '$coord_file': $!\n";
    while (<COORD>) {
        chomp;
        if ( /\d+\t\d+\t\d+\t\d+\t(\d+)\t(\d+)\t(\d+\.\d+)\t(\d+)\t\d+\t(\D+\d+)\t(\w*\d+\.\d)$/ ) {
        # since the ID lines of the fastas were shortened in temp file 'ecoli_multi.fasta', the acc#s should be the last element of each line in the coords file
        # $1 = length of ref alignment, $2 = length of query alignment, $3 = identity percentage,
        # $4 = length of ref seq, $5 = allele-nr. (e.g. ADK15), $6 = accession-nr.
        $lines++; # after 'if' to skip headers
            if ($1 != $4) { # write error to @detail_err for error file
        push (@detail_err, "$genome_desc{$6}, $5: Reference (allele) alignment doesn't have a correct length, therefore allele is not included in output!\n");
        $discarded++;
        next;
            } elsif ($2 != $4) {
        push (@detail_err, "$genome_desc{$6}, $5: Query (genome) alignment doesn't have a correct length, therefore allele is not included in output!\n");
        $discarded++;
        next;
            } elsif ($3 != '100.00') {
        push (@detail_err, "$genome_desc{$6}, $5: Identity is not 100\%, therefore is not included in output!\n");
        $discarded++;
        next;
            }
        $strain_allele{$allele}{$6} = $5; # associate allele# with acc# and store as hash in hash in %strain_allele
    }
    }
    close COORD;
    # error identifications for later print out in error.txt
    if ($discarded == 0) {
    if ($lines < $genome_number) {
        $error = 1; # switch $error from 0 to 1 to indicate error was found
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
    return 1;
}


### Shorten the genome descriptions of the ID headers
sub shorten_desc {
    my $desc = shift;
    $desc =~ s/Escherichia coli/Ecoli/;
    $desc =~ s/\'//g;
    $desc =~ s/( DNA, complete genome| chromosome, complete genome|, complete genome| chromosome, complete sequence| complete genome|, complete sequence|, strain (\S+)|, whole genome shotgun sequence)$//;
    $desc =~ s/\s/_/g;
    return $desc;
}
