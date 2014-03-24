#!/usr/bin/perl
use warnings;
use strict;
use File::Find; # module to traverse directory trees

my $usage = "\n".
   "\t#####################################################################\n".
   "\t# $0 genbank|refseq (y|n)                    #\n".
   "\t#                                                                   #\n".
   "\t# Unpacks and concatenates all draft and complete genomes (in       #\n".
   "\t# genbank and fasta format) downloaded from the NCBI ftp server     #\n".
   "\t# (ftp://ftp.ncbi.nih.gov/). The script traverses the dowloaded     #\n".
   "\t# NCBI ftp-folder structure and thus has to be called from the      #\n".
   "\t# top level (containing the folder 'ftp.ncbi.nih.gov').             #\n".
   "\t# For complete genomes PLASMIDS are concatenated to the CHROMOSOMES #\n".
   "\t# to create multi-genbank files. In draft genomes, scaffold and     #\n".
   "\t# contig files are controlled for annotation. The ones with         #\n".
   "\t# annotation are then used to create multi-genbank files.           #\n".
   "\t# Multi-fasta files are created for the corresponding genbank-      #\n".
   "\t# file or if no annotation exists for the file which contains       #\n".
   "\t# more sequence information (either contigs or scaffolds).          #\n".
   "\t# Use option 'genbank' to use Genbank genomes and 'refseq' for      #\n".
   "\t# RefSeq genomes. Concatenated files are stored in the result       #\n".
   "\t# folders './genbank' and './refseq' respectively. Set yes/no       #\n".
   "\t# in the program call to delete current result folder and           #\n".
   "\t# create a new one.                                                 #\n".
   "\t#                                                                   #\n".
   "\t# version 0.1                                      Andreas Leimbach #\n".
   "\t# 15.09.2012                                  aleimba[at]gmx[dot]de #\n".
   "\t#####################################################################\n\n";


### Print usage if -h|--h|--help is given as argument or options are not given
my ($db, $ask) = @ARGV; # $ask = if not given ask if ok to proceed to delete the current result folder and write new
if (!defined $db) {
    die $usage;
} elsif ($db =~ m/-h/) {
    die $usage;
}


### Set the correct directories to start traversion of the folder structure
my $dir_complete; # for complete genomes
my $dir_draft; # for draft genomes
if ($db =~ /genbank/i) {
    $db = 'genbank';
    $dir_complete = './ftp.ncbi.nih.gov/genbank/genomes/Bacteria';
    $dir_draft = './ftp.ncbi.nih.gov/genbank/genomes/Bacteria_DRAFT';
    ask_exit($db, $ask); # subroutine to ask if the previous result folder should be deleted or exit the script
    rm_dir($db); # subroutine to remove the result directory and all its contents prior filling it with new files (in case files have changed)
} elsif ($db =~ /refseq/i) {
    $db = 'refseq';
    $dir_complete = './ftp.ncbi.nih.gov/genomes/Bacteria';
    $dir_draft = './ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT';
    ask_exit($db, $ask);
    rm_dir($db);
}


### Concatenate all genbank and fasta files for complete genomes in subdirectories and write to the genbank/refseq folders
find ({wanted => \&completes, no_chdir => 1}, $dir_complete); # the function 'find' from File::Find takes a reference to a subroutine and a starting directory to walk through recursively (e.g. 'find (\$concat, $dirname)')
# For each file/directory/link in the tree (also in the current directory) it calls the referenced function and changes to that directory with 'chdir()' --> this can be stopped with 'no_chdir => 1', parameters to find are passed as a hash reference


### Unzip and unpack all draft files with *.tgz in subdirectories and concat in result folders
my %draft_files; # store the filenames of the concatenated files in here
find ({wanted => \&drafts, no_chdir => 1}, $dir_draft);


### Keep only gbk files with annotation (either only in contig, scaffold) and the corresponding fastas for each strain, if no annotation exists keep the fasta file with the most sequence information (contig or scaffold), delete the other files
my ($A, $C, $G, $T); # to calculate the total base count (subs 'seq_size' and 'base_count')
my $skip = ''; # not possible to delete elements from a hash while iterating through it, thus need a variable to skip files according to a regex
foreach my $file (sort{$b cmp $a} keys %draft_files) { # reverse sort to get scaf files first if existent
    $skip =~ s/(\S+)\_draft\_(scaf|con)\.\w*$/$1/; # get rid of file-extension
    if ($file =~ /$skip/) { # skip all files that fit to regex, if result file was determined previously
        next;
    }
    if ($file =~ /\_draft\_scaf\.gbk$/) { # if scaffold files exist
        my $scaf_gbkfile = my $scaf_fasfile =$file;
        $scaf_fasfile =~ s/gbk/fasta/;
        my $con_gbkfile = $scaf_gbkfile;
        $con_gbkfile =~ s/\_scaf\./\_con\./;
        my $con_fasfile = $con_gbkfile;
        $con_fasfile =~ s/gbk/fasta/;
        my $scaf_size = -s "./$db/$scaf_gbkfile"; # file size of the gbk file
        my $con_size = -s "./$db/$con_gbkfile";
        my $scaf_gene = count_gene($scaf_gbkfile); # subroutine to count genes in genbank files
        my $con_gene = count_gene($con_gbkfile);
        if ($scaf_gene > $con_gene && $scaf_size > $con_size && $con_gene == 0) { # So many conditions needed? Do they always work?
            unlink "./$db/$con_gbkfile" or die "Can't remove file \'$con_gbkfile\': $!\n"; # remove the contig-gbk without annotation
            $skip = $file; # scaffold file with annotation found, skip all further files in %draft_files
            my $scaf_gbkcount = seq_size($scaf_gbkfile); # subroutine to calculate total bases (A, C, G, and T) of genbank
            unlink "./$db/$con_fasfile" or die "Can't remove file \'$con_fasfile\': $!\n"; # also delete the corresponding fasta
            my $scaf_fascount = seq_size($scaf_fasfile); # count total bases of fasta
            warn_seq($scaf_fasfile, $scaf_gbkcount, $scaf_fascount); # subroutine to warn if the genbank and the fasta file have differing sequence sizes
            next;
        } elsif ($con_gene > $scaf_gene && $con_size > $scaf_size && $scaf_gene == 0) {
            unlink "./$db/$scaf_gbkfile" or die "Can't remove file \'$scaf_gbkfile\': $!\n";
            $skip = $file;
            my $con_gbkcount = seq_size($con_gbkfile);
            unlink "./$db/$scaf_fasfile" or die "Can't remove file \'$scaf_fasfile\': $!\n";
            my $con_fascount = seq_size($con_fasfile);
            warn_seq($con_fasfile, $con_gbkcount, $con_fascount);
            next;
        } elsif ($con_gene == 0 && $scaf_gene == 0) { # no annotation in both genbanks, just keep one fasta
            unlink "./$db/$scaf_gbkfile" or die "Can't remove file \'$scaf_gbkfile\': $!\n";
            unlink "./$db/$con_gbkfile" or die "Can't remove file \'$con_gbkfile\': $!\n";
            $skip = $file;
            my $scaf_fascount = seq_size($scaf_fasfile);
            my $con_fascount = seq_size($con_fasfile);
            if ($scaf_fascount > $con_fascount || $scaf_fascount == $con_fascount) { # keep scaffold file if both same base count
                unlink "./$db/$con_fasfile" or die "Can't remove file \'$con_fasfile\': $!\n";
            } elsif ($con_fascount > $scaf_fascount) {
                unlink "./$db/$scaf_fasfile" or die "Can't remove file \'$scaf_fasfile\': $!\n";
            }
            next;
        }
    } elsif ($file =~ /\_draft\_con\.gbk$/) { # no scaffold files exist
        my $con_gene = count_gene($file);
        $skip = $file;
        if ($con_gene == 0) { # no annotation exists, only keep the fasta
            unlink "./$db/$file" or die "Can't remove file \'$file\': $!\n";
        }
    }
}


### Print result folder
print "All genome files have been written to the folder \'./$db\'!\n";

exit;


###############
# Subroutines #
###############

### Subroutine to ask if the script should proceed with deleting the previous result folder or exit
sub ask_exit {
    my ($db, $ask) = @_;
    if (!defined($ask)) { # $ask can also be defined as an option in the script call
        print "The script will delete an already existing result folder \'./$db\', proceed [y|n]? ";
        $ask = <STDIN>;
    }
    if ($ask =~ /n/) {
        print "Program is stopping without execution!\n";
        exit;
    }
    return 1;
}


### Subroutine that counts all bases
sub base_count {
    my $seq = shift;
    $A = ($seq =~ tr/[aA]//) + $A; # pattern match modifier like 'i' don't work with transliterations
    $C = ($seq =~ tr/[cC]//) + $C;
    $G = ($seq =~ tr/[gG]//) + $G;
    $T = ($seq =~ tr/[tT]//) + $T;
    return 1;
}


### Subroutine to concatenate all genbank and fasta files for complete genomes and store in the result folder
sub completes { # subroutine gets one argument $_, the file/directory seen by 'find'
    if (/^\.+$/ || /.*\/Bacteria$/) { # skip system folders (. and ..) and first folder
        return 1;
    } elsif (-d) {
        my $file = $_; # because of 'no_chdir => 1', $_ contains the full path
        $file =~ s/.*\/(\S*)$/$1/; # get the folder name to create the filename
        $file =~ s/Escherichia_coli_/Ecoli_/;
        # $file =~ s/Shigella_/S/; # 'Shigella_D9' ends up as 'SD9'
        dir_check($db); # subroutine to test if result directory exists, if not create
        system("cat $_/*.gbk > ./$db/$file.gbk"); # unix system call for concatenation
        system("cat $_/*.fna > ./$db/$file.fasta");
    }
    return 1;
}


### Count the genes in genbank files
sub count_gene {
    my $file = shift;
    my $gene_count = 0;
    open(IN, "<./$db/$file") or die "Can't open file \'$file\': $!\n";
    while (<IN>) {
        $gene_count++ if /\s{3,}gene\s{3,}/; # Some Genbank files only annotated with 'gene' tags (instead of 'gene' and 'CDS'); at least three whitespaces
    }
    close IN;
    return $gene_count;
}


### Subroutine to check for result directories (./genbank and ./refseq), if they don't exist create; ALWAYS as the prior directory with its contents is deleted at the start of the script with 'rm_dir'
sub dir_check {
    my $dir = shift;
    if (!-e $db) { # mkdir if not existent
        mkdir $db or die "Couldn't create result directory \'$db\': $!\n";
    }
    return 1;
}


### Subroutine to work through all possible draft files and call next subroutine 'unpack_concat'
sub drafts {
    my @result;
    push(@result, unpack_concat('contig.gbk')); # subroutine to unpack and concat resp. files, returns nothing if argument doesn't fit ($result[x] undefined)
    push(@result, unpack_concat('scaffold.gbk'));
    push(@result, unpack_concat('contig.fna'));
    push(@result, unpack_concat('scaffold.fna'));
    foreach (@result) {
        $draft_files{$_} = 1; # save filename in hash '%draft_files' if array element '$result[x]' defined
    }
    return 1;
}


### Subroutine to remove the prior result directories (./genbank and ./refseq) and all their contents, before new files are written in
sub rm_dir {
    my $directory = shift;
    if (-e $directory) {
        opendir (DIR, $directory) or die "Can't opendir $directory: $!"; # system unix call 'rm -rf' is too dangerous
        while (defined(my $file = readdir(DIR))) {
            if ($file =~ /^\.+$/) { # skip system folders (. and ..)
                next;
            }
            unlink "./$directory/$file";
        }
        close DIR;
        rmdir $directory or die "Can't remove directory \'$directory\': $!\n"; # delete empty directory
        return 1;
    }
    return 0;
}


### Subroutine to remove the unpacked files and clean up
sub rm_files {
    my ($directory, $assembly) = @_;
    opendir (DIR, $directory) or die "Can't opendir \'$directory\' to remove the unpacked files: $!";
    $assembly =~ s/.*(\..*)$/$1/; # get the file-extension
    while (defined(my $file = readdir(DIR))) {
        if ($file =~ /.*$assembly$/) { # Why does it not work to include: '-f $file &&'?
            unlink "$directory/$file";
        }
    }
    close DIR;
    return 1;
}


### Calculate the total sequence in multi-genbanks and -fastas
sub seq_size {
    my $file = shift;
    my $seq_count = 0;
    $A = $C = $G = $T = 0;
    open(IN, "<./$db/$file") or die "Can't open file \'$file\': $!\n";
    if ($file =~ /.*\.gbk/) {
        while (my $seq = <IN>) {
            if ($seq =~ /^ORIGIN\s*$/) {
                $seq = <IN>; # Don't want the ORIGIN line, but the next one onwards until '//'
                while ($seq !~ /\/\//) {
                    base_count($seq);
                    $seq = <IN>;
                }
            }
        }
    } elsif ($file =~ /.*\.fasta$/) {
        while (my $seq = <IN>) {
            if ($seq =~ /^>/) {
                next;
            }
            base_count($seq);
        }
    }
    close IN;
    $seq_count = $A + $C + $G + $T;
    return $seq_count;
}


### Subroutine to unzip, unpack all draft files with *.tgz, and concatenate the contigs/scaffolds in the result folder
sub unpack_concat {
    my $assembly = shift;
    if (-f && /.*\.$assembly\.tgz/) { # because of 'no_chdir => 1', $_ contains the full path to the file
        my $file = $_;
        my $directory = $File::Find::dir; # internal variable for File::Find, contains the directory path to the file
        system("tar xfz $file --directory=$directory"); # system call to unpack the $file.tgz
        $file = $File::Find::dir;
        $file =~ s/.*\/(\S*)$/$1/; # get the folder name to create the filename
        $file =~ s/Escherichia_coli_/Ecoli_/;
        # $file =~ s/Shigella_/S/; # 'Shigella_D9' ends up as 'SD9'
        dir_check($db); # subroutine to test if result directory exists, if not create; ALWAYS since above 'rm_dir'
        if ($assembly =~ /contig.gbk/) {
            $file = $file . '_draft_con.gbk';
            system("cat $directory/*.gbk > ./$db/$file");
            rm_files($directory, $assembly); # subroutine to remove unpacked files
            return $file;
        } elsif ($assembly =~ /scaffold.gbk/) {
            $file = $file . '_draft_scaf.gbk';
            system("cat $directory/*.gbk > ./$db/$file");
            rm_files($directory, $assembly);
            return $file;
        } elsif ($assembly =~ /contig.fna/) {
            $file = $file . '_draft_con.fasta';
            system("cat $directory/*.fna > ./$db/$file");
            rm_files($directory, $assembly);
            return $file;
        } elsif ($assembly =~ /scaffold.fna/) {
            $file = $file . '_draft_scaf.fasta';
            system("cat $directory/*.fna > ./$db/$file");
            rm_files($directory, $assembly);
            return $file;
        }
    }
    return; # return undefined to scan @results for the concatenated file in subroutine 'drafts'
}


### Subroutine to warn if the genbank and corresponding fasta file have differing sequence
sub warn_seq {
    my ($file, $gbk_size, $fasta_size) = @_;
    if ($gbk_size != $fasta_size) {
        $file =~ s/\.fasta//; # get rid of file-extension
        print "There is a difference in the sequence length of corresponding files \’$file\.gbk: $gbk_size\’ and \’$file\.fasta: $fasta_size\’!\n"; # ALSO write to an error file???
    }
    return 1;
}
