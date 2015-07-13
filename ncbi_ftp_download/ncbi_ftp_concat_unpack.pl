#!/usr/bin/perl
use warnings;
use strict;
use File::Find; # module to traverse directory trees

my $usage = << "USAGE";

  #######################################################################
  # $0 genbank|refseq [y]                        #
  #                                                                     #
  # Unpacks and concatenates all draft and complete genomes (in genbank #
  # and fasta format) downloaded from NCBI's FTP server                 #
  # (ftp://ftp.ncbi.nlm.nih.gov/). The script traverses the downloaded  #
  # NCBI FTP folder structure and thus has to be called from the top    #
  # level (containing the folder './ftp.ncbi.nlm.nih.gov').             #
  # Therefore, use the bash-shell wrapper script 'ncbi_ftp_download.sh',#
  # which employs 'wget' to download the genomes and mirrors the NCBI   #
  # FTP server folder structure locally (with 'ftp.ncbi.nlm.nih.gov'    #
  # being the top folder). Afterwards, 'ncbi_ftp_download.sh' runs      #
  # 'ncbi_ftp_concat_unpack.pl' with both 'genbank' and 'refseq'        #
  # options, as well as option 'y' to overwrite the old result folders  #
  # (see below). Both scripts have to be in the same directory (or in   #
  # the path) to run 'ncbi_ftp_download.sh'.                            #
  # For COMPLETE genomes PLASMIDS are concatenated to the CHROMOSOMES to#
  # create multi-genbank/-fasta files (script 'split_multi-seq_file.pl' #
  # can be used to split the multi-seq file to single-seq files).       #
  # In DRAFT genomes, SCAFFOLD and/or CONTIG files, designated by       #
  # 'draft_scaf' or 'draft_con', are controlled for annotation (i.e. if #
  # gene primary feature tags exist). Usually only one file contains    #
  # annotations. The one with annotation is then used to create multi-  #
  # genbank files. Multi-fasta files are created for the corresponding  #
  # genbank-file or, if no annotation exists, for the file which        #
  # contains more sequence information (either contigs or scaffolds),   #
  # and if equal, scaffolds are preferred.                              #
  #                                                                     #
  # Use option 'genbank' to extract/copy GenBank genomes and 'refseq'   #
  # for RefSeq genomes. Concatenated files are stored in the result     #
  # folders './genbank' and './refseq', respectively. Set option 'y' in #
  # the program call to delete previous result folders and create new   #
  # ones (otherwise, the script will ask user if to proceed).           #
  # If sequence size discrepancies between a genbank and its            #
  # corresponding fasta file are found, error file 'seq_errors.txt' will#
  # contain warnings.                                                   #
  #                                                                     #
  # version 0.2.1, update: 13.07.2015                  Andreas Leimbach #
  # 15.09.2012                                    aleimba[at]gmx[dot]de #
  #######################################################################

USAGE
;


### Print usage if -h|--h|--help is given as argument or options are not given
my ($db, $ask) = @ARGV; # if $ask not given, run subroutine 'ask_exit' (below) to get input from STDIN
if (!defined $db) {
    die $usage;
} elsif ($db =~ m/-h/) {
    die $usage;
}
my $err = 'seq_errors.txt'; # Error file will be created if sequence discrepancies between genbank and corresponding fasta files will be found (see sub 'warn_seq')
if (-e $err) { # remove error file from previous run, as it will be opened/created in append mode
    unlink $err or die "Can't remove previous error file \'$err\' for new script run: $!\n";
}


### Set the correct directories to start traversion of the folder structure
my $dir_complete; # for complete genomes
my $dir_draft; # for draft genomes
if ($db =~ /genbank/i) {
    $db = 'genbank';
    $dir_complete = './ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria';
    $dir_draft = './ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT';
    ask_exit($db, $ask); # subroutine to ask if the current result folder should be deleted and created new, or exit the script
    rm_dir($db); # subroutine to remove the result directory and all its contents prior filling it with new files (in case files have changed)
} elsif ($db =~ /refseq/i) {
    $db = 'refseq';
    $dir_complete = './ftp.ncbi.nlm.nih.gov/genomes/Bacteria';
    $dir_draft = './ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT';
    ask_exit($db, $ask);
    rm_dir($db);
} else {
    die "\nFatal error: wrong database name, give either 'genbank' or 'refseq' as first argument!\n\n";
}


### Concatenate all genbank and fasta files for COMPLETE genomes in subdirectories and write to the genbank/refseq folders
find ({wanted => \&completes, no_chdir => 1}, $dir_complete); # the function 'find' from File::Find takes a reference to a subroutine and a starting directory to walk through recursively (e.g. 'find (\$concat, $dirname)')
# For each file/directory/link in the tree (also in the current directory) it calls the referenced function and changes to that directory with 'chdir()' --> this can be stopped with 'no_chdir => 1', parameters to find are passed as a hash reference


### Unzip and unpack all DRAFT files with *.tgz in subdirectories and concat in result folders
my %draft_files; # store the filenames of the concatenated files in here (in subroutine 'drafts')
find ({wanted => \&drafts, no_chdir => 1}, $dir_draft);


### For DRAFTS go through all the files in the result folder and keep only gbk files with annotation (either only in contig or scaffold) and the corresponding fastas for each strain; if no annotation exists keep the fasta file with the most sequence information (contig or scaffold); delete the other files (gbk and corresponding fasta) in the result folder
my ($A, $C, $G, $T); # to calculate the total base count (subs 'seq_size' and 'base_count')
my $skip = ''; # not possible to delete elements from a hash while iterating through it, thus need a variable to skip files
foreach my $file (sort{$b cmp $a} keys %draft_files) { # reverse sort to get scaf files first if existent
    $skip =~ s/(\S+)\_draft\_(scaf|con)\.\w*$/$1/; # get rid of file-extension to skip for scaffold genomes corresponding contig gbk
    if ($file =~ /$skip/) { # $skip stores $file from previous foreach
        next;
    }
    if ($file =~ /\_draft\_scaf\.gbk$/) { # if scaffold file exists
        print "\nProcessing DRAFT genome: $file\n"; # status message
        my $scaf_gbkfile = my $scaf_fasfile = $file;
        $scaf_fasfile =~ s/gbk$/fasta/;
        my $con_gbkfile = $scaf_gbkfile;
        $con_gbkfile =~ s/\_scaf\./\_con\./;
        my $con_fasfile = $con_gbkfile;
        $con_fasfile =~ s/gbk$/fasta/;
        my $scaf_size = -s "./$db/$scaf_gbkfile"; # file size of the scaffold gbk file
        my $con_size;
        if (-e "./$db/$con_gbkfile") { # only a scaffold file might exist and not a contig file
            $con_size = -s "./$db/$con_gbkfile";
        } else {
            $con_size = 0;
        }
        my $scaf_gene = count_gene($scaf_gbkfile); # subroutine to count genes in genbank files
        # print "scaf file size: $scaf_size\tscaf_gene: $scaf_gene\n"; # print statement to test output
        my $con_gene = count_gene($con_gbkfile); # 0 if '$con_gbkfile' doesn't exist
        # print "con file size: $con_size\tcon_gene: $con_gene\n"; # print statement to test output
        if ($scaf_gene > $con_gene && $scaf_size > $con_size && $con_gene == 0) { # so many conditions needed? Do they always work?
            print "Annotation in scaffolds: Keeping \'$scaf_gbkfile\' and \'$scaf_fasfile\' and removing contig files (if they exist)!\n"; # status message
            rm_exist_file($db, $con_gbkfile); # remove the contig-gbk without annotation; subroutine, if only a scaffold file exists and not a contig file
            $skip = $file; # scaffold file with annotation found, skip a corresponding contig gbk in %draft_files (if existent)
            my $scaf_gbkcount = seq_size($scaf_gbkfile); # subroutine to calculate total bases (A, C, G, T) of multi-genbank or -fasta files
            rm_exist_file($db, $con_fasfile); # also delete the corresponding fasta
            my $scaf_fascount = seq_size($scaf_fasfile); # count total bases of fasta
            # print "scaf_gbkcount: $scaf_gbkcount\tscaf_fascount: $scaf_fascount\n"; # print statement to test output
            warn_seq($scaf_fasfile, $scaf_gbkcount, $scaf_fascount); # subroutine to warn if the genbank and the fasta file have differing sequence sizes
            next;
        } elsif ($con_gene > $scaf_gene && $con_size > $scaf_size && $scaf_gene == 0) {
            print "Annotation in contigs: Keeping \'$con_gbkfile\' and \'$con_fasfile\' and removing scaffold files!\n";
            unlink "./$db/$scaf_gbkfile" or die "Can't remove file \'$scaf_gbkfile\': $!\n";
            $skip = $file;
            my $con_gbkcount = seq_size($con_gbkfile);
            unlink "./$db/$scaf_fasfile" or die "Can't remove file \'$scaf_fasfile\': $!\n";
            my $con_fascount = seq_size($con_fasfile);
            # print "con_gbkcount: $con_gbkcount\tcon_fascount: $con_fascount\n"; # print statement to test output
            warn_seq($con_fasfile, $con_gbkcount, $con_fascount);
            next;
        } elsif ($con_gene == 0 && $scaf_gene == 0) { # no annotation in both genbanks, just keep one fasta
            print "No annotation in either scaffold or contig file, but ";
            unlink "./$db/$scaf_gbkfile" or die "Can't remove file \'$scaf_gbkfile\': $!\n";
            rm_exist_file($db, $con_gbkfile);
            $skip = $file;
            my $scaf_fascount = seq_size($scaf_fasfile);
            my $con_fascount = seq_size($con_fasfile); # 0 if '$con_gbkfile' doesn't exist
            # print "\nscaf_fascount: $scaf_fascount\tcon_fascount: $con_fascount\n"; # print statement to test output
            if ($scaf_fascount > $con_fascount || $scaf_fascount == $con_fascount) { # keep scaffold file if both same base count
                print "more/equal sequence information in scaffolds (or contigs don't exist). Thus keeping \'$scaf_fasfile\' and removing contig fasta file!\n";
                rm_exist_file($db, $con_fasfile);
            } elsif ($con_fascount > $scaf_fascount) {
                print "more sequence information in contigs. Thus keeping \'$con_fasfile\' and removing scaffold fasta file!\n";
                unlink "./$db/$scaf_fasfile" or die "Can't remove file \'$scaf_fasfile\': $!\n";
            }
            next;
        }
    } elsif ($file =~ /\_draft\_con\.gbk$/) { # if no scaffold file exists, contig file exists
        print "\nProcessing DRAFT genome: $file\n"; # status message
        print "Only contigs exist ";
        $skip = $file;
        my $con_gene = count_gene($file);
        my $con_fasfile = $file;
        $con_fasfile =~ s/gbk$/fasta/;
        if ($con_gene == 0) { # no annotation exists, only keep the fasta
            unlink "./$db/$file" or die "Can't remove file \'$file\': $!\n";
            print "but not annotated, thus keeping only fasta file \'$con_fasfile\'!\n";
        } else {
            print "and annotated, thus keeping contig genbank, \'$file\', and fasta files, \'$con_fasfile\'!\n";
            my $con_gbkcount = seq_size($file);
            my $con_fascount = seq_size($con_fasfile);
            # print "con_gbkcount: $con_gbkcount\tcon_fascount: $con_fascount\n"; # print statement to test output
            warn_seq($con_fasfile, $con_gbkcount, $con_fascount);
        }
    }
}


### Print result folder
print "\n#### All genome files have been written to the folder \'./$db\'!\n";
if (-e $err) {
    print "#### Sequence discrepancies between genbank and corresponding fasta files have been found, see file \'$err\'!\n";
}

exit;


###############
# Subroutines #
###############

### Subroutine to ask if the script should proceed with deleting the previous result folder or exit, if $ask not given as ARGV
sub ask_exit {
    my ($db, $ask) = @_;
    if (!defined($ask) && -e $db) {
        print "The script will delete an already existing result folder \'./$db\', proceed [y|n]? ";
        $ask = <STDIN>;
    } elsif (!defined($ask) && !-e $db) {
        $ask = 'y';
    }
    if ($ask !~ /y/i) {
        die "Script aborted!\n";
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


### Subroutine to concatenate all genbank and fasta files for COMPLETE genomes and store in the result folder
sub completes { # subroutine gets one argument $_, the file/directory seen by 'find'
    if (/^\.+$/ || /.*\/Bacteria$/) { # skip system folders (. and ..) and 'Bacteria' folder
        return 1;
    } elsif (-d) { # '-d' = only directories
        my $file = $_; # because of 'no_chdir => 1', $_ contains the full path to the file/dir
        $file =~ s/.*\/(\S*)$/$1/; # get the folder name to create the filename
        print "\nProcessing COMPLETE genome: $file\n"; # status message
        $file =~ s/Escherichia_coli_/Ecoli_/; # shorten final E. coli file names
        # $file =~ s/Shigella_/S/; # 'Shigella_D9' ends up as 'SD9'
        my $gbk_file = "$file.gbk";
        my $fas_file = "$file.fasta";
        dir_check($db); # subroutine to test if result directory exists, if not create; ALWAYS for first file since above 'rm_dir'
        system("cat $_/*.gbk > ./$db/$gbk_file"); # unix system call for concatenation
        system("cat $_/*.fna > ./$db/$fas_file");
        my $gbk_count = seq_size($gbk_file); # subroutine to calculate total bases (A, C, G, and T) of genbank or fasta files
        my $fas_count = seq_size($fas_file);
        # print "gbk_count: $gbk_count\tfas_count: $fas_count\n"; # print statement to test output
        warn_seq($fas_file, $gbk_count, $fas_count); # subroutine to warn if the genbank and the fasta file have differing sequence sizes
    }
    return 1;
}


### Count the genes in genbank files
sub count_gene {
    my $file = shift;
    my $gene_count = 0;
    if (-e "./$db/$file") { # if the file doesn't exist (e.g. only scaffold file and not contig file) $gene_count will remain 0
        open(IN, "<./$db/$file") or die "Can't open file \'$file\': $!\n";
        while (<IN>) {
            $gene_count++ if /\s{3,}gene\s{3,}/; # Some Genbank files only annotated with 'gene' tags (instead of 'gene' and 'CDS'); '\s{3,}' at least three whitespaces
        }
        close IN;
    }
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


### Subroutine to work through all possible DRAFT files and call next subroutine 'unpack_concat'
sub drafts {
    my @result; # stores results of sub 'unpack_concat'
    push(@result, unpack_concat('contig.gbk')); # subroutine to unpack and concat respective files; returns nothing if argument (here 'contig.gbk') doesn't fit to filename and $result[x] will be undefined
    push(@result, unpack_concat('scaffold.gbk'));
    push(@result, unpack_concat('contig.fna'));
    push(@result, unpack_concat('scaffold.fna'));
    foreach (@result) { # for each succesful subroutine run @result contains one element, i.e. the concatenated filename
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
            unlink "./$directory/$file" or die "Can't remove file \'$file\': $!\n";
        }
        close DIR;
        rmdir $directory or die "Can't remove directory \'$directory\': $!\n"; # delete empty directory
        return 1;
    }
    return 0;
}


### Test for file existence and remove file
sub rm_exist_file {
    my ($db, $file) = @_;
    if (-e "./$db/$file") {
        unlink "./$db/$file" or die "Can't remove file \'$file\': $!\n";
    }
    return 1;
}


### Subroutine to remove the unpacked files and clean up
sub rm_files {
    my ($directory, $assembly) = @_;
    opendir (DIR, $directory) or die "Can't opendir \'$directory\' to remove the unpacked files: $!";
    $assembly =~ s/.*(\..*)$/$1/; # get the file-extension to remove files
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
    if (-e "./$db/$file") { # if the file doesn't exist (e.g. only scaffold file and not contig file) $seq_count will remain 0
        $A = $C = $G = $T = 0;
        open(IN, "<./$db/$file") or die "Can't open file \'$file\': $!\n";
        if ($file =~ /.*\.gbk/) {
            while (my $seq = <IN>) {
                if ($seq =~ /^ORIGIN\s*$/) {
                    $seq = <IN>; # don't want the ORIGIN line, but the next one onwards until '//'
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
    }
    return $seq_count;
}


### Subroutine to unzip, unpack all DRAFT files with *.tgz, and concatenate the contigs/scaffolds in the result folder
sub unpack_concat {
    my $assembly = shift;
    if (-f && /.*\.$assembly\.tgz/) { # only if $_ is a file (-f) and fits to $assembly.tgz; because of 'no_chdir => 1' $_ contains the full path to the file (with filename)
        my $file = $_;
        my $directory = $File::Find::dir; # internal variable for File::Find, contains the directory path to the file
        system("tar xfz $file --directory=$directory"); # system call to unpack the $file.tgz
        $file = $File::Find::dir; # replace path+file to just directory path
        $file =~ s/.*\/(\S*)$/$1/; # get the folder name to create the filename
        $file =~ s/Escherichia_coli_/Ecoli_/; # shorten final E. coli file names
        # $file =~ s/Shigella_/S/; # 'Shigella_D9' ends up as 'SD9'
        dir_check($db); # subroutine to test if result directory exists, if not create; ALWAYS for first file since above 'rm_dir'
        if ($assembly =~ /contig.gbk/) {
            $file = $file . '_draft_con.gbk';
            system("cat $directory/*.gbk > ./$db/$file");
            rm_files($directory, $assembly); # subroutine to remove unpacked files
            return $file; # @results in sub 'drafts' defined
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
        open (ERR, ">>$err");
        print ERR "There is a difference in the sequence length of corresponding files \'$file\.gbk: $gbk_size\' and \'$file\.fasta: $fasta_size\'!\n\n";
        close ERR;
    }
    return 1;
}
