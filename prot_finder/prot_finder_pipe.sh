#!/bin/bash
set -e

#############
# Functions #
#############

usage () {
    cat 1>&2 << EOF # ${0##*/} parameter expansion substitution with variable '0' to get shell script filename without path
Usage: ${0##*/} [OPTION] -q query.faa -f (embl|gbk) > blast_hits.tsv
or:    ${0##*/} [OPTION] -q query.faa -s subject.faa -d result_dir \\
       > result_dir/blast_hits.tsv

Bash wrapper script to run a pipeline consisting of optional
'cds_extractor.pl' (with its options '-p -f'), BLASTP, 'prot_finder.pl',
and optional Clustal Omega. 'cds_extractor.pl' (only for shell script
option '-f') and 'prot_finder.pl' either have to be installed in the
global PATH or present in the current working directory. BLASTP is run
with disabled query filtering, locally optimal Smith-Waterman alignments,
and increasing the number of database sequences to show alignments
to 500 for BioPerl parsing (legacy: '-F F -s T -b 500', plus: '-seg
no -use_sw_tback -num_alignments 500').

The script ends with the STDERR message 'Pipeline finished!', if this
is not the case have a look at the log files in the result directory
for errors.

Mandatory options:
    -q <str>           Path to query protein multi-FASTA file (*.faa)
                       with unique FASTA IDs
    -f <str>           File extension for files in the current working
                       directory to use for 'cds_extractor.pl' (e.g.
                       'embl' or 'gbk'); excludes shell script option '-s'
    or
    -s <str>           Path to subject protein multi-FASTA file (*.faa)
                       already created with 'cds_extractor.pl' (and its
                       options '-p -f'), will not run 'cds_extractor.pl';
                       excludes shell script option '-f'

Optional options:
    -h                 Print usage
    -d <str>           Path to result folder [default = results_i#_cq#]
    -p (legacy|plus)   BLASTP suite to use [default = plus]
    -e <real>          E-value for BLASTP [default = 1e-10]
    -t <int>           Number of threads to be used for BLASTP and
                       Clustal Omega [default = all processors on
                       system]
    -i <int>           Query identity cutoff for significant hits
                       [default = 70]
    -c <int>           Query coverage cutoff [default = 70]
    -k <int>           Subject coverage cutoff [default = 0]
    -b                 Give only best hit (highest identity) for each
                       subject sequence
    -a                 Multiple alignment of each multi-FASTA result
                       file with Clustal Omega
    -o <str>           Path to executable Clustal Omega binary if not
                       in global PATH; requires shell script option '-a'
    -m                 Clean up all non-essential files

Author: Andreas Leimbach <aleimba[at]gmx[dot]de>
EOF
}


### Check external dependencies
check_commands () {
    which "$1" > /dev/null || err "Required executable '$1' not found in global PATH, please install.$2"
}

### Check cutoff options input
check_cutoff_options () {
    local message="Option '-$2' requires an integer number >= 0 or <= 100 as value, not '$1'!"
    [[ $1 =~ ^[0-9]+$ ]] || err "$message"
    (( $1 <= 100 )) || err "$message" # arithmetic expression (can only handle integer math, not float)
}


### Error messages
err () {
    echo -e "\n### Fatal error: $*" 1>&2
    exit 1
}


### Run status of script to STDERR instead of STDOUT
msg () {
    echo -e "# $*" 1>&2
}


########
# MAIN #
########

shopt -s extglob # enable extended globs for bash

Cmdline="$*"

### Getopts
Blastp_Suite="plus"
Evalue="1e-10"
Threads="$(nproc --all)" # get max number of processors on system
Ident_Cut=70
Cov_Query_Cutoff=70
Cov_Subject_Cutoff=0

while getopts ':q:f:s:d:p:e:t:i:c:k:bao:mh' opt; do # beginning ':' indicates silent mode, trailing ':' after each option requires value
    case $opt in
        q) Query_File=$OPTARG
           [[ -r $Query_File ]] || err "Cannot read query file '$Query_File'!"
           ;;
        f) Subject_Ext=$OPTARG
           [[ -n "$(find . -maxdepth 1 -name "*.${Subject_Ext}" -print -quit)" ]] || err "No files with the option '-f' specified file extension '$Subject_Ext' found in the current working directory!"
           ;;
        s) Subject_File=$OPTARG
           [[ -r $Subject_File ]] || err "Cannot read subject file '$Subject_File'!"
           ;;
        d) Result_Dir=$OPTARG;; # checked below
        p) Blastp_Suite=$OPTARG
           [[ $Blastp_Suite = @(plus|legacy) ]] || err "Option '-p' only allows 'plus' for BLASTP+ or 'legacy' for legacy BLASTP as value, not '$Blastp_Suite'!" # extended glob (regex more expensive)
           ;;
        e) Evalue=$OPTARG
           [[ $Evalue =~ ^([0-9][0-9]*|[0-9]+e-[0-9]+)$ ]] || err "Option '-e' requires a real number (either integer or scientific exponential notation) as value, not '$Evalue'!"
           ;;
        t) Threads=$OPTARG
           [[ $Threads =~ ^[1-9][0-9]*$ ]] || err "Option '-t' requires an integer > 0 as value, not '$Threads'!"
           ;;
        i) Ident_Cut=$OPTARG
           check_cutoff_options "$Ident_Cut" "i"
           ;;
        c) Cov_Query_Cutoff=$OPTARG
           check_cutoff_options "$Cov_Query_Cutoff" "c"
           ;;
        k) Cov_Subject_Cutoff=$OPTARG
           check_cutoff_options "$Cov_Subject_Cutoff" "k"
           ;;
        b) Opt_Best_Hit=1;;
        a) Opt_Align=1;;
        o) Clustal_Path=$OPTARG
           [[ -x $Clustal_Path ]] || err "Option '-o' requires the path to an executable Clustal Omega binary as value, not '$Clustal_Path'!"
           ;;
        m) Opt_Clean_Up=1;;
        h) usage; exit;; # usage function, exit code zero
        \?) err "Invalid option '-$OPTARG'. See usage with '-h'!";;
        :) err "Option '-$OPTARG' requires a value. See usage with '-h'!";;
    esac
done


### Check options and enforce mandatory options
[[ $Query_File && ($Subject_Ext || $Subject_File) ]] || err "Mandatory options '-q' and '-f' or '-s' are missing!"

[[ $Subject_Ext && $Subject_File ]] && err "Options '-f' and '-s' given which exclude themselves. Choose either '-f' OR '-s'!"

(( Threads <= $(nproc) )) || err "Number of threads for option '-t', '$Threads', exceeds the maximum $(nproc) processors on the system!"

[[ ! $Opt_Align && $Clustal_Path ]] && Opt_Align=1 && msg "Option '-o' requires option '-a', forcing option '-a'!"


### Check external dependencies
echo 1>&2 # newline
msg "Checking pipeline dependencies"
[[ $Opt_Align && ! $Clustal_Path ]] && check_commands "clustalo" " Or use option '-o' to give the path to the binary!"

for exe in cds_extractor.pl formatdb blastall makeblastdb blastp prot_finder.pl; do
    [[ $Subject_File && $exe == cds_extractor.pl ]] && continue
    [[ $Blastp_Suite == legacy && $exe = @(makeblastdb|blastp) ]] && continue # extended glob
    [[ $Blastp_Suite == plus && $exe = @(formatdb|blastall) ]] && continue
    if [[ $exe = *.pl ]]; then # glob
        if [[ -r "./$exe" ]]; then # present in current wd
            [[ $exe =~ ^cds ]] && Cds_Extractor_Cmd="perl cds_extractor.pl"
            [[ $exe =~ ^prot ]] && Prot_Finder_Cmd="perl prot_finder.pl"
            continue
        else
            [[ $exe =~ ^cds ]] && Cds_Extractor_Cmd="cds_extractor.pl"
            [[ $exe =~ ^prot ]] && Prot_Finder_Cmd="prot_finder.pl"
            check_commands "$exe" " Or copy the Perl script in the current working directory."
        fi
        continue
    fi
    check_commands "$exe"
done

msg "Script call command: ${0##*/} $Cmdline"


### Create result folder
if [[ ! $Result_Dir ]]; then # can't give default before 'getopts' in case cutoffs are set by the user
    Result_Dir="results_i${Ident_Cut}_cq${Cov_Query_Cutoff}"
else
    Result_Dir="${Result_Dir%/}" # parameter expansion substitution to get rid of a potential '/' at the end of Result_Dir path
fi

if [[ -d $Result_Dir ]]; then # make possible to redirect STDOUT output into result_dir (corresponding to option '-f' in 'protein_finder.pl' script)
    skip=0
    for file in "$Result_Dir"/*; do
        if [[ -s $file || $skip -eq 1 ]]; then # die if a file with size > 0 or more than one file already in result_dir
            err "Result directory '$Result_Dir' already exists! You can use option '-d' to set a different result directory name."
        fi
        skip=1
    done
else
    mkdir -pv "$Result_Dir" 1>&2
fi


### Run cds_extractor.pl
if [[ $Subject_Ext ]]; then
    msg "Running cds_extractor.pl on all '*.$Subject_Ext' files in the current working directory"
    for file in *."$Subject_Ext"; do
        file_no_ext="${file%.${Subject_Ext}}.faa" # parameter expansion substitution to get rid of file extension and replace with new one (*.faa are the output files from cds_extractor)
        File_Names+=("$file_no_ext") # append to array
        eval "$Cds_Extractor_Cmd -i $file -p -f &>> $Result_Dir/cds_extractor.log" # '&>' instead of '/dev/null' for error catching
    done
    Subject_File="$Result_Dir/prot_finder.faa" # for creating BLASTP db below
    cat "${File_Names[@]}" > "$Subject_File" # concatenate files stored in the array, "${array[@]}" expands to list of array elements (words)
fi


### Run BLASTP
msg "Running BLASTP '$Blastp_Suite' with subject '$Subject_File', query '$Query_File', evalue '$Evalue', and $Threads threads"
Blast_Report="$Result_Dir/prot_finder.blastp"
if [[ $Blastp_Suite == legacy ]]; then
    formatdb -p T -i "$Subject_File" -n prot_finder_db
    blastall -p blastp -d prot_finder_db -i "$Query_File" -o "$Blast_Report" -e "$Evalue" -F F -s T -b 500 -a "$Threads"
elif [[ $Blastp_Suite == plus ]]; then
    makeblastdb -in "$Subject_File" -input_type fasta -dbtype prot -out prot_finder_db &> "$Result_Dir/makeblastdb.log" # '&>' instead of '/dev/null' for error catching
    blastp -db prot_finder_db -query "$Query_File" -out "$Blast_Report" -evalue "$Evalue" -seg no -use_sw_tback -num_alignments 500 -num_threads "$Threads"
fi


### Run prot_finder.pl
msg "Running prot_finder.pl with identity cutoff '$Ident_Cut', query coverage cutoff '$Cov_Query_Cutoff', and subject coverage cutoff '$Cov_Subject_Cutoff'"
Cmd="$Prot_Finder_Cmd -d $Result_Dir -f -q $Query_File -s $Subject_File -r $Blast_Report -i $Ident_Cut -cov_q $Cov_Query_Cutoff -cov_s $Cov_Subject_Cutoff"
[[ $Opt_Best_Hit ]] && Cmd="$Cmd -b" # append to command
[[ $Opt_Align ]] && Cmd="$Cmd -a -t $Threads"
[[ $Clustal_Path ]] && Cmd="$Cmd -p $Clustal_Path"
eval "$Cmd" 2> "$Result_Dir/prot_finder.log" # '2>' instead of '/dev/null' for error catching

msg "All result files stored in directory '$Result_Dir'"


### Clean up non-essential files
if [[ $Opt_Clean_Up ]]; then
    msg "Removing non-essential output files, option '-m'"
    for file in "${File_Names[@]}"; do # remove output files from cds_extractor
        rm -v "$file" 1>&2
    done
    [[ $Subject_Ext ]] && rm -v "$Subject_File" 1>&2 # 'cat' from cds_extractor
    if [[ $Blastp_Suite == legacy ]]; then
        rm -v formatdb.log 1>&2
        [[ -r error.log ]] && rm -v error.log 1>&2 # no idea where this guy is coming from or what is its trigger
    fi
    rm -v prot_finder_db.p* "$Blast_Report" "$Result_Dir"/*.log "${Subject_File}.idx" 1>&2
fi

msg "Pipeline finished!"
