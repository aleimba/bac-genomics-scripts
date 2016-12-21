#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<trunc_seq.pl> - truncate sequence files

=head1 SYNOPSIS

C<perl trunc_seq.pl 20 3500 seq-file.embl E<gt>
seq-file_trunc_20_3500.embl>

B<or>

C<perl trunc_seq.pl file_of_filenames_and_coords.tsv>

=head1 DESCRIPTION

This script truncates sequence files according to the given
coordinates. The features/annotations in RichSeq files (e.g. EMBL or
GENBANK format) will also be adapted accordingly. Use option B<-o> to
specify a different output sequence format. Input can be given directly
as a file and truncation coordinates to the script, with the start
position as the first argument, stop as the second and (the path to)
the sequence file as the third. In this case the truncated sequence
entry is printed to C<STDOUT>. Input sequence files should contain only
one sequence entry, if a multi-sequence file is used as input only the
B<first> sequence entry is truncated.

Alternatively, a file of filenames (fof) with respective coordinates
and sequence files in the following B<tab-separated> format can be
given to the script (the header is optional):

 #start\tstop\tseq-file
 300\t9000\t(path/to/)seq-file
 50\t1300\t(path/to/)seq-file2

With a fof the resulting truncated sequence files are printed into a
results directory. Use option B<-r> to specify a different results
directory than the default.

It is also possible to truncate a RichSeq sequence file loaded into the
L<Artemis|http://www.sanger.ac.uk/science/tools/artemis> genome browser
from the Sanger Institute: Select a subsequence and then go to Edit
-E<gt> Subsequence (and Features)

=head1 OPTIONS

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-o>=I<str>, B<-outformat>=I<str>

Specify different sequence format for the output (files) [fasta, embl,
or gbk]

=item B<-r>=I<str>, B<-result_dir>=I<str>

Path to result folder for fof input [default = './trunc_seq_results']

=item B<-v>, B<-version>

Print version number to C<STDOUT>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

If a single sequence file is given to the script the truncated sequence
file is printed to C<STDOUT>. Redirect or pipe into another tool as
needed.

=back

B<or>

=over 20

=item F<./trunc_seq_results>

If a fof is given to the script, all output files are stored in a
results folder

=item F<./trunc_seq_results/seq-file_trunc_start_stop.format>

Truncated output sequence files are named appended with 'trunc' and the
corresponding start and stop positions

=back

=head1 EXAMPLES

=over

=item C<perl trunc_seq.pl -o gbk 120 30000 seq-file.embl E<gt>
seq-file_trunc_120_3000.gbk>

=back

B<or>

=over

=item C<perl trunc_seq.pl -o fasta 5300 18500 seq-file.gbk | perl
revcom_seq.pl -i fasta E<gt> seq-file_trunc_revcom.fasta>

=back

B<or>

=over

=item C<perl trunc_seq.pl -r path/to/trunc_embl_dir -o embl
file_of_filenames_and_coords.tsv>

=back

=head1 DEPENDENCIES

=over

=item B<L<BioPerl|http://www.bioperl.org>>

Tested with BioPerl version 1.007001

=back

=head1 VERSION

 0.2                                               update: 2015-12-07
 0.1                                                       2013-08-02

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 (GPLv3) of the
License, or (at your option) any later version.

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
use autodie;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO; # bioperl module to handle sequence input/output
#use Bio::Seq; # bioperl module to handle sequences with features ### apparently not needed, methods inherited
#use Bio::SeqUtils; # bioperl module with additional methods (including features) for Bio::Seq objects ### apparently not needed, methods inherited

### Get options with Getopt::Long
my $Out_Format_Opt; # optional different output seq file format
my $Result_Dir = 'trunc_seq_results'; # path to result folder for fof input
my $VERSION = 0.2;
my ($Opt_Version, $Opt_Help);
GetOptions ('outformat=s' => \$Out_Format_Opt,
            'result_dir=s' => \$Result_Dir,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
if ($Opt_Version) {
    print "$0 $VERSION\n";
    exit;
}



### Check input (@ARGV); didn't include STDIN as input option, too complicated here with fof etc.
my $Fof; # file of filenames (fof) with truncation coords
my $Start;
my $Stop;
my $Seq_File;
if (@ARGV < 1 || @ARGV == 2 || @ARGV > 3) {
    my $warning = "\n### Fatal error: Give either three arguments,\n$0\tstart\tstop\tseq-file\nor one file of sequence filenames with truncation coordinates as argument! Please see the usage with option '-h' if unclear!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
} elsif (@ARGV == 1) { # fof
    check_file_exists($ARGV[0]); # subroutine to check for file existence
    $Fof = shift;
} elsif (@ARGV == 3) {
    check_file_exists($ARGV[2]); # subroutine
    if ($ARGV[0] !~ /^\d+$/ || $ARGV[1] !~ /^\d+$/) {
        my $warning = "\n### Fatal error: With a single sequence file input the first and second arguments are the start and stop positions for truncation, and need to include ONLY digits:\n$0\tstart\tstop\tseq-file\nPlease see the usage with option '-h' if unclear!\n";
        pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
    }
    ($Start, $Stop, $Seq_File) = @ARGV;
}



### Truncate the sequence and write either to STDOUT for single seq file input or output files for fof
if ($Fof) {
    open (my $fof_fh, "<", "$Fof");

    # create result folder
    $Result_Dir =~ s/\/$//; # get rid of a potential '/' at the end of $Result_Dir path
    if (-e $Result_Dir) {
        empty_dir($Result_Dir); # subroutine to empty a directory with user interaction
    } else {
        mkdir $Result_Dir;
    }

    while (my $line = <$fof_fh>) {
        chomp $line;
        next if ($line =~ /^\s*$/ || $line =~ /^#/); # skip empty or comment lines

        die "\n### Fatal error: Line '$.' of the '$Fof' file of sequence filenames plus truncation coordinates does not include the mandatory tab-separated two NUMERICAL start and stop truncation positions, and the sequence file (without any other whitespaces):\nstart\tstop\tpath/to/seq-file\n" if ($line !~ /^\d+\t\d+\t\S+$/);
        ($Start, $Stop, $Seq_File) = split(/\t/, $line);
        check_file_exists($Seq_File); # subroutine

        my ($seqin, $truncseq) = trunc_seq($Start, $Stop, $Seq_File); # subroutine to create a Bio::SeqIO input object and truncate the respective Bio::Seq object
        my $seqout = seq_out($seqin, $Start, $Stop, $Seq_File); # subroutine to create a Bio::SeqIO output object, $seqin needed for format guessing, $Start/$Stop/$Seq_File needed for output filenames
        $seqout->write_seq($truncseq);
    }
    close $fof_fh;

} else { # single seq file, @ARGV == 3
    my ($seqin, $truncseq) = trunc_seq($Start, $Stop, $Seq_File); # subroutine
    my $seqout = seq_out($seqin); # subroutine, without $Start/$Stop/$Seq_file for STDOUT output
    $seqout->write_seq($truncseq);
}

exit;



###############
# Subroutines #
###############

### Subroutine to check if file exists
sub check_file_exists {
    my $file = shift;
    die "\n### Fatal error: File '$file' does not exist: $!\n" if (!-e $file);
}



### Subroutine to empty a directory with user interaction
sub empty_dir {
    my $dir = shift;
    print STDERR "\nDirectory '$dir' already exists! You can use either option '-r' to set a different output result directory name, or do you want to replace the directory and all its contents [y|n]? ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
        unlink glob "$dir/*"; # remove all files in results directory
    } else {
        die "\nScript abborted!\n";
    }
    return 1;
}



### Subroutine to create a Bio::SeqIO output object
sub seq_out {
    my ($seqin, $start, $stop, $seq_file) = @_;

    my $out_format; # need to keep $Out_Format_Opt for several seq files with fof
    if ($Out_Format_Opt) {
        $Out_Format_Opt = 'genbank' if ($Out_Format_Opt =~ /(gbk|gb)/i); # allow shorter input for GENBANK format
        $out_format = $Out_Format_Opt;
    } else { # same format as input file
        if (ref($seqin) =~ /Bio::SeqIO::(genbank|embl|fasta)/) { # from bioperl guessing
            $out_format = $1;
        } else {
            die "\n### Fatal error: Could not determine input file format, please set an output file format with option '-o'!\n";
        }
    }

    my $seqout; # Bio::SeqIO object
    if ($seq_file) { # fof
        $seq_file =~ s/\S+(\/|\\)//; # remove input filepaths, aka 'basename' ('/' for Unix and '\' for Windows)
        my $file_ext;
        if ($out_format eq 'genbank') {
            $file_ext = 'gbk'; # back to shorter file extension for GENBANK format
        } else {
            $file_ext = $out_format;
        }
        $seq_file =~ s/^(\S+)\.\w+$/$Result_Dir\/$1\_trunc_$start\_$stop\.$file_ext/; # append also result directory to output filename
        $seqout = Bio::SeqIO->new(-file => ">$seq_file", -format => $out_format);

    } else { # single seq file input
        $seqout = Bio::SeqIO->new(-fh => \*STDOUT, -format => $out_format); # printing to STDOUT requires '-format'
    }

    return $seqout;
}



### Subroutine create a Bio::SeqIO input object and truncate the respective Bio::Seq object
sub trunc_seq {
    my ($start, $stop, $seq_file) = @_;
    print STDERR "\nTruncating \"$seq_file\" to coordinates $start..$stop ...\n";
    my $seqin = Bio::SeqIO->new(-file => "<$seq_file"); # Bio::SeqIO object; no '-format' given, leave it to bioperl guessing
    my $count = 0;
    my $truncseq;
    while (my $seq_obj = $seqin->next_seq) { # Bio::Seq object
        $count++;
        if ($count > 1) {
            warn "\n### Warning: More than one sequence entry in sequence file '$seq_file', but only the FIRST sequence entry will be truncated and printed to STDOUT or a result file!\n\n";
            last;
        }
        $truncseq = Bio::SeqUtils->trunc_with_features($seq_obj, $start, $stop);
    }
    return ($seqin, $truncseq); # $seqin needed for outformat guessing in subroutine seqout
}
