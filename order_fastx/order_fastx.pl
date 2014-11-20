#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<order_fastx.pl> - order sequences in FASTA or FASTQ files

=head1 SYNOPSIS

C<perl order_fastx.pl -i infile.fasta -l order_id_list.txt
E<gt> ordered.fasta>

=head1 DESCRIPTION

Order sequence entries in FASTA or FASTQ sequence files according to
an ID list with a given order. Beware, the IDs in the order list
have to be B<identical> to the entire IDs in the sequence file.

However, the ">" or "@" ID identifiers of FASTA or FASTQ files,
respectively, can be omitted in the ID list.

The file type is detected automatically. But, you can set the file
type manually with option B<-f>. FASTQ format assumes B<four> lines
per read, if this is not the case run the FASTQ file through
L<C<fastx_fix.pl>|/fastx_fix> or use Heng Li's L<C<seqtk
seq>|https://github.com/lh3/seqtk>:

C<seqtk seq -l 0 infile.fq E<gt> outfile.fq>

The script can also be used to pull a subset of sequences in the ID
list from the sequence file. Probably best to set option flag B<-s>
in this case, see L<"Optional options"> below. But, rather use
L<C<filter_fastx.pl>|/filter_fastx>.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input FASTA or FASTQ file

=item B<-l>=I<str>, B<-list>=I<str>

List with sequence IDs in specified order

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-f>=I<fasta|fastq>, B<-file_type>=I<fasta|fastq>

Set the file type manually

=item B<-e>, B<-error_files>

Write missing IDs in the seq file or the order ID list without an
equivalent in the other to error files instead of C<STDERR> (see
L<"OUTPUT"> below)

=item B<-s>, B<-skip_errors>

Skip missing ID error statements, excludes option B<-e>

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The newly ordered sequences are printed to C<STDOUT>. Redirect or
pipe into another tool as needed.

=item (F<order_ids_missing.txt>)

If IDs in the order list are missing in the sequence file with
option B<-e>

=item (F<seq_ids_missing.txt>)

If IDs in the sequence file are missing in the order ID list with
option B<-e>

=back

=head1 EXAMPLES

=over

=item C<perl order_fastx.pl -i infile.fq -l order_id_list.txt -s -f
fastq E<gt> ordered.fq>

=item C<perl order_fastx.pl -i infile.fasta -l order_id_list.txt -e
E<gt> ordered.fasta>

=back

=head1 VERSION

 0.1                                                       20-11-2014

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

### Get the options with Getopt::Long
my $Seq_File; # sequence file to order sequences in
my $Order_List; # order ID list for seq file
my $File_Type; # set file type; otherwise detect file type by file extension
my $Opt_Error_Files; # print missing IDs not found in order list or seq file to error files instead of STDERR
my $Opt_Skip_Errors; # skip missing IDs error statements/files
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Seq_File,
            'list=s' => \$Order_List,
            'file_type=s' => \$File_Type,
            'error_files' => \$Opt_Error_Files,
            'skip_errors' => \$Opt_Skip_Errors,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Seq_File || !$Order_List) {
    my $warning = "\n### Fatal error: Options '-i' or '-l' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Enforce mandatory or optional options
die "\n### Fatal error:\nUnknown file type '$File_Type' given with option '-f'. Please choose from either 'fasta' or 'fastq'!\n" if ($File_Type && $File_Type !~ /(fasta|fastq)/i);
warn "\n### Warning:\nIgnoring option flag '-e', because option '-s' set at the same time!\n\n" if ($Opt_Error_Files && $Opt_Skip_Errors);



### Order input FASTA/FASTQ file according to a given list
open (my $Order_List_Fh, "<", "$Order_List");
open (my $Input_Fh, "<", "$Seq_File"); # pipe from STDIN not working because of 'seek' on filehandle
get_file_type() if (!$File_Type); # subroutine to determine file type by file extension

my %Order_List_IDs; # store order IDs and indicate if found in seq file
my %Seq_File_IDs; # store seq file IDs and indicate if present in order list

my $Next_Fasta_ID; # for multi-line FASTA input files to store next entry header/ID line while parsing in subroutine 'get_fastx_entry'
my $Parse_Run = 1; # indicate FIRST parsing cycle through seq file to collect all seq IDs

while (my $ord_id = <$Order_List_Fh>) {
    chomp $ord_id;
    next if ($ord_id =~ /^\s*$/); # skip emtpy lines in order list
    $ord_id =~ s/^(>|@)//; # remove ">/@" for WHOLE string regex match -> ID in order list can be given with ">/@" or without (will be appended again in print)

    if ($Order_List_IDs{$ord_id}) {
        die "\n### Fatal error:\n'$ord_id' exists several times in '$Order_List' and IDs should be unique!\n";
    } else {
        $Order_List_IDs{$ord_id} = 1; # changes to 2 if ID was found in seq file
    }

    while (<$Input_Fh>) {
        if (/^\s*$/) { # skip empty lines in input
            warn "\n### Warning:\nFASTQ file includes empty lines, which is unusual. Parsing the FASTQ reads might fail so check the output file afterwards if the script didn't quit with a fatal error. However, consider running the input FASTQ file through 'fix_fastx.pl'!\n\n" if ($File_Type =~ /fastq/i);
            next;
        }
        chomp;

        # FASTA file
        if ($File_Type =~ /fasta/i) {
            $_ = get_fastx_entry($_); # subroutine to read one FASTA sequence entry (seq in multi-line or not), returns anonymous array

        # FASTQ file
        } elsif ($File_Type =~ /fastq/i) {
            $_ = get_fastx_entry($_); # subroutine to read one FASTQ read composed of FOUR mandatory lines, returns reference to array
        }

        if ($Seq_File_IDs{$_->[0]} && $Parse_Run == 1) { # only for first parse cycle, subsequent parsings of course will find the same IDs
            die "\n### Fatal error:\n'$_->[0]' exists several times in '$Seq_File' and IDs should be unique!\n";
        } elsif (!$Seq_File_IDs{$_->[0]}) {
            $Seq_File_IDs{$_->[0]} = 1; # changes to 2 if present in order list
        }

        if ($ord_id =~ /^$_->[0]$/) { # order ID hit in seq file with the WHOLE string; de-reference array
            $Order_List_IDs{$ord_id} = 2; # set to ID found
            $Seq_File_IDs{$_->[0]} = 2;

            print ">$_->[0]\n$_->[1]\n\n" if ($File_Type =~ /fasta/i); # print seq entry to STDOUT
            print "\@$_->[0]\n$_->[1]\n$_->[2]\n$_->[3]\n" if ($File_Type =~ /fastq/i);

            next if ($Parse_Run == 1); # parse the complete seq file once (skip 'last' below) to collect all seq IDs
            last; # jump out of seq file 'while'
        }
    }

    # rewind seq file for next order list ID
    $Next_Fasta_ID = '';
    seek $Input_Fh, 0, 0;
    $. = 0; # set line number of seq file to 0 (seek doesn't do it automatically)
    $Parse_Run = 0;
}
close $Input_Fh;
close $Order_List_Fh;



### Print order and seq IDs that were not found in seq file or in order list, resp.
if (!$Opt_Skip_Errors) {
    # order IDs not found in seq file
    missing_IDs(\%Order_List_IDs, 'order_ids_missing.txt', 'order'); # subroutine to identify and print missing IDs

    # seq file IDs not found in order list
    missing_IDs(\%Seq_File_IDs, 'seq_ids_missing.txt', 'sequence'); # subroutine
}

exit;


#############
#Subroutines#
#############

### Test for output file existence and give warning to STDERR
sub file_exist {
    my $file = shift;
    if (-e $file) {
        warn "\n### Warning:\nThe error file '$file' exists already, the current errors will be appended to the existing file!\n";
        return 1;
    }
    return 0;
}



### Get sequence entries from FASTA/Q file
sub get_fastx_entry {
    my $line = shift;

    # possible multi-line seq in FASTA
    if ($File_Type =~ /fasta/i) {
        my ($seq, $header);
        if ($. == 1) { # first line of file
            die "\n### Fatal error:\nNot a FASTA input file, first line of file should be a FASTA ID/header line and start with a '>':\n$line\n" if ($line !~ /^>/);
            $header = $line;
        } elsif ($Next_Fasta_ID) {
            $header = $Next_Fasta_ID;
            $seq = $line;
        }
        while (<$Input_Fh>) {
            chomp;
            $Next_Fasta_ID = $_ if (/^>/); # store ID/header for next seq entry
            $header =~ s/^>//; # remove '>' for WHOLE string regex match in MAIN
            return [$header, $seq] if (/^>/); # return anonymous array with current header and seq
            $seq .= $_; # concatenate multi-line seq
        }
        $header =~ s/^>//; # see above
        return [$header, $seq] if eof;

    # FASTQ: FOUR lines for each FASTQ read (seq-ID, sequence, qual-ID [optional], qual)
    } elsif ($File_Type =~ /fastq/i) {
        my @fastq_read;

        # read name/ID, line 1
        my $seq_id = $line;
        die "\n### Fatal error:\nThis read doesn't have a sequence identifier/read name according to FASTQ specs, it should begin with a '\@':\n$seq_id\n" if ($seq_id !~ /^@.+/);
        $seq_id =~ s/^@//; # remove '@' to make comparable to $qual_id and for WHOLE string regex match in MAIN
        push(@fastq_read, $seq_id);

        # sequence, line 2
        chomp (my $seq = <$Input_Fh>);
        die "\n### Fatal error:\nRead '$seq_id' has a whitespace in its sequence, which is not allowed according to FASTQ specs:\n$seq\n" if ($seq =~ /\s+/);
        die "\n### Fatal error:\nRead '$seq_id' has a IUPAC degenerate base (except for 'N') or non-nucleotide character in its sequence, which is not allowed according to FASTQ specs:\n$seq\n" if ($seq =~ /[^acgtun]/i);
        push(@fastq_read, $seq);

        # optional quality ID, line 3
        chomp (my $qual_id = <$Input_Fh>);
        die "\n### Fatal error:\nThe optional sequence identifier/read name for the quality line of read '$seq_id' is not according to FASTQ specs, it should begin with a '+':\n$qual_id\n" if ($qual_id !~ /^\+/);
        push(@fastq_read, $qual_id);
        $qual_id =~ s/^\+//; # if optional ID is present check if equal to $seq_id in line 1
        die "\n### Fatal error:\nThe sequence identifier/read name of read '$seq_id' doesn't fit to the optional ID in the quality line:\n$qual_id\n" if ($qual_id && $qual_id ne $seq_id);

        # quality, line 4
        chomp (my $qual = <$Input_Fh>);
        die "\n### Fatal error:\nRead '$seq_id' has a whitespace in its quality values, which is not allowed according to FASTQ specs:\n$qual\n" if ($qual =~ /\s+/);
        die "\n### Fatal error:\nRead '$seq_id' has a non-ASCII character in its quality values, which is not allowed according to FASTQ specs:\n$qual\n" if ($qual =~ /[^[:ascii]]/);
        die "\n### Fatal error:\nThe quality line of read '$seq_id' doesn't have the same number of symbols as letters in the sequence:\n$seq\n$qual\n" if (length $qual != length $seq);
        push(@fastq_read, $qual);

        return \@fastq_read; # return array-ref
    }
    return 0;
}



### Determine file type via file extension (FASTA or FASTQ)
sub get_file_type {
    if ($Seq_File =~ /.+\.(fa|fas|fasta|ffn|fna|frn|fsa)$/) { # use "|fsa)(\.gz)*$" if unzip inside script
        $File_Type = 'fasta';
    } elsif ($Seq_File =~ /.+\.(fastq|fq)$/) {
        $File_Type = 'fastq';
    }

    die "\n### Fatal error:\nFile type could not be automatically detected. Sure this is a FASTA/Q file? If yes, you can force the file type by setting option '-f' to either 'fasta' or 'fastq'!\n" if (!$File_Type);
    print STDERR "Detected file type: $File_Type\n";
    return 1;
}



### Identify and print IDs with no hit in order list or seq file
sub missing_IDs {
    my ($hash_ref, $error_file, $mode) = @_;

    my @missed = grep ($hash_ref->{$_} == 1, keys %$hash_ref); # set to 2 if hit, 1 if "only" present

    if (@missed) {
        file_exist($error_file) if ($Opt_Error_Files); # subroutine
        open (my $error_fh, ">>", $error_file) if ($Opt_Error_Files);

        print STDERR "\n### Warning:\nSome $mode IDs were not found in '";
        if ($mode eq 'order') {
            print STDERR "$Seq_File";
        } elsif ($mode eq 'sequence') {
            print STDERR "$Order_List";
        }
        print STDERR "', listed ";
        print STDERR "below:\n" if (!$Opt_Error_Files);
        print STDERR "in error file '$error_file'!\n" if ($Opt_Error_Files);

        foreach (sort @missed) {
            if (!$Opt_Error_Files) {
                print STDERR "$_\t"; # separated by tab
            } elsif ($Opt_Error_Files) {
                print $error_fh "$_\n"; # separated by newline
            }
        }
        print STDERR "\n" if (!$Opt_Error_Files); # final newline for STDERR print

        close $error_fh if ($Opt_Error_Files);
    }
    return 1;
}
