#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<sample_fastx-txt.pl> - random subsampling of FASTA, FASTQ, or TEXT files

=head1 SYNOPSIS

C<perl sample_fastx-txt.pl -i infile.fasta -n 100 E<gt> subsample.fasta>

B<or>

C<zcat reads.fastq.gz | perl sample_fastx-txt.pl -i - -n 100000
E<gt> subsample.fastq>

=head1 DESCRIPTION

Randomly subsample FASTA, FASTQ, and TEXT files.

Empty lines in the input files will be skipped and not included in
sampling. Format TEXT assumes one entry per single line. FASTQ
format assumes B<four> lines per read, if this is not the case run
the FASTQ file through L<C<fastx_fix.pl>|/fastx_fix> or use Heng
Li's L<C<seqtk seq>|https://github.com/lh3/seqtk>:

C<seqtk seq -l 0 infile.fq E<gt> outfile.fq>

The file type is detected automatically. However, if automatic
detection fails, TEXT format is assumed. As a last resort, you can
set the file type manually with option B<-f>.

This script is an implementation of the I<reservoir sampling>
algorithm (or I<Algorithm R (3.4.2)>) described in Donald Knuth's
L<I<The Art of Computer Programming>|https://en.wikipedia.org/wiki/The_Art_of_Computer_Programming>.
It is designed to randomly pull a small sample size from a
(potential) huge input file of indeterminate size, which
(potentially) doesn't fit into main memory. The beauty of reservoir
sampling is that it requires only one pass through the input file.
The memory consumption of the algorithm is proportional to the
sample size, thus large sample sizes will consume lots of memory as
the whole sample will be held in memory. On the other hand, the size
of the initial file is irrelevant.

An alternative tool, which is a lot faster, is C<seqtk sample> from
the L<I<seqtk toolkit>|https://github.com/lh3/seqtk>.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input FASTA/Q or TEXT file, or piped C<STDIN> (-)

=item B<-n>=I<int>, B<-num>=I<int>

Number of entries/reads to subsample

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-f>=I<fasta|fastq|text>, B<-file_type>=I<fasta|fastq|text>

Set the file type manually

=item B<-s>=I<int>, B<-seed>=I<int>

Set starting random seed. For B<paired-end> read data use the B<same
random seed> for both FASTQ files with option B<-s> to retain
pairing (see L<"EXAMPLES"> below).

=item B<-t>=I<int>, B<-title_skip>=I<int>

Skip the specified number of header lines in TEXT files before
subsampling and append them again afterwards. If you want to get rid
of the header as well, pipe the subsample output to
L<C<sed>|https://www.gnu.org/software/sed/manual/sed.html> (see
C<man sed> and L<"EXAMPLES"> below).

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The subsample of the input file is printed to C<STDOUT>. Redirect or
pipe into another tool as needed.

=back

=head1 EXAMPLES

=over

=item C<perl sample_fastx-txt.pl -i read-pair_1.fq -n 1000000 -s 123
E<gt> sub-pair_1.fq>

=item C<perl sample_fastx-txt.pl -i read-pair_2.fq -n 1000000 -s 123
E<gt> sub-pair_2.fq>

=item C<perl sample_fastx-txt.pl -i infile.txt -n 100 -f text -t 3
E<gt> subsample.txt>

=item C<perl sample_fastx-txt.pl -i infile.txt -n 350 -t 2 | sed
'1,2d' E<gt> sub_no-header.txt>

=back

=head1 VERSION

 0.1                                                       18-11-2014

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

=head1 ACKNOWLEDGEMENTS

I got the idea for reservoir sampling from Sean Eddy's keynote at
the Janelia meeting on L<I<High Throughput Sequencing for
Neuroscience>|http://cryptogenomicon.wordpress.com/2014/11/01/high-throughput-sequencing-for-neuroscience/>
which he posted in his blog
L<I<Cryptogenomicon>|http://cryptogenomicon.wordpress.com/>. The
L<I<Wikipedia
article>|https://en.wikipedia.org/wiki/Reservoir_sampling> and the
L<I<PerlMonks>|http://www.perlmonks.org/index.pl?node_id=177092>
implementation helped a lot, as well.

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

### Get options with Getopt::Long
my $Input_File; # input file to subsample from
my $Sample_Num; # number of sequence entries/reads/lines to sample
my $File_Type; # set file type; otherwise detect file type by file extension, or by looking at the first line of the file if input is piped from STDIN
my $Seed; # give starting seed number for 'rand' to make results repeatable; also needed for paired FASTA/Q files
my $Title_Skip; # optionally skip given number of header lines of TEXT files
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_File,
            'num=i' => \$Sample_Num,
            'file_type=s' => \$File_Type,
            'seed=i' => \$Seed,
            'title_skip=i' => \$Title_Skip,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_File || !$Sample_Num) {
    my $warning = "\n### Fatal error: Options '-i' or '-n' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}
die "\n### Fatal error:\nUnknown file type '$File_Type' given with option '-f'. Please choose from either 'fasta', 'fastq', or 'text'!\n" if ($File_Type && $File_Type !~ /(fasta|fastq|text)/i);



### Pipe input from STDIN or open input file
my $Input_Fh;
if ($Input_File eq '-') { # file input via STDIN
    $Input_Fh = *STDIN; # capture typeglob of STDIN
} else { # input via input file
    open ($Input_Fh, "<", "$Input_File");
    get_file_type('ext') if (!$File_Type); # subroutine to determine file type; mode 'ext' (extension) for input files
}



### Reservoir sampling
my @Rsvr_Array; # sample reservoir (array for TEXT files, array-in-array for FASTA/Q files)
my @Title; # store header for TEXT files to append to the file later on, see option flag '-title_skip'
my $Next_Fasta_ID; # for multi-line FASTA input files to store next entry header/ID line while parsing in subroutine 'get_fastx_entry'

# 1. fill reservoir
while (<$Input_Fh>) {
    if (/^\s*$/) { # skip empty lines in input
        warn "\n### FASTQ file includes empty lines, which is unusual. Sampling FASTQ reads might fail so check the output file afterwards if the script didn't quit with a fatal error. However, consider running the input FASTQ file through 'fix_fastx.pl'!\n\n" if ($File_Type =~ /fastq/i);
        next;
    }
    chomp;
    get_file_type('stdin', $_) if (!$File_Type && $Input_File eq '-' && $. == 1); # subroutine to determine file type; mode 'stdin' piped input with first line of file
    warn "\n### Ignoring option '-t|-title_skip' as it has no effect on FASTA/Q files! Use the option only for TEXT files!\n\n" if($Title_Skip && $File_Type !~ /text/i && $. == 1);

    # FASTA file
    if ($File_Type =~ /fasta/i) {
        $_ = get_fastx_entry($_); # subroutine to read one FASTA sequence entry (seq in multi-line or not), returns anonymous array

    # FASTQ file
    } elsif ($File_Type =~ /fastq/i) {
        $_ = get_fastx_entry($_); # subroutine to read one FASTQ read composed of FOUR mandatory lines, returns reference to array

    # "single-line" TEXT file
    } elsif ($File_Type =~ /text/i) {
        warn "\n### Sure this is a TEXT file? The first line suspiciously looks like a FASTA ID/header line:\n$_\nBut, proceeding with file type 'text' ...\n\n" if ($. == 1 && $_ =~ /^>.*/);
        warn "\n### Sure this is a TEXT file? The first line suspiciously looks like a FASTQ sequence ID line:\n$_\nBut, proceeding with file type 'text' ...\n\n" if ($. == 1 && $_ =~ /^@.+/);
        if ($Title_Skip) { # skip possible header lines of TEXT file
            while ($Title_Skip) {
                push(@Title, $_); # store header line
                chomp($_ = <$Input_Fh>);
                $Title_Skip--;
            }
        }
    }

    push(@Rsvr_Array, $_); # fill reservoir
    last if (@Rsvr_Array == $Sample_Num); # reservoir is filled
}
die "\n### Fatal error:\nInsufficient records in input for sample size '$Sample_Num'!\n" if (@Rsvr_Array < $Sample_Num);

# 2. randomly replace elements in the reservoir with a decreasing probability
srand($Seed) if ($Seed); # force seed value
while (<$Input_Fh>) {
    if (/^\s*$/) {
        warn "\n### FASTQ file includes empty lines, which is unusual. Sampling FASTQ reads might fail so check the output file afterwards if the script didn't quit with a fatal error. However, consider running the input FASTQ file through 'fix_fastx.pl'!\n\n" if ($File_Type =~ /fastq/i);
        next;
    }
    my $rand_num = int(rand($.)); # choose an integer between 0 and eof-1; inclusive because array zero-based

    # replace elements in reservoir array
    if ($rand_num < @Rsvr_Array) {
        chomp;

        # FASTA file
        if ($File_Type =~ /fasta/i) {
            $_ = get_fastx_entry($_); # subroutine

        # FASTQ file
        } elsif ($File_Type =~ /fastq/i) {
            $_ = get_fastx_entry($_); # subroutine
        }
        $Rsvr_Array[$rand_num] = $_; # TEXT files are single-line based

    } elsif ($rand_num >= @Rsvr_Array) { # skip residual entry/read lines of FASTA/Q files
        if ($File_Type =~ /fasta/i) {
            get_fastx_entry($_); # subroutine without storing returning anonymous array (to skip FASTA multi-line seq)

        } elsif ($File_Type =~ /fastq/i) { # skip second, third, and fourth line of FASTQ file
            <$Input_Fh>; <$Input_Fh>; <$Input_Fh>;
        }
    }
}
close $Input_Fh;



### Print subsample to STDOUT
if (@Title) { # put header back in TEXT output
    foreach (@Title) {
        print "$_\n";
    }
}
foreach (@Rsvr_Array) {
    if ($File_Type =~ /fast/i) { # for both FASTA/Q file types
        foreach (@$_) { # de-reference array-ref for iteration
            print "$_\n";
        }

    # single-line TEXT file
    } else {
        print "$_\n";
    }
}

exit;


###############
# Subroutines #
###############

### Get sequence entries/reads from FASTA/Q file
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
            return [$header, $seq] if (/^>/); # return anonymous array with current header and seq
            $seq .= $_; # concatenate multi-line seq
        }
        return [$header, $seq] if eof;

    # FASTQ: FOUR lines for each FASTQ read (seq-ID, sequence, qual-ID [optional], qual)
    } elsif ($File_Type =~ /fastq/i) {
        my @fastq_read;

        # read name/ID, line 1
        my $seq_id = $line;
        die "\n### Fatal error:\nThis read doesn't have a sequence identifier/read name according to FASTQ specs, it should begin with a '\@':\n$seq_id\n" if ($seq_id !~ /^@.+/);
        push(@fastq_read, $seq_id);
        $seq_id =~ s/^@//; # remove '@' to make comparable to $qual_id

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
}



### Determine file type (FASTA, FASTQ, or TEXT)
sub get_file_type {
    my ($mode, $line) = @_; # mode either 'ext' or 'stdin' ('stdin' needs first line of input file)

    # determine file type via file extension
    if ($mode eq 'ext') {
        if ($Input_File =~ /.+\.(fa|fas|fasta|ffn|fna|frn|fsa)$/) { # use "|fsa)(\.gz)*$" if unzip inside script
            $File_Type = 'fasta';
        } elsif ($Input_File =~ /.+\.(fastq|fq)$/) {
            $File_Type = 'fastq';
        } else {
            $File_Type = 'text';
        }

    # determine by looking at first line of file
    } elsif ($mode eq 'stdin') {
        if ($line =~ /^>.*/) {
            $File_Type = 'fasta';
        } elsif ($line =~ /^@.+/) {
            $File_Type = 'fastq';
        } else {
            $File_Type = 'text';
        }
    }

    print STDERR "Detected file type: $File_Type\n";
    return 1;
}
