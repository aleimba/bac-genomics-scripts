#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<revcom_seq.pl> - reverse complement (multi-)sequence files

=head1 SYNOPSIS

C<perl revcom_seq.pl seq-file.embl E<gt> seq-file_revcom.embl>

B<or>

C<perl cat_seq.pl multi-seq_file.embl | perl revcom_seq.pl -i embl
E<gt> seq_file_cat_revcom.embl>

=head1 DESCRIPTION

This script reverse complements (multi-)sequence files. The
features/annotations in RichSeq files (e.g. EMBL or GENBANK format)
will also be adapted accordingly. Use option B<-o> to specify a
different output sequence format. Input files can be given directly via
C<STDIN> or as a file. If C<STDIN> is used, the input sequence file
format has to be given with option B<-i>. Be careful to set the correct
input format.

=head1 OPTIONS

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-o>=I<str>, B<-outformat>=I<str>

Specify different sequence format for the output [fasta, embl, or gbk]

=item B<-i>=I<str>, B<-informat>=I<str>

Specify the input sequence file format, only needed for C<STDIN> input

=item B<-v>, B<-version>

Print version number to C<STDOUT>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The reverse complemented sequence file is printed to C<STDOUT>.
Redirect or pipe into another tool as needed.

=back

=head1 EXAMPLES

=over

=item C<perl revcom_seq.pl -o gbk seq-file.embl E<gt>
seq-file_revcom.gbk>

=back

B<or>

=over

=item C<for file in *.embl; do perl revcom_seq.pl -o fasta "$file"
E<gt> "${file%.embl}"_revcom.fasta; done>

=back

=head1 DEPENDENCIES

=over

=item B<L<BioPerl|http://www.bioperl.org>>

Tested with BioPerl version 1.007001

=back

=head1 VERSION

 0.2                                               update: 2015-12-10
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
my $In_Format; # input seq file format needed for STDIN
my $Out_Format; # optional different output seq file format
my $VERSION = 0.2;
my ($Opt_Version, $Opt_Help);
GetOptions ('informat=s' => \$In_Format,
            'outformat=s' => \$Out_Format,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help)
            or pod2usage(-verbose => 1, -exitval => 2);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
if ($Opt_Version) {
    print "$0 $VERSION\n";
    exit;
}



### Check input (@ARGV and STDIN)
if (-t STDIN && ! @ARGV) {
    my $warning = "\n### Fatal error: No STDIN and no input file given as argument, please supply one of them and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
} elsif (!-t STDIN && @ARGV) {
    my $warning = "\n### Fatal error: Both STDIN and an input file given as argument, please supply only either one and/or see help with '-h'!\n";
    pod2usage(-verbose => 0, -message => $warning, -exitval => 2);
}
die "\n### Fatal error: Too many arguments given, only STDIN or one input file allowed as argument! Please see the usage with option '-h' if unclear!\n" if (@ARGV > 1);
die "\n### Fatal error: File '$ARGV[0]' does not exist!\n" if (@ARGV && $ARGV[0] ne '-' && !-e $ARGV[0]);



### Bio::SeqIO objects for input and output
print STDERR "\nReverse complementing";
my $Seqin; # Bio::SeqIO object
if (-t STDIN) { # input from file
    warn "\n### Warning: Ignoring input file format ('-i $In_Format'), because input file given and not STDIN!\n\n" if ($In_Format);
    my $seq_file = shift;
    $Seqin = Bio::SeqIO->new(-file => "<$seq_file"); # Bio::SeqIO object; no '-format' given, leave it to bioperl guessing
    print STDERR " '$seq_file' ";
} elsif (!-t STDIN) { # input from STDIN
    die "\n### Fatal error: Sequence file given as STDIN requires an input file format, please set one with option '-i' and/or see help with '-h'!\n" if (!$In_Format);
    $In_Format = 'genbank' if ($In_Format =~ /(gbk|gb)/i); # allow shorter format string for 'genbank'
    $Seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => $In_Format); # capture typeglob of STDIN, requires '-format'
    print STDERR " input file ";
}
print STDERR "...\n";

my $Seqout; # Bio::SeqIO object
if ($Out_Format) {
    $Out_Format = 'genbank' if ($Out_Format =~ /(gbk|gb)/i);
} else { # same format as input file
    if (!-t STDIN) {
        $Out_Format = $In_Format;
    } else {
        if (ref($Seqin) =~ /Bio::SeqIO::(genbank|embl|fasta)/) { # from bioperl guessing
            $Out_Format = $1;
        } else {
            die "\n### Fatal error: Could not determine input file format, please set an output file format with option '-o'!\n";
        }
    }
}
$Seqout = Bio::SeqIO->new(-fh => \*STDOUT, -format => $Out_Format); # printing to STDOUT requires '-format'


### Write reverse complemented sequence (and its features) to STDOUT
while (my $seq_obj = $Seqin->next_seq) { # Bio::Seq object; for multi-seq files
    my $revcom = Bio::SeqUtils->revcom_with_features($seq_obj);
    $Seqout->write_seq($revcom);
}

exit;
