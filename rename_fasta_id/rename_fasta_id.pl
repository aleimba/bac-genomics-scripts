#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<rename_fasta_id.pl> - rename fasta IDs according to regular expressions

=head1 SYNOPSIS

C<perl rename_fasta_id.pl -i file.fasta -p "NODE_.+$" -r "K-12_" -n -a c E<gt> out.fasta>

B<or>

C<zcat file.fasta.gz | perl rename_fasta_id.pl -i - -p "coli" -r "" -o E<gt> out.fasta>

=head1 DESCRIPTION

This script uses the built-in Perl substitution operator C<s///> to
replace strings in FASTA IDs. To do this, a B<pattern> and a
B<replacement> have to be provided (Perl regular expression syntax
can be used). The leading '>' character for the FASTA ID will be
removed before the substitution and added again afterwards. FASTA
IDs will be searched for matches with the B<pattern>, and if found
the B<pattern> will be replaced by the B<replacement>.

B<IMPORTANT>: Enclose the B<pattern> and the B<replacement> in
quotation marks (' or ") if they contain characters that would be
interpreted by the shell (e.g. pipes '|', brackets etc.).

For substitutions without any appendices in a UNIX OS you can of
course just use the great
L<C<sed>|https://www.gnu.org/software/sed/manual/sed.html> (see
C<man sed>), e.g.:

C<sed 's/^E<gt>pattern/E<gt>replacement/' file.fasta>

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Input FASTA file or piped STDIN (-) from a gzipped file

=item B<-p>=I<str>, B<-pattern>=I<str>

Pattern to be replaced in FASTA ID

=item B<-r>=I<str>, B<-replacement>=I<str>

Replacement to replace the pattern with. To entirely remove the
pattern use '' or "" as input for B<-r>.

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-c>, B<-case-insensitive>

Match pattern case-insensitive

=item B<-g>, B<-global>

Replace pattern globally in the string

=item B<-n>, B<-numerate>

Append a numeration/the count of the pattern hits to the
replacement. This is e.g. useful to number contigs consecutively in
a draft genome.

=item B<-a>=I<str>, B<-append>=I<str>

Append a string after the numeration, e.g. 'c' for chromosome

=item B<-o>, B<-output>

Verbose output of the substitutions that were carried out, printed
to C<STDERR>

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The FASTA file with substituted ID lines is printed to C<STDOUT>.
Redirect or pipe into another tool as needed.

=back

=head1 EXAMPLES

=over

=item C<perl rename_fasta_id.pl -i file.fasta -p "T" -r "a" -c -g -o>

=back

=head1 VERSION

 0.1                                                 09-11-2014

=head1 AUTHOR

 Andreas Leimbach                                    aleimba[at]gmx[dot]de

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
use autodie;
use Getopt::Long;
use Pod::Usage;

### Get the options with Getopt::Long
my $Input_File; # input fasta file
my $Pattern; # pattern to search for in the FASTA IDs
my $Replacement; # regex to replace pattern with
my $Opt_Case; # substitute case-insensitive
my $Opt_Global; # substitute pattern globally in string
my $Opt_Numerate; # append the count of the performed substitions to each replacement regex
my $Append; # append an additional string after $Opt_Numerate
my $Opt_Output; # print substitutions to STDERR
my $VERSION = 0.1;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_File,
            'pattern=s' => \$Pattern,
            'replacement=s' => \$Replacement,
            'case-insensitive' => \$Opt_Case,
            'global' => \$Opt_Global,
            'numerate' => \$Opt_Numerate,
            'append:s' => \$Append,
            'output' => \$Opt_Output,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_File || !$Pattern) {
    my $warning = "\n### Fatal error: Options '-i' or '-p' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Pipe input from STDIN or open input file
my $Input_Fh;
if ($Input_File eq '-') { # file input via STDIN
    $Input_Fh = *STDIN; # capture typeglob of STDIN
} else { # input via input file
    open ($Input_Fh, "<", "$Input_File");
}



### Parse FASTA file
my $Substitution_Count = 0; # count substitutions
while (<$Input_Fh>) {
    chomp;

    # only substitute in FASTA ID lines
    if (/^>/) {
        # only substitute if pattern found, case-sensitive or case-INsensitive
        if (/$Pattern/ || (/$Pattern/i && $Opt_Case)) {
            $_ = substitute_string($_); # subroutine

        # "reprint" FASTA IDs, which don't fit the pattern
        } else {
            print "$_\n";
        }

    # "reprint" sequence/non-ID lines of FASTA files
    } else {
        print "$_\n";
    }
}
print STDERR "$Substitution_Count substitutions have been carried out\n";

exit;


#############
#Subroutines#
#############

### Subroutine to rename headers/ID lines of the FASTA file
sub substitute_string {
    my $string = shift;
    $string =~ s/^>//; # get rid of '>', append afterwards

    print STDERR "$string " if ($Opt_Output); # optional verbose output to STDERR
    $Substitution_Count++; # count occurences of carried out substitutions

    # substitutions
    if ($Opt_Global && $Opt_Case) {
        $string =~ s/$Pattern/$Replacement/gi;
    } elsif ($Opt_Case) {
        $string =~ s/$Pattern/$Replacement/i;
    } elsif ($Opt_Global) {
        $string =~ s/$Pattern/$Replacement/g;
    } else {
        $string =~ s/$Pattern/$Replacement/;
    }

    # output to STDOUT, optionally STDERR
    print ">$string";
    print STDERR "-> $string" if ($Opt_Output);
    if ($Opt_Numerate) {
        print "$Substitution_Count";
        print STDERR "$Substitution_Count" if ($Opt_Output);
    }

    if ($Append) {
        print "$Append";
        print STDERR "$Append" if ($Opt_Output);
    }

    print "\n";
    print STDERR "\n" if ($Opt_Output);

    return 1;
}
