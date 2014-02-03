#!/usr/bin/perl

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Bio::SeqIO; # bioperl module to handle sequence input/output

my $usage = << "USAGE";

  ##################################################################
  # $0 -i seq_file -f in_format -o out_format #
  #                                                                #
  # Converts a (multi-)sequence file of a specific format to a     #
  # differently formatted output file, with the help of BioPerl    #
  # (www.bioperl.org).                                             #
  # Formats are e.g. embl, fasta, gbk.                             #
  #                                                                #
  # Mandatory options:                                             #
  # -i, -input       input sequence file                           #
  # -f, -format      input format                                  #
  # -o, -out_format  output format                                 #
  # Optional options:                                              #
  # -h, -help        print usage                                   #
  # -v, -version     print version number                          #
  #                                                                #
  # Adjust unix loop to run the script with all files in the       #
  # current working directory, e.g.:                               #
  # for i in *.gbk; do perl seq_format_converter.pl -i \$i -f gbk \\ #
  # -o embl; done                                                  #
  #                                                                #
  # version 0.2, update: 03-02-2014                     A Leimbach #
  # 10-11-2011                               aleimba[at]gmx[dot]de #
  ##################################################################

USAGE
;

### Get options with Getopt::Long
my $infile; # input sequence file
my $in_format; # input sequence file format
my $out_format; # desired output file format
my $version = 0.2;
my ($opt_version, $opt_help);
GetOptions ('input=s' => \$infile,
            'format=s' => \$in_format,
            'out_format=s' => \$out_format,
            'version' => \$opt_version,
            'help|?' => \$opt_help);


### Print usage
if ($opt_help) {
    die $usage;
} elsif ($opt_version) {
    die "$0 $version\n";
} elsif (!$infile || !$in_format || !$out_format) {
    die $usage, "### Fatal error: Option(s) or argument(s) for \'-i\', \'-f\', \'-o\' are missing!\n\n";
}


### Allow shorter format string for 'genbank'
$in_format = 'genbank' if ($in_format =~ /gbk/i);
my $outfile = $infile;
$outfile =~ s/\.\w+$/\.$out_format/; # remove file extension from infile and append out_format
$out_format = 'genbank' if ($out_format =~ /gbk/i);


### SeqIO objects for input and output
my $seq_in = Bio::SeqIO->new(-file => "<$infile", -format => $in_format); # a Bio::SeqIO object
my $seq_out = Bio::SeqIO->new(-file => ">$outfile", -format => $out_format); # a Bio::SeqIO object


### Write sequence to different format
while (my $seqobj = $seq_in->next_seq) { # a Bio::Seq object
    $seq_out->write_seq($seqobj);
}
print "\n\tCreated new file $outfile!\n\n";

exit;
