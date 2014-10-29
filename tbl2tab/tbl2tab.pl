#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<tbl2tab.pl> - convert tbl to tab-separated format and back

=head1 SYNOPSIS

C<perl tbl2tab.pl -m tbl2tab -i feature_table.tbl -s -l locus_prefix>

B<or>

C<perl tbl2tab.pl -m tab2tbl -i feature_table.tab -g -l locus_prefix
-p "gnl|dbname|">

=head1 DESCRIPTION

NCBI's feature table (B<tbl>) format is needed for the submission of
genomic data to GenBank with the NCBI tools
L<Sequin|http://www.ncbi.nlm.nih.gov/Sequin/> or
L<tbl2asn|http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2>. tbl files
can be created with automatic annotation systems like
L<Prokka|http://www.vicbioinformatics.com/software.prokka.shtml>.
C<tbl2tab.pl> can convert a tbl file to a tab-separated format (tab)
and back to the tbl format. The tab-delimited format is useful to
manipulate the data more comfortably in a spreadsheet software (e.g.
LibreOffice or MS Excel). For a conversion back to tbl format save
the file in the spreadsheet software as a tab-delimited text file.
The script is intended for microbial genomes, but might also be
useful for eukaryotes.

Regular expressions are applied in mode B<tbl2tab> to correct gene
names and words in '/product' values  to lowercase initials (with
the exception of 'Rossman' and 'Willebrand'). The resulting tab file
can then be used to check for possible errors.

The first four header columns of the B<tab> format are mandatory,
'seq_id' for the SeqID, and for each primary tag/feature (e.g. CDS,
RNAs, repeat_region etc.), 'start', 'stop', and 'primary_tag'. These
mandatory columns have to be filled in every row in the tab file.
All the following columns will be included as tags/qualifiers (e.g.
'/locus_tag', '/product', '/EC_number', '/note' etc.) in the
conversion to the tbl file if a value is present.

There are three special cases:

B<First>, '/pseudo' will be included as a tag if I<any> value (the
script uses 'T' for true) is present in the B<tab> format. If a
primary tag is indicated as pseudo both the primary tag and the
accessory 'gene' primary tag (for CDS/RNA features with option
B<-g>) will include a '/pseudo' qualifier in the resulting B<tbl>
file. B<Pseudo-genes> are indicated by 'pseudo' in the 'primary_tag'
column, thus the 'pseudo' column is ignored in these cases.

B<Second>, tag '/gene_desc' is reserved for the 'product' values of
pseudo-genes, thus a 'gene_desc' column in a tab file will be
ignored in the conversion to tbl.

B<Third>, column 'protein_id' in a tab file will also be ignored in
the conversion. '/protein_id' values are created from option B<-p>
and the locus_tag for each CDS primary feature.

Furthermore, with option B<-s> G2L-style spreadsheet formulas
(L<Goettingen Genomics
Laboratory|http://appmibio.uni-goettingen.de/>) can be included with
additional columns, 'spreadsheet_locus_tag', 'position', 'distance',
'gene_number', and 'contig_order'. These columns will not be
included in a conversion to the tbl format. Thus, if you want to
include e.g. the locus_tags from the formula in column
'spreadsheet_locus_tag' in the resulting tbl file copy the B<values>
to the column 'locus_tag'!

To illustrate the process two example files are included in the
repository, F<example.tbl> and F<example2.tab>, which are
interconvertible (see L</"EXAMPLES"> below).

B<Warning>, be aware of possible errors introduced by automatic
format conversions using a spreadsheet software like MS Excel, see
e.g. Zeeberg et al. 2004
(L<http://www.ncbi.nlm.nih.gov/pubmed/15214961>).

For more information regarding the feature table and the submission
process see NCBI's L<prokaryotic annotation
guide|http://www.ncbi.nlm.nih.gov/genbank/genomesubmit> and the
L<bacterial genome submission
guide|http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation>.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-m>=I<tbl2tab|tab2tbl>, B<-mode>=I<tbl2tab|tab2tbl>

Conversion mode, either 'tbl2tab' or 'tab2tbl' [default = 'tbl2tab']

=item B<-i>=I<str>, B<-input>=I<str>

Input tbl or tab file to be converted to the other format

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head3 Mode B<tbl2tab>

=over 20

=item B<-l>=I<str>, B<-locus_prefix>=I<str>

Only in combination with option B<-s> and there mandatory to include
the locus_tag prefix in the formula for column 'spreadsheet_locus_tag'

=item B<-c>, B<-concat>

Concatenate values of identical tags within one primary tag with '~'
(e.g. several '/EC_number' or '/inference' tags)

=item B<-e>=I<str>, B<-empty>=I<str>

String used for primary features without value for a tag [default = '']

=item B<-s>, B<-spreadsheet>

Include formulas for spreadsheet editing

=item B<-f>=I<e|g>, B<-formula_lang>=I<e|g>

Syntax language of the spreadsheet formulas, either 'English' or
'German'. If you're still encountering problems with the formulas
set the decimal and thousands separator manually in the options of
the spreadsheet software (instead of using the operating system
separators). [default = 'e']

=back

=head3 Mode B<tab2tbl>

=over 20

=item B<-l>=I<str>, B<-locus_prefix>=I<str>

Prefix to the SeqID if not present already in the SeqID

=item B<-g>, B<-gene>

Include accessory 'gene' primary tags (with '/gene', '/locus_tag'
and possibly '/pseudo' tags) for 'CDS/RNA' primary tags; NCBI standard

=item B<-t>, B<-tags_full>

Only in combination with option B<-g>, include '/gene' and
'/locus_tag' tags additionally in primary tag, not only in accessory
'gene' primary tag

=item B<-p>=I<str>, B<-protein_id_prefix>=I<str>

Prefix for '/protein_id' tags; don't forget the double quotes for
the string, otherwise the shell will intepret as pipe [default =
'gnl|goetting|']

=back

=head1 OUTPUT

=over 20

=item F<*.tab|tbl>

Result file in the opposite format

=item (F<hypo_putative_genes.txt>)

Created in mode 'tab2tbl', indicates if CDSs are annotated as
'hypothetical/putative/predicted protein' but still have a gene name

=back

=head1 EXAMPLES

=over

=item C<perl tbl2tab.pl -m tbl2tab -i example.tbl -s -l EPE>

=item C<perl tbl2tab.pl -m tab2tbl -i example2.tab -g -l EPE>

=back

=head1 VERSION

 0.2                                               update: 29-10-2014
 0.1                                                       24-06-2014

=head1 AUTHOR

 Andreas Leimbach                               aleimba[at]gmx[dot]de

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



### Get options with Getopt::Long, works also abbreviated and with two "--": -i, --i, -input ...
my $Input_File; # input file
my $Mode = 'tbl2tab'; # mode of script, i.e. either convert from tbl2tab or from tab2tbl; default 'tbl2tab'
my $Locus_Prefix = ''; # required for option 'spreadsheet' in mode 'tbl2tab', in mode 'tab2tbl' optional
my $Opt_Concat; # optionally, concatenate values of the same tag within one primary tag in a tbl file in one column in the resulting tab file with '~' (e.g. several 'EC_number' tags etc.)
my $Empty = ''; # optionally, set what should be used for tags without a value in resulting tab file; default is nothing
my $Opt_Spreadsheet; # optionally, include formulas for spreadsheet editing (e.g. Libre Office, MS Excel)
my $Formula_Lang_Spreadsheet = 'e'; # optionally, either German or English formulas in Spreadsheet option; default 'e' for English
my $Opt_Gene; # optionally, include accessory gene primary tags (with '/gene' and '/locus_tag' [and '/pseudo'] tags) for CDS|RNA primary tags
my $Opt_Tags_Full; # optionally, include '/gene' and '/locus_tag' additionally in primary tag not only accessory 'gene' primary tag
my $Protein_Id_Prefix = 'gnl|goetting|'; # optionally give a different string to prefix the '/protein_id' tags
my $VERSION = 0.2;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_File,
            'mode=s' => \$Mode,
            'locus_prefix:s' => \$Locus_Prefix,
            'concat' => \$Opt_Concat,
            'empty:s' => \$Empty,
            'spreadsheet' => \$Opt_Spreadsheet,
            'formula_lang:s' => \$Formula_Lang_Spreadsheet,
            'gene' => \$Opt_Gene,
            'tags_full' => \$Opt_Tags_Full,
            'protein_id_prefix:s' => \$Protein_Id_Prefix,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_File) {
    my $warning = "\n### Fatal error: Option '-i' or its argument is missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
} elsif (!($Mode =~ /tbl2tab/i || $Mode =~ /tab2tbl/i)) { # case-insensitive
    my $warning = "\n### Fatal error: Incorrect run mode with option '-m' given! Please choose either 'tbl2tab' or 'tab2tbl' for '-m'!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Enforce mandatory or optional options
if ($Mode =~ /tbl2tab/i && ($Opt_Gene || $Opt_Tags_Full || $Protein_Id_Prefix ne 'gnl|goetting|')) {
    warn "\nIncompatible option(s) '-g', '-p', or '-t' set with mode 'tbl2tab'. Ignoring the option(s)!\n";
} elsif ($Mode =~ /tab2tbl/i && ($Opt_Concat || $Empty || $Opt_Spreadsheet || $Formula_Lang_Spreadsheet ne 'e')) {
    warn "\nIncompatible option(s) '-c', '-e', '-f', or '-s' set with mode 'tab2tbl'. Ignoring the option(s)!\n";
    $Formula_Lang_Spreadsheet = 'e' if ($Formula_Lang_Spreadsheet ne 'e'); # avoid die with error below if option not set correctly
}

if ($Mode =~ /tbl2tab/i && ($Locus_Prefix || $Formula_Lang_Spreadsheet ne 'e') && !$Opt_Spreadsheet) {
    warn "\nOption(s) '-l' or '-f' set, but not '-s'. Forcing option '-s'!\n";
    $Opt_Spreadsheet = 1;
}
if ($Mode =~ /tbl2tab/i && !$Locus_Prefix && $Opt_Spreadsheet) {
    warn "\nOption '-s' set, but not '-l'. Please give a prefix for the locus tags: ";
    chomp($Locus_Prefix = <>);
}

if ($Formula_Lang_Spreadsheet !~ /^(e|g)/i) {
    die "\n### Fatal error: Incorrect language for option '-f' given! Please choose either 'e|eng' or 'g|ger' for '-f'!\n";
}

if ($Mode =~ /tab2tbl/i && !$Opt_Gene && $Opt_Tags_Full) {
    warn "\nOption '-t' set, but not '-g'. Forcing option '-g'!\n";
    $Opt_Gene = 1;
}



### Read in tbl or tab-separated data and write to result file in the opposite format
my $Out_File = $Input_File;
$Out_File =~ s/^(.+)\.\w+$/$1/; # strip filename extension
my $Error_File = 'hypo_putative_genes.txt';
if ($Mode =~ /tbl2tab/i) {
    my ($data_hash_ref, $tags_max_count_hash_ref) = read_tbl(); # subroutine
    $Out_File .= '.tab';
    write_tab($data_hash_ref, $tags_max_count_hash_ref); # subroutine

} elsif ($Mode =~ /tab2tbl/i) {
    $Out_File .= '.tbl';
    read_tab_write_tbl(); # subroutine
}



### Message which file was created
if ($Mode =~ /tbl2tab/i) {
    print "Input tbl file '$Input_File' was converted to tab output file '$Out_File'!\n";
} elsif ($Mode =~ /tab2tbl/i) {
    print "Input tab file '$Input_File' was converted to tbl output file '$Out_File'!\n";
}

if (-e $Error_File) {
    if (-s $Error_File >= 40) { # smaller than just the header, which should be 27 bytes
        warn "\n### Warning: CDSs found that are annotated with 'hypothetical|putative|predicted protein' but still include a '/gene' tag, see file '$Error_File'!\n";
    } elsif (-s $Error_File < 40) {
        unlink $Error_File;
    }
}


exit;

###############
# Subroutines #
###############

### Subroutine to test for file existence and give warning to STDERR
sub file_exist {
    my $file = shift;
    if (-e $file) {
        warn "\nThe result file \'$file\' exists already and will be overwritten!\n\n";
        return 1;
    }
    return 0;
}



### Print tag values to tab result file
sub print_tag2tab {
    my ($tag, $data_hash_ref, $seq_id, $pos, $tag_max_count) = @_;

    if ($Opt_Concat) { # values concatenated by '~' in $data_hash_ref
        if ($data_hash_ref->{$seq_id}->{$pos}->{$tag}) {
            print "\t$data_hash_ref->{$seq_id}->{$pos}->{$tag}";
            return 1;
        } else {
            print "\t$Empty";
            return 1;
        }

    } elsif (!$Opt_Concat) { # split concatenated values in individual values
        my @values = split(/~/, $data_hash_ref->{$seq_id}->{$pos}->{$tag}) if ($data_hash_ref->{$seq_id}->{$pos}->{$tag});
        if (@values) {
            foreach (@values) {
                print "\t$_";
            }
            print "\t$Empty" x ($tag_max_count - @values); # fill residual columns till maximum occurrence
        } else {
            print "\t$Empty" x $tag_max_count;
        }
        return 1;
    }

    return 0;
}



### Print tag values to tbl result file
sub print_tag2tbl {
    my ($tag, $value) = @_;
    return 0 if ($value =~ /^$Empty$/);
    if ($tag =~ /pseudo/) {
        print "\t\t\t$tag\n";
        return 1;
    }

    ### remove quotations from values introduced by Excel by saving as tab-separated file:
    ### https://office.microsoft.com/en-001/excel-help/excel-formatting-and-features-that-are-not-transferred-to-other-file-formats-HP010014105.aspx
    ### - if a cell contains a comma, the cell contents are enclosed in double quotation marks
    ### - if the data contains a quotation mark, double quotation marks will replace the quotation mark, and the cell contents are also enclosed in double quotation marks
    $value =~ s/""/"/g;
    $value =~ s/^"//;
    $value =~ s/"$//;

    foreach (split(/~/, $value)) {
        print "\t\t\t$tag\t$_\n";
    }
    return 1;
}


### Read in data from tab input file and write it to tbl output file
sub read_tab_write_tbl {
    file_exist($Error_File); # subroutine
    open (my $error_file_fh, ">", $Error_File);
    print $error_file_fh "row\tlocus_tag\tgene\tproduct\n";

    open (my $input_file_fh, "<", $Input_File);
    my $header = <$input_file_fh>;
    $header =~ s/\R/\012/; # convert line to unix-style line endings
    chomp $header;
    if ($header !~ /^seq_id\tstart\tstop\tprimary_tag\t/) { # check if tbl file starts with mandatory header fields or quit
        die "\n### Fatal error: Input tab file '$Input_File' doesn't start with the mandatory 'seq_id', 'start', 'stop', and 'primary_tag' tab-separated header fields. Sure this is a valid tab file?\nExiting program!\n\n";
    }

    my @tags;
    foreach (split(/\t/, $header)) {
        last if (/spreadsheet_locus_tag/); # skip all optional extra spreadsheet columns
        push(@tags, $_); # store all header fields/columns to associate with each field in each line
    }

    file_exist($Out_File); # subroutine
    open (my $out_file_fh, ">", $Out_File);
    select $out_file_fh; # select fh for standard print/f output

    my $row = 1; # count row numbers of tab input file for $Error_File (start with '1' as header already parsed above)
    my $seq_id = ''; # store previous SeqID for multi-contig/replicon tab files
    while (<$input_file_fh>) {
        $row++;
        $_ =~ s/\R/\012/; # convert line to unix-style line ending
        chomp;
        next if ($_ =~ /^\s+$/ || $_ =~ /^$/); # skip empty lines

        my ($locus_tag, $gene, $hypo_putative) = ('', '', ''); # needed for $Error_File

        my @cells = split(/\t/, $_);
        for (my $i = 0; $i < 4; $i++) { # check each row for mandatory fields
            if ($cells[$i] =~ /^$/) {
                close $out_file_fh;
                unlink $Out_File;
                die "\n### Fatal error: Row $row of input tab file '$Input_File' is missing a value for one of the mandatory fields 'seq_id', 'start', 'stop', or 'primary_tag'!\nExiting program!\n\n";
            }
        }

        ### print SeqID
        if ($cells[0] ne $seq_id) { # print new contig for multi-contig/replicon tab files
            $seq_id = $cells[0];
            $cells[0] = $Locus_Prefix."_".$cells[0] if ($cells[0] !~ /$Locus_Prefix/ && $Locus_Prefix); # append locus_tag prefix only if not present already
            print ">Feature $cells[0]\n";
        }

        ### print accessory 'gene' primary tags with '/locus_tag', '/gene', and potential '/pseudo' tags
        if ($Opt_Gene && $cells[3] =~ /CDS|RNA/) { # accessory 'gene' primary tags only for CDS and RNA (rRNA, tRNA ...) features
            print "$cells[1]\t$cells[2]\tgene\n";
            my $column_count = 0;
            foreach my $tag (@tags) {
                if ($tag =~ /locus_tag/ && $cells[$column_count] =~ /^$/) { # CDSs|RNAs mandatory need a '/locus_tag' with option '-g'
                    close $out_file_fh;
                    unlink $Out_File;
                    die "\n### Fatal error: Row $row of input tab file '$Input_File' is missing a 'locus_tag' which is mandatory for option '-g'!\nExiting program!\n\n";
                }
                print_tag2tbl($tag, $cells[$column_count]) if ($tag =~ /locus_tag|^gene$|pseudo/); # subroutine; '^gene$' needed, so 'gene_desc' isn't hit (see below)
                $column_count++;
            }
        }

        ### print primary tag
        print "$cells[1]\t$cells[2]"; # start\tstop
        if ($cells[3] =~ /pseudo/) { # pseudo-gene should include '/pseudo' tag
            print "\tgene\n";
            print "\t\t\tpseudo\n";
        } else {
            print "\t$cells[3]\n";
        }

        ### print tags with values
        for (my $i = 4; $i < @tags; $i++) { # start with field 5 of array with header fields/columns (the first 4 are mandatory, see above)
            next if ($tags[$i] =~ /gene_desc/); # skip 'gene_desc' fields in tab file, reserved for pseudo-genes (see below)

            ### enforce mandatory tags for CDS primary tags
            if ($tags[$i] =~ /product/ && $cells[3] =~ /CDS/) {
                if ($cells[$i] =~ /(hypothetical|putative|predicted) protein/) { # needed for $Error_File
                    $hypo_putative = $cells[$i];
                } elsif ($cells[$i] =~ /^$/) { # CDSs mandatory need a value for '/product'
                    close $out_file_fh;
                    unlink $Out_File;
                    die "\n### Fatal error: Row $row of input tab file '$Input_File' is missing a 'product' value which is mandatory for CDS primary tags!\nExiting program!\n\n";
                }
            }

            if ($tags[$i] =~ /locus_tag/ && $cells[3] =~ /CDS/) { # '/protein_id' mandatory for 'CDS' primary tags
                $locus_tag = $cells[$i]; # needed for $Error_File
                my $protein_id = "$Protein_Id_Prefix".$cells[$i];
                print_tag2tbl('protein_id', $protein_id); # subroutine
            }
            next if ($tags[$i] =~ /protein_id/); # skip 'protein_id' field in tab file as they should be created from the 'locus_tag' column

            $gene = $cells[$i] if ($tags[$i] =~ /^gene$/ && $cells[3] =~ /CDS/); # needed for $Error_File

            ### enforce mandatory tags for CDS/RNA primary tags
            next if ($tags[$i] =~ /locus_tag|^gene$/ && $Opt_Gene && !$Opt_Tags_Full && $cells[3] =~ /CDS|RNA/); # skip '/locus_tag' and '/gene' tags if accessory gene primary tags are present for CDS|RNA features (except option '-t' is set)
            if ($tags[$i] =~ /product/ && $cells[3] =~ /RNA/ && $cells[$i] =~ /^$/) { # RNAs mandatory need a value for '/product'
                close $out_file_fh;
                unlink $Out_File;
                die "\n### Fatal error: Row $row of input tab file '$Input_File' is missing a 'product' value which is mandatory for RNA primary tags!\nExiting program!\n\n";
            }

            ### enforce mandatory tag for pseudo-genes (have 'pseudo' as primary tag in tab file)
            if ($tags[$i] =~ /locus_tag/ && $Opt_Gene && $cells[3] =~ /pseudo/ && $cells[$i] =~ /^$/) { # pseudo-genes mandatory need a '/locus_tag' with option '-g'
                close $out_file_fh;
                unlink $Out_File;
                die "\n### Fatal error: Row $row of input tab file '$Input_File' is missing a 'locus_tag' which is mandatory for option '-g'!\nExiting program!\n\n";
            }

            ### the rest
            if ($tags[$i] =~ /product/ && $cells[3] =~ /pseudo/) { # write 'product' values for pseudo-genes to '/gene_desc' tags
                print_tag2tbl('gene_desc', $cells[$i]); # subroutine
            } elsif ($tags[$i] =~ /pseudo/ && $cells[3] =~ /pseudo/) { # skip 'pseudo' tag if pseudo-gene
                next;
            } else {
                print_tag2tbl($tags[$i], $cells[$i]); # subroutine
            }

        }
        print $error_file_fh "$row\t$locus_tag\t$gene\t$hypo_putative\n" if ($hypo_putative && $gene);
    }

    select STDOUT;
    close $input_file_fh;
    close $error_file_fh;
    close $out_file_fh;
    return 1;
}



### Read in data from tbl input file
sub read_tbl {
    open (my $input_file_fh, "<", $Input_File);

    my $seq_id = <$input_file_fh>; # SeqID
    $seq_id =~ s/\R/\012/; # convert line to unix-style line endings
    chomp $seq_id;
    if ($seq_id !~ /^>Feature/) { # check if tbl file starts with mandatory '>Feature' and get first SeqID
        die "\n### Fatal error: tbl file doesn't start with a '>Feature SeqID' line. Sure this is a valid tbl file?\nExiting program!\n\n";
    } else {
        $seq_id =~ s/>Feature (\S+)\s*$/$1/; # only use non-whitespace characters as SeqID
    }

    my %data; # hash-in-hash-in-hash to store tbl input data
    my $pos_key; # store start..stop for each primary tag and use as key in %data
    my $primary_tag; # store previous primary tag to determine if values for repeatedly occuring tags should be concatenated
    my %tags_max_count; # hash to store all occuring tags with maximal number of presence (within a single primary tag) in the tbl file for final tab column headers
    my @tags; # array to store all tags of each primary tag, supplement to %tags_max_count

    while (<$input_file_fh>) {
        $_ =~ s/\R/\012/; # convert line to unix-style line ending
        chomp;
        next if ($_ =~ /^\s+$/); # skip empty lines
        my @fields = split(/\t/, $_);

        ### get next SeqID from '>Feature' line for multi-contig/replicon tbl files
        if ($fields[0] =~ /^>Feature (\S+)\s*$/) {
            $seq_id = $1;

        ### get primary tags/features and fill %tags_max_count
        } elsif ($fields[0] =~ /^\d+$/) { # $fields[2] with primary tag
            foreach my $tag (@tags) { # fill %tags_max_count for previous primary tag
                if ($tags_max_count{$tag}) {
                    $tags_max_count{$tag} = grep(/$tag/, @tags) if ($tags_max_count{$tag} < grep(/$tag/, @tags));
                } elsif (!$tags_max_count{$tag}) {
                    $tags_max_count{$tag} = grep(/$tag/, @tags);
                }
            }
            @tags = (); # empty tags array for new current primary tag

            $pos_key = "$fields[0]..$fields[1]"; # position of primary tag used as key for %data, "start..stop"
            $primary_tag = $fields[2];
            if (!$data{$seq_id}->{$pos_key}->{'primary_tag'} || $data{$seq_id}->{$pos_key}->{'primary_tag'} =~ /gene/) { # if primary tag not present or overwrite accessory 'gene' primary tag
                $data{$seq_id}->{$pos_key}->{'primary_tag'} = $primary_tag; # store data in anonymous hash-in-hash-in-hash
                $data{$seq_id}->{$pos_key}->{'start'} = $fields[0]; # to be able to sort afterwards via the start position
            } elsif ($data{$seq_id}->{$pos_key}->{'primary_tag'} =~ /pseudo/) { # 'gene' primary tag with '/pseudo' tag will be replaced by 'pseudo' primary tag for pseudo-genes (see below), however if 'gene' primary tag is ACCESSORY to CDS|RNA primary tag replace by this primary tag and include '/pseudo' tag
                $data{$seq_id}->{$pos_key}->{'primary_tag'} = $primary_tag;
                $data{$seq_id}->{$pos_key}->{'pseudo'} = 'T'; # value 'T' for true
            }

        ### get tags/qualifiers
        } elsif ($fields[3] =~ /^\w+/) {
            push(@tags, $fields[3]) if ($fields[3] !~ /gene_desc/); # store tags for current primary tag; skip '/gene_desc' as reserved for pseudo-genes (replaced by '/product' see below)
            if ($fields[3] =~ /pseudo/) {
                if ($data{$seq_id}->{$pos_key}->{'primary_tag'} =~ /gene/) { # change 'gene' primary tag of pseudo-genes to 'pseudo' (if accessory 'gene' primary tag will be replaced by *actual* primary tag, see above)
                    $data{$seq_id}->{$pos_key}->{'primary_tag'} = 'pseudo';
                } else { # else include a '/pseudo' tag with value 'T' for true
                    $data{$seq_id}->{$pos_key}->{'pseudo'} = 'T';
                }
                next; # next line
            }

            ### remove quotations from values introduced by Excel by saving as tab-separated file (see above)
            $fields[4] =~ s/""/"/g;
            $fields[4] =~ s/^"//;
            $fields[4] =~ s/"$//;

            ### adjust '/gene' and '/product' values to NCBI standard
            if ($fields[3] =~ /gene/) {
                $fields[4] =~ s/(\w+)/\l$1/; # first character of gene name should be lower case
                $fields[4] =~ s/^(\w)$/\u$1/; # one letter phage genes should be upper case
            } elsif ($fields[3] =~ /product/) {
                $fields[4] =~ s/\b([A-Z][a-z]{3,})/\l$1/g; # lower the case for '/protein' value initials
                $fields[4] =~ s/(rossman|willebrand)/\u$1/; # exception to the rule
            }

            if ($fields[3] =~ /gene_desc/) { # '/gene_desc' tags from pseudo-genes replaced by '/product' for resulting tab file
                $data{$seq_id}->{$pos_key}->{'product'} = $fields[4];
                next;
            }
            if ($data{$seq_id}->{$pos_key}->{$fields[3]} && $data{$seq_id}->{$pos_key}->{'primary_tag'} =~ /$primary_tag/) { # tag already exists for this position (e.g. several EC_numbers), concatenate the additional values with '~' as separator only WITHIN the same primary tag (OTHERWISE overwrite in 'else' below)
                $data{$seq_id}->{$pos_key}->{$fields[3]} .= '~'.$fields[4];
            } else { # tag doesn't exist yet or overwrite if current primary tag at the same position of previous (e.g. accessory 'gene' primary tag to CDS/RNA primary tag)
                $data{$seq_id}->{$pos_key}->{$fields[3]} = $fields[4];
            }
        }
    }
    foreach my $tag (@tags) { # fill %tags_max_count for last primary tag
        if ($tags_max_count{$tag}) {
            $tags_max_count{$tag} = grep(/$tag/, @tags) if ($tags_max_count{$tag} < grep(/$tag/, @tags));
        } elsif (!$tags_max_count{$tag}) {
            $tags_max_count{$tag} = grep(/$tag/, @tags);
        }
    }
    close $input_file_fh;
    return \%data, \%tags_max_count;
}



### Write data to tab output file
sub write_tab {
    my ($data_hash_ref, $tags_max_count_hash_ref) = @_;
    file_exist($Out_File); # subroutine
    open (my $out_file_fh, ">", $Out_File);
    select $out_file_fh; # select fh for standard print/f output

    ### print header for tab result file
    print "seq_id\tstart\tstop\tprimary_tag\tlocus_tag"; # mandatory columns/fields in tab file
    if ($Opt_Concat) {
        foreach (sort keys %{$tags_max_count_hash_ref}) { # print residual tags
            print "\t$_" if (!/locus_tag/);
        }
    } elsif (!$Opt_Concat) {
        foreach (sort keys %{$tags_max_count_hash_ref}) {
            print "\t$_" x $tags_max_count_hash_ref->{$_} if (!/locus_tag/); # print max occurrence (in tbl) of each residual tag
        }
    }

    print "\tspreadsheet_locus_tag\tposition\tdistance\tgene_number\tcontig_order" if ($Opt_Spreadsheet); # print optional spreadsheet header columns
    print "\n";

    ### variables for optional spreadsheet formulas
    my @spread_columns = ("A".."AZ") if ($Opt_Spreadsheet); # columns in spreadsheet software for formulas
    my ($tags_column_count, $spread_row_count, $spread_contig_order) = (0, 1, 1) if ($Opt_Spreadsheet);
    if ($Opt_Spreadsheet) {
        if ($Opt_Concat) {
            $tags_column_count = (scalar keys %{$tags_max_count_hash_ref}) - scalar grep($_ =~ /locus_tag/, keys %{$tags_max_count_hash_ref}); # subtract tags for correct spreadsheet formulas
        } elsif (!$Opt_Concat) {
            foreach (keys %{$tags_max_count_hash_ref}) {
                next if ($_ =~ /locus_tag/);
                $tags_column_count += $tags_max_count_hash_ref->{$_};
            }
        }
    }

    ### print data from hash into tab result file, optional with G2L-style spreadsheet formulas
    foreach my $seq_id (sort keys %{$data_hash_ref}) {
        foreach my $pos (sort {$data_hash_ref->{$seq_id}->{$a}->{'start'} <=> $data_hash_ref->{$seq_id}->{$b}->{'start'}} keys $data_hash_ref->{$seq_id}) { # sort each position entry in %data via start position
            print "$seq_id";
            my ($start, $stop) = split(/\.\./, $pos);
            print "\t$start\t$stop";
            print "\t$data_hash_ref->{$seq_id}->{$pos}->{'primary_tag'}"; # primary_tag should always be present
            print_tag2tab('locus_tag', $data_hash_ref, $seq_id, $pos, 1); # subroutine; locus_tag should occur always just one time per primary tag
            foreach (sort keys %{$tags_max_count_hash_ref}) { # print residual tags
                print_tag2tab($_, $data_hash_ref, $seq_id, $pos, $tags_max_count_hash_ref->{$_}) if (!/locus_tag/); # subroutine
            }

            ### G2L-style spreadsheet formulas
            if ($Opt_Spreadsheet) {
                $spread_row_count++;
                if ($Formula_Lang_Spreadsheet =~ /^e/i) {
                    print "\t=\"$Locus_Prefix\"", '&"_"&A', "$spread_row_count&TEXT(", $spread_columns[$tags_column_count+8], $spread_row_count, ',"0000")&"0"'; # spreadsheet column 'spreadsheet_locus_tag'
                } elsif ($Formula_Lang_Spreadsheet =~ /^g/i) {
                    print "\t=\"$Locus_Prefix\"", '&"_"&A', "$spread_row_count&TEXT(", $spread_columns[$tags_column_count+8], $spread_row_count, ';"0000")&"0"';
                }
                print "\t=MIN(B$spread_row_count:C$spread_row_count)"; # spreadsheet column 'position'
                print "\t=", $spread_columns[$tags_column_count+6], $spread_row_count + 1, "-MAX(B$spread_row_count:C$spread_row_count)"; # spreadsheet column 'distance'
                if ($Formula_Lang_Spreadsheet =~ /^e/i) {
                    print "\t=IF(", $spread_columns[$tags_column_count+9], $spread_row_count, '=', $spread_columns[$tags_column_count+9], $spread_row_count - 1, ",$spread_columns[$tags_column_count+8]", $spread_row_count - 1, '+1,1)'; # spreadsheet column 'gene_number'
                } elsif ($Formula_Lang_Spreadsheet =~ /^g/i) {
                    print "\t=WENN(", $spread_columns[$tags_column_count+9], $spread_row_count, '=', $spread_columns[$tags_column_count+9], $spread_row_count - 1, ";$spread_columns[$tags_column_count+8]", $spread_row_count - 1, '+1;1)';
                }
                print "\t$spread_contig_order"; # spreadsheet column 'contig_order'
            }

            print "\n";
        }
        $spread_contig_order++; # next contig/replicon (SeqID) in tbl file
    }

    select STDOUT;
    close $out_file_fh;
    return 1;
}
