#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<genomes_feature_table.pl> - create a feature table for genomes in EMBL and GENBANK format

=head1 SYNOPSIS

C<perl genomes_feature_table.pl path/to/genome_dir E<gt> feature_table.tsv>

=head1 DESCRIPTION

A genome feature table lists basic stats/info (e.g. genome size, GC
content, coding percentage, accession number(s)) and the numbers of
annotated primary features (e.g. CDS, genes, RNAs) of genomes. It
can be used to have an overview of these features in different
genomes, e.g. in comparative genomics publications.

C<genomes_feature_table.pl> is designed to extract (or calculate)
these basic stats and B<all> annotated primary features from RichSeq
files (B<EMBL> or B<GENBANK> format) in a specified directory (with
the correct file extension, see option B<-e>). The B<default> directory
is the current working directory. The primary features are counted
and the results for each genome printed in tab-separated format. It
is a requirement that each file contains B<only one> genome
(complete or draft, with or without plasmids).

The most important features will be listed first, like genome
description, genome size, GC content, coding percentage (calculated
based on non-pseudo CDS annotation), CDS and gene numbers, accession
number(s) (first..last in the sequence file), RNAs (rRNA, tRNA,
tmRNA, ncRNA), and unresolved bases (IUPAC code 'N'). If plasmids are
annotated in a sequence file, the number of plasmids are counted
and listed as well (needs a I</plasmid="plasmid_name"> tag in the
I<source> primary tag, see e.g. Genbank accession number
L<CP009167|http://www.ncbi.nlm.nih.gov/nuccore/CP009167>). Use
option B<-p> to list plasmids as separate entries (lines) in the
feature table.

For draft genomes the number of contigs/scaffolds are counted. All
contigs/scaffolds of draft genomes should be marked with the I<WGS>
keyword (see e.g. draft NCBI Genbank entry
L<JSAY00000000|http://www.ncbi.nlm.nih.gov/nuccore/JSAY00000000>).
If this is not the case for your file(s) you can add those keywords
to each sequence entry with the following Perl one-liners (will edit
files in place). For files in B<GENBANK> format if 'KEYWORDS    .' is
present

C<perl -i -pe 's/^KEYWORDS(\s+)\./KEYWORDS$1WGS\./' file>

or if 'KEYWORDS' isn't present at all

C<perl -i -ne 'if(/^ACCESSION/){ print; print "KEYWORDS    WGS.\n";} else{ print;}' file>

For files in B<EMBL> format if 'KW   .' is present

C<perl -i -pe 's/^KW(\s+)\./KW$1WGS\./' file>

or if 'KW' isn't present at all

C<perl -i -ne 'if(/^DE/){ $dw=1; print;} elsif(/^XX/ && $dw){ print; $dw=0; print "KW   WGS.\n";} else{ print;}' file>

=head1 OPTIONS

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-e>=I<str>, B<-extensions>=I<str>

File extensions to include in the analysis (EMBL or GENBANK format),
either comma-separated list or multiple occurences of the option
[default = ebl,emb,embl,gb,gbf,gbff,gbank,gbk,genbank]

=item B<-p>, B<-plasmids>

Optionally list plasmids as extra entries in the feature table, if
they are annotated with a I</plasmid="plasmid_name"> tag in the
I<source> primary tag

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The resulting feature table is printed to C<STDOUT>. Redirect or
pipe into another tool as needed (e.g. C<cut>, C<grep>, or C<head>).

=back

=head1 EXAMPLES

=over

=item C<perl genomes_feature_table.pl -p -e gb,gbk E<gt> feature_table_plasmids.tsv>

=item C<perl genomes_feature_table.pl path/to/genome_dir/ -e gbf -e embl E<gt> feature_table.tsv>

=back

=head1 DEPENDENCIES

=over

=item B<L<BioPerl|http://www.bioperl.org>)>

Tested with BioPerl version 1.006923

=back

=head1 VERSION

 0.5                                               update: 14-09-2015
 0.1                                                       25-11-2011

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

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO; # BioPerl module to handle sequence input/output
use Bio::SeqFeatureI; # BioPerl module to handle features in a sequence

### Get the options with Getopt::Long
my @Extensions; # optionally give list of file-extensions to search in (otherwise defaults will be set below); comma-separated or several occurences of the option
my $Opt_Plasmids; # optionally list plasmids in multi-seq files as extra entries in final feature table
my $VERSION = 0.5;
my ($Opt_Version, $Opt_Help);
GetOptions ('extensions=s' => \@Extensions,
            'plasmids' => \$Opt_Plasmids,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
warn "\n### Warning:\nToo many arguments given, only the path to the directory with the RichSeq files (EMBL or GENBANK format) is needed!\n" if (@ARGV > 1);
my $Genome_Dir = '.'; # directory with RichSeq genome files in EMBL or GENBANK format; default current directory
$Genome_Dir = shift if (@ARGV);
if (!-d $Genome_Dir) {
    my $warning = "\n### Fatal error: Directory '$Genome_Dir' does not exist: $!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Get optional file-extensions or set defaults
if (@Extensions) {
    @Extensions = split(/,/,join(',', @Extensions)); # split (potential comma-separated lists) and join (potential several occurences) of file-extensions
} else {
    @Extensions = qw(ebl emb embl gb gbf gbff gbank gbk genbank); # default file extensions
}
my $Extensions = join('|', @Extensions); # join for regex below



### Save all primary features from all seq-files (to include all possibilities) and count them for each genome/replicon individually
my $File_Count = 0; # Count number of files
my %Strain_Features; # hash of hash; store all primary features of each strain, counted; additional genome/plasmid description (DEFINITION line), sequence length, gc_content, coding percentage, accession numbers, contigs/scaffolds (for drafts), unresolved bases ('n/N's), and plasmids (without option '-p')

opendir (my $Dir_Fh, $Genome_Dir);
while (my $file = readdir($Dir_Fh)) {
    if ($file =~ /.+\.($Extensions)$/) { # only include lines with given files extensions
        $File_Count++;
        my $Id; # filename used as unique ID (optionally with appended plasmid name) for hashes; need to declare here for checks after 'while next_seq' loop

        my $seqio_obj = Bio::SeqIO->new(-file => "<$Genome_Dir/$file");
        die "\n### Fatal error: File '$file' is not a RichSeq file in EMBL or GENBANK format!\n" if (ref($seqio_obj) !~ /Bio\:\:SeqIO\:\:[genbank|embl]/); # check if correct file format; Bio::SeqIO::genbank or Bio::SeqIO::embl object

        my %concat_seq; # store the concatenated sequences for all replicons or contigs/scaffolds of a genome with key '$Id'; only GC-content and unresolved Ns saved in '%Strain_Features', thus extra hash declared here to remove content and save memory for each file
        my %coding_bases; # base count to subsequently calculate coding percentage with key '$Id'; only coding percentage saved in '%Strain_Features', thus extra hash
        my %plasmid_count; # count the number of plasmids for current genome; hash with key as plasmid name needed for draft plasmids with SEVERAL contigs for ONE plasmid (with the same plasmid name of course) in one multi seq-file
        my $desc; # '$seq_obj->desc' in case of only-plasmid seq-file but not option '-p' use the (last) plasmid description for feature table

        my $seq_count = 0;
        while (my $seq_obj = $seqio_obj->next_seq) { # for multi-seq files with several seq entries (e.g. several replicons or contigs/scaffolds)
            $seq_count++;
            $Id = $file; # set unique $Id, based on filename; has to be set for every seq entry in case of plasmids and option '-p' (see below)
            warn "\n### Warning:\nMore than 10 sequence entries in file '$file', but no 'WGS' keyword found. Sure this is not a draft genome? Otherwise please see the help with option '-h'!\n\n" if ($seq_count == 10 && !$Strain_Features{$Id}{'contigs_scaffolds'});

            my ($plasmid) = grep { $_->primary_tag eq 'source' && $_->has_tag('plasmid') } $seq_obj->get_SeqFeatures; # check if current seq entry is a plasmid, needs '/plasmid="plasmid_name"' tag in primary tag 'source'; $plasmid contains now a primary tag 'source' feat_obj
            ($plasmid) = $plasmid->get_tag_values('plasmid') if ($plasmid); # get plasmid name; replace the 'feat_obj' in $plasmid with the tag value of '/plasmid'
            $plasmid_count{$plasmid} = 1 if ($plasmid); # count number of plasmids

            $Id .= $plasmid if ($Opt_Plasmids && $plasmid); # for plasmids with option '-p' concatenate $Id with plasmid name

            $desc = $seq_obj->desc;
            $desc =~ s/(\s*(DNA|chromosome)*, complete genome\.*|\s*(chromosome)*, complete sequence\.*|\s*complete genome, strain.*|\s*, whole genome shotgun sequence\.*|\.)$// if ($desc); # shorten strain description
            $Strain_Features{$Id}{'desc'} = $desc if ((!$Opt_Plasmids && !$plasmid) || $Opt_Plasmids); # don't use a plasmid description without '-p'; with '-p' take all descs; overwrite until last seq entry

            $Strain_Features{$Id}{'contigs_scaffolds'}++ if ($seq_obj->keywords =~ /WGS/); # count contigs/scaffolds and indicate that draft

            if (!$concat_seq{$Id}) { # first seq entry for $Id in potential multi-seq file ('!$Strain_Features{$Id}{'length'}' or '!$Strain_Features{$Id}{'first_acc'}' can also be used)
                $Strain_Features{$Id}{'first_acc'} = $seq_obj->accession_number;
                print STDERR "\nProcessing ";
                if ($Strain_Features{$Id}{'contigs_scaffolds'}) {
                    print STDERR "draft contigs/scaffolds ";
                } else {
                    print STDERR "complete replicon(s) ";
                }
                print STDERR "in '$file': $desc; $Strain_Features{$Id}{'first_acc'}\n";
                $Strain_Features{$Id}{'length'} = $seq_obj->length;
                $concat_seq{$Id} = $seq_obj->seq; # save the sequence to concatenate to potential following seq entries of same $Id in multi-seq files
                $coding_bases{$Id} = 0; # need to declare for addition below
            } else { # further contigs/scaffolds of drafts or multi-seq complete genome without option '-p'
                $Strain_Features{$Id}{'last_acc'} = $seq_obj->accession_number; # give the range of acc#s for multi-seq files; overwrite until last sequence
                $Strain_Features{$Id}{'length'} += $seq_obj->length; # add all previous lengths to the current for '$Id'
                $concat_seq{$Id} .= $seq_obj->seq; # concatenate all sequences for '$Id'
            }


            ($Strain_Features{$Id}{'gc'}, $Strain_Features{$Id}{'unresolved_n'}) = gc_content($concat_seq{$Id}); # subroutine to calculate the GC-content and count 'N's in seq; unfortunately has to be calculated every time (don't know what is last seq entry for each '$Id' or plasmid/chromosome etc. in file)


            foreach my $feat_obj ($seq_obj->get_SeqFeatures) {
                if ($feat_obj->primary_tag eq 'source') { # exclude source primary tag, has no feature annotation
                    next;
                } elsif ($feat_obj->primary_tag eq 'gene') { # count 'gene' primary tags
                    $Strain_Features{$Id}{'gene'}++;
                    $Strain_Features{$Id}{'pseudo_gene'}++ if ($feat_obj->has_tag('pseudo'));
                } elsif ($feat_obj->primary_tag eq 'CDS') {
                    $Strain_Features{$Id}{'CDS'}++;
                    if ($feat_obj->has_tag('pseudo')) {
                        $Strain_Features{$Id}{'pseudo_CDS'}++;
                    } elsif (!$feat_obj->has_tag('pseudo')) {
                        $coding_bases{$Id} = ($feat_obj->location->end - $feat_obj->location->start) + 1 + $coding_bases{$Id}; # coding_perc only calculated based on non-pseudo CDS annotation
                    }
                } else { # get ALL the other primary tags, shouldn't have 'pseudo' tags
                    $Strain_Features{$Id}{$feat_obj->primary_tag}++;
                }
            }

            $Strain_Features{$Id}{'coding_perc'} = ($coding_bases{$Id}/$Strain_Features{$Id}{'length'})*100;
            $Strain_Features{$Id}{'coding_perc'} = sprintf("%.2f", $Strain_Features{$Id}{'coding_perc'}); # round percentage to two decimals
            $Strain_Features{$Id}{'coding_perc'} = 100 if ($Strain_Features{$Id}{'coding_perc'} > 100); # > 100 can happen with very small replicons/contigs having overlapping CDSs
        }

        $Strain_Features{$Id}{'plasmids'} = keys %plasmid_count if (!$Opt_Plasmids && keys %plasmid_count > 0); # number of plasmids in the current genome, only needed without option '-p'
        $Strain_Features{$Id}{'desc'} = $desc if (!$Strain_Features{$Id}{'desc'}); # in case of plasmid-only seq-files and not option '-p'

        #die "\n### Fatal error:\nFile '$file' contains a complete genome with several replicons. However, more than one replicon (which should only be the chromosome) doesn't have a '/plasmid=\"plasmid_name\"' tag in primary tag 'source'.\nIf it is a complete genome include these tags, or if it is a draft genome include the 'WGS' keyword in ALL contigs/scaffolds as described in the help (option '-h')!\n" if ($seq_count - keys %plasmid_count > 1 && !$Strain_Features{$Id}{'contigs_scaffolds'}); # too strict? Also, overlaps with warn 'Warning:\nMore than 10 sequence' above!
        #die "\n### Fatal error:\nFile '$file' contains a draft genome with only one contig.\nIf it is a complete genome remove the 'WGS' keyword from the only present replicon (see '-h')!\n" if ($seq_count == 1 && $Strain_Features{$Id}{'contigs_scaffolds'}); # too strict? e.g. uncircularized single-contig draft seq files
    }
}
closedir $Dir_Fh;
die "\n### Fatal error: No EMBL or GENBANK format files found!\n" if (!$File_Count);



### Get all annotated primary features present in all files
my %Primary_Features;
foreach my $id (keys %Strain_Features) {
    foreach (keys %{$Strain_Features{$id}}) { # de-references hash of hash
        $Primary_Features{$_} = 1;
    }
}



### Print header of output with all existent primary features
print "# name (last contig/scaffold for drafts)\tsize [bp]\tGC-content [%]\tcoding percentage [%]\tCDS (pseudo)\tgenes (pseudo)\trRNA\ttRNA\ttmRNA\tncRNA\taccession (first..last)\tcontigs/scaffolds\tunresolved bases (Ns)\t"; # specific print order for the most common primary features, '#' indicates comment line
my $Common_Primary_Feats = 'desc|length|gc|coding_perc|CDS|pseudo_CDS|gene|pseudo_gene|rRNA|tRNA|tmRNA|ncRNA|first_acc|last_acc|contigs_scaffolds|unresolved_n'; # the common primaries for the regexs and 'split' for print out below

my $Common_Primary_Count = grep(/^($Common_Primary_Feats)$/, keys %Primary_Features); # count occurence of common primary features (not all present in all seq entries) for tab or newline print below
my $i = 0;
foreach (sort keys %Primary_Features) { # print the residual primary features
    if ($_ =~ /^($Common_Primary_Feats)$/) { # exclude the print order ones
        next;
    } else {
        $i++;
        print "$_";
        print "\t" if ($i < (keys %Primary_Features) - $Common_Primary_Count); # don't print tab for last element (but newline below)
    }
}
print "\n";

print "# NA = tag, key, or qualifier not existent\n"; # comment line



### Print the primary features for each strain
foreach my $id (sort {$Strain_Features{$a}{'desc'} cmp $Strain_Features{$b}{'desc'}} keys %Strain_Features) {
    foreach (split(/\|/, $Common_Primary_Feats)) { # split '$Common_Primary_Feats' string to print according to print order above
        print " (" if (/^pseudo/); # 'pseudo_CDS' and 'pseudo_gene' are printed in parantheses
        if ($Strain_Features{$id}{$_}) { # print value only if defined
            print ".." if (/last_acc/); # for 'first_acc..last_acc'
            if ($Strain_Features{$id}{$_} =~ /unknown/) { # BioPerl returns 'unknown' for non-existent acc-nr
                print "NA";
            } else {
                print "$Strain_Features{$id}{$_}";
            }
        } else { # feature not defined for this '$Id'
            next if ($_ =~ /last_acc/); # skip 'last_acc' if not defined (otherwise will print 'NA' for it)
            print "NA";
        }
        print ")" if (/^pseudo/); # 'pseudo_CDS' and 'pseudo_gene' are printed in parantheses
        print "\t" unless (($_ =~ /^(CDS|gene)$/) || ($_ =~ /first_acc/ && $Strain_Features{$id}{'last_acc'})); # CDS and gene are followed by parantheses for possible pseudos, 'first_acc' might be followed by '..last_acc'; thus no tab for these
    }

    $i = 0;
    foreach (sort keys %Primary_Features) { # print the residual primary features
        if ($_ =~ /^($Common_Primary_Feats)$/) { # exclude the above ones
            next;
        } elsif ($Strain_Features{$id}{$_}) {
            $i++;
            print "$Strain_Features{$id}{$_}";
        } else {
            $i++;
            print "NA";
        }
        print "\t" if ($i < (keys %Primary_Features) - $Common_Primary_Count); # don't print tab for last element
    }

    print "\n"; # newline for each line in the feature table
}

print STDERR "\n$File_Count RichSeq file(s) was/were processed!\n";

exit;



###############
# Subroutines #
###############

### Subroutine to calculate GC-content
sub gc_content {
    my $seq = shift;
    my $a = ($seq =~ tr/[aA]//); # transliterations don't accept modifiers like case-insensitive 'i'
    my $c = ($seq =~ tr/[cC]//);
    my $g = ($seq =~ tr/[gG]//);
    my $t = ($seq =~ tr/[tT]//);
    my $n = ($seq =~ tr/[nN]//);
    my $gc_content = (($c + $g)/($a + $c + $g + $t))*100;
    $gc_content = sprintf("%.2f", $gc_content); # round percentage to two decimals
    return ($gc_content, $n);
}
