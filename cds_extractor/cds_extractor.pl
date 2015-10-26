#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<cds_extractor.pl> - extract protein or DNA sequences from CDS features

=head1 SYNOPSIS

C<perl cds_extractor.pl -i seq_file.[embl|gbk] -p>

=head1 DESCRIPTION

Extracts protein or DNA sequences of CDS features from a
(multi)-RichSeq file (e.g. EMBL or GENBANK format) and writes them to a
multi-FASTA file. The FASTA headers for each CDS include either the
locus tag, if that's not available, protein ID, gene, or an internal
CDS counter as identifier (in this order). The organism info
includes also possible plasmid names. Pseudogenes (tagged by
B</pseudo>) are not included (except in the CDS counter).

In addition to the identifier, FASTA headers include gene (B<g=>),
product (B<p=>), organism (B<o=>), and EC numbers (B<ec=>), if these
are present for a CDS. Individual EC numbers are separated by
B<semicolons>. The location/position (B<l=>start..stop) of a CDS will
always be included. If gene is used as FASTA header ID
'B<g=>gene' will only be included with option B<-f>.

Fuzzy locations in the feature table of a sequence file are not
taken into consideration for B<l=>. If you set options B<-u>
and/or B<-d> and the feature location overlaps a B<circular>
replicon boundary, positions are marked with '<' or '>' in the
direction of the exceeded boundary. Features with overlapping
locations in B<linear> sequences (e.g. contigs) will be skipped and
are B<not> included in the output! A CDS feature is on the lagging
strand if start > stop in the location. In the special case of
overlapping circular sequence boundaries this is reversed.

Of course, the B<l=> positions are separate for each sequence in a
multi-sequence file. Thus, if you want continuous positions for the
CDSs run these files first through L<C<cat_seq.pl>|/cat_seq>.

Optionally, a file with locus tags can be given to extract only
these CDS features with option B<-l> (each locus tag in a new line).

=head1 OPTIONS

=head2 Mandatory options

=over 23

=item B<-i>=I<str>, B<-input>=I<str>

Input RichSeq sequence file including CDS annotation (e.g. EMBL or
GENBANK format)

=item B<-p>, B<-protein>

Extract B<protein> sequence for each CDS feature, excludes option B<-n>

B<or>

=item B<-n>, B<-nucleotide>

Extract B<nucleotide> sequence for each CDS feature, excludes option
B<-p>

=back

=head2 Optional options

=over 28

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-u>=I<int>, B<-upstream>=I<int>

Include given number of flanking nucleotides upstream of each CDS
feature, forces option B<-n>

=item B<-d>=I<int>, B<-downstream>=I<int>

Include given number of flanking nucleotides downstream of each CDS
feature, forces option B<-n>

=item B<-c>=I<str>, B<-cds_prefix>=I<str>

Prefix for the internal CDS counter [default = 'CDS']

=item B<-l>=I<str>, B<-locustag_list>=I<str>

List of locus tags to extract only those (each locus tag on a new line)

=item B<-f>, B<-full_header>

If gene is used as ID include additionally 'B<g=>gene' in FASTA
headers, so downstream analyses can recognize the gene tag (e.g.
L<C<prot_finder.pl>|/prot_finder>).

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 23

=item F<*.faa>

Multi-FASTA file of CDS protein sequences

B<or>

=item F<*.ffn>

Multi-FASTA file of CDS DNA sequences

=item (F<no_annotation_err.txt>)

Lists input files missing CDS annotation, script exited with B<fatal error> i.e. no FASTA output file

=item (F<double_id_err.txt>)

Lists input files with ambiguous FASTA IDs, script exited with B<fatal error> i.e. no FASTA output file

=item (F<locus_tag_missing_err.txt>)

Lists CDS features without locus tags

=item (F<linear_seq_cds_overlap_err.txt>)

Lists CDS features overlapping sequence border of a B<linear>
molecule, which are not included in the result multi-FASTA file

=back

=head1 EXAMPLES

=over

=item C<perl cds_extractor.pl -i seq_file.gbk -p -l locus_tags.txt>

=item C<perl cds_extractor.pl -i seq_file.embl -n -l locus_tags.txt -u 100 -d 20>

=item C<perl cds_extractor.pl -i Ecoli_MG1655.gbk -p -f -c MG1655>

=back

=head1 DEPENDENCIES

=over

=item B<BioPerl (L<http://www.bioperl.org>)>

Tested with BioPerl version 1.006923

=back

=head1 VERSION

 0.7.1                                             update: 26-10-2015
 0.1                                                       24-05-2012

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
use Bio::SeqIO; # bioperl module to handle sequence input/output
# use Bio::Seq; # bioperl module to play with the sequence and its features ### apparently not needed, methods inherited
# use Bio::SeqFeatureI; # bioperl module to handle features in a sequence ### apparently not needed, methods inherited



### Get options with Getopt::Long, works also abbreviated and with two "--": -i, --i, -input ...
my $Seq_File; # RichSeq sequence file including feature annotation
my $Opt_Protein; # extract protein sequences for each CDS feature; excludes option '-n'
my $Opt_Nucleotide; # extract nucleotide sequences for each CDS feature; excludes option '-p'
my $Upstream = 0; # include given number of flanking nucleotides upstream of each CDS feature; forces option '-n'
my $Downstream = 0; # include given number of flanking nucleotides downstream of each CDS feature; forces option '-n'
my $Prefix; # prefix for the internal CDS counter
my $Locustag_List; # list of locus_tags to extract only those
my $Opt_Full_Header; # include a full FASTA header for downstream 'prot_finder.pl' analysis (needed to identify /gene feature tag if used as identifier)
my $VERSION = '0.7.1';
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Seq_File,
            'protein' => \$Opt_Protein,
            'nucleotide' => \$Opt_Nucleotide,
            'upstream:i' => \$Upstream,
            'downstream:i' => \$Downstream,
            'cds_prefix:s' => \$Prefix,
            'locustag_list:s' => \$Locustag_List,
            'full_header' => \$Opt_Full_Header,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Seq_File) {
    my $warning = "\n### Fatal error: Option '-i' or its argument is missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}



### Enforce mandatory or optional options
if (($Upstream || $Downstream) && !$Opt_Nucleotide) {
    warn "Option '-d' and/or '-u' set, but not '-n'. Forcing option '-n'!\n";
    $Opt_Nucleotide = 1;
    undef $Opt_Protein;
}
if (!$Opt_Protein && !$Opt_Nucleotide) {
    die "\n### Fatal error: None of the mandatory options ('-p' or '-n') given! Please choose one of the extraction methods!\n";
} elsif ($Opt_Protein && $Opt_Nucleotide) {
    die "\n### Fatal error: Both mandatory options ('-p' and '-n') given! Choose only one of the extraction methods!\n";
}
$Prefix = 'CDS' if (!$Prefix); # default for internal CDS



### Create a bioperl Bio::SeqIO object with $Seq_File
my $seqio_object = Bio::SeqIO->new(-file => "<$Seq_File"); # no '-format' to leave to bioperl guessing



### Print the input and output file which are processed and written, respectively
print "\nInput: $Seq_File\t";
my $Ori_Filename = $Seq_File; # store original filename to give error statements
$Seq_File =~ s/(.+)\.\w+$/$1/; # strip filename extension
my $Out_Fh; # filehandle for output file
if ($Opt_Protein) {
    $Seq_File = $Seq_File.'.faa';
    print "Output: $Seq_File\n";
    open ($Out_Fh, ">", "$Seq_File");
} elsif ($Opt_Nucleotide) {
    $Seq_File = $Seq_File.'.ffn';
    print "Output: $Seq_File\n";
    open ($Out_Fh, ">", "$Seq_File");
}



### Get the list of locus_tags if $Locustag_List defined
my %Locus_Tags; # store locus_tags and indicate if found in $Seq_File
if ($Locustag_List) {
    open (my $locus_fh, "<", "$Locustag_List");
    while (<$locus_fh>) {
        chomp;
        next if (/^$/); # skip empty lines
        $Locus_Tags{$_} = 0; # changes to 1 if the locus tag was found in $Seq_File, see below
    }
    close $locus_fh;
}



### Write output FASTA with respective header lines and either protein seq from /translation feature tag or nucleic acid subseq
my $Anno_Present = 0; # switches to true if CDS primary features present
my $Organism; # store /organism tag value to include in each header, include plasmid name
my $Cds_Count = sprintf("%04d", 0); # internal CDS counter, counts also 'pseudo' CDSs
my %Double_Id; # control if a locus_tag or protein_id is ambiguous in $Seq_File and die (they should be unique); see subroutine 'control_double'

my $No_Locus_Tag = 0; # count how many CDS features don't have a locus_tag
my $Locus_Tag_Miss_Err = 'locus_tag_missing_err.txt'; # error file that contains all alternative IDs (see subroutine 'locus_tag_missing')
my $Locus_Tag_Miss_Err_Fh; # filehandle

my $Linear_Seq_Overlap = 0; # count CDS features which overlap linear sequence borders (see subroutine 'linear_overlap')
my $Linear_Overlap_Err = 'linear_seq_cds_overlap_err.txt'; # error file that contains all affected IDs
my $Linear_Seq_Overlap_Fh; # filehandle

while (my $seq_object = $seqio_object->next_seq) { # a Bio::Seq object
    my @feat_objects = $seq_object->get_SeqFeatures; # slurp all features from the seq to check CDS primary feature annotation

    if (grep ($_->primary_tag eq 'CDS', @feat_objects) == 0) { # scalar context
        next; # if no CDS features present skip to next $seq_object in (multi-)seq file
    }
    $Anno_Present = 1; # if not skipped CDS annotation present

    $Organism = ''; # empty $Organism for each $seq_object
    foreach my $feat_object (@feat_objects) { # a Bio::SeqFeatureI object
        if ($Locustag_List) { # skip residual FEATURE objects if all locus_tags in $Locustag_List found
            if (grep ($Locus_Tags{$_} == 1, keys %Locus_Tags) == keys %Locus_Tags) { # scalar context
                last;
            }
        }

        if ($feat_object->primary_tag eq 'source') {
            $Organism = feature_value_eval($feat_object, 'organism'); # subroutine to evaluate existence of the tag and return value
            if ($feat_object->has_tag('plasmid')) { # subroutine 'feature_value_eval' also possible, but method 'has_tag' more elegantly
                my ($plasmid) = $feat_object->get_tag_values('plasmid'); # values always returned as ARRAYS
                $plasmid =~ s/\s/_/g;
                $Organism = $Organism.'-plasmid_'.$plasmid;
            }
        }

        if ($feat_object->primary_tag eq 'CDS') {
            $Cds_Count++;
            if ($feat_object->has_tag('pseudo')) { # skip pseudogenes, they don't include '/translation'
                next;
            }
            my ($start, $stop, $strand, $seq_stop, $location_stop) = get_location($feat_object, $seq_object); # subroutine to get start, stop and strand of feature; $seq_stop needed for sub 'print_seq' and with $location_stop needed for sub 'print_location'

            if ($feat_object->has_tag('locus_tag')) { # LOCUS_TAG-ELSIF
                my ($locus_tag) = $feat_object->get_tag_values('locus_tag');
                control_double($locus_tag, 'locus_tag'); # subroutine to control uniqueness of identifier (here value of /locus_tag)

                if ($Locustag_List) { # if '-locustag_list' is set only get those CDSs
                    my ($locus) = grep (/$locus_tag/i, keys %Locus_Tags); # case-insensitive for typing mistakes
                    if ($locus) {
                        $Locus_Tags{$locus} = 1;
                        if ($Opt_Nucleotide && ($start eq 'not_circular' || ($location_stop && $location_stop eq 'not_circular'))) { # CDS feature overlaps seq border of linear molecule
                            linear_overlap($locus_tag); # subroutine to write overlap CDSs to '$Linear_Overlap_Err'
                            next; # jump to the next feature object
                        }
                        print_fasta_ID($feat_object, $locus_tag); # subroutine to print the identifier (here /locus_tag), /gene (g=) and /product (p=) (both if present) to the FASTA header
                        print_location($start, $stop, $strand, $seq_stop, $location_stop); # subroutine to print feature location/position (l=) to FASTA header
                        print_org_ec($feat_object); # subroutine to print /organism (o=) and /EC_numbers (ec=) (both if present) to FASTA header
                        print_seq($feat_object, $start, $stop, $strand, $seq_stop, $seq_object); # subroutine to print the protein or nucleic sequence
                    }
                    next; # jump to the next feature object
                }

                ### '$Locustag_List' not defined
                if ($Opt_Nucleotide && ($start eq 'not_circular' || ($location_stop && $location_stop eq 'not_circular'))) {
                    linear_overlap($locus_tag); # subroutine
                    next; # jump to the next feature object
                }
                print_fasta_ID($feat_object, $locus_tag); # subroutine

            } elsif ($Locustag_List) { # in case a list of locus_tags is given, no need to look at CDSs that don't have a /locus_tag
                next; # jump to the next feature object

            } elsif ($feat_object->has_tag('protein_id')) { # PROTEIN_ID-ELSIF
                my ($protein_id) = $feat_object->get_tag_values('protein_id');
                control_double($protein_id, 'protein_id'); # subroutine
                locus_tag_missing($protein_id, 'protein_id'); # subroutine to inform no /locus_tag is present for this CDS
                $No_Locus_Tag++;
                if ($Opt_Nucleotide && ($start eq 'not_circular' || ($location_stop && $location_stop eq 'not_circular'))) {
                    linear_overlap($protein_id); # subroutine
                    next;
                }
                print_fasta_ID($feat_object, $protein_id); # subroutine

            } elsif ($feat_object->has_tag('gene')) { # GENE-ELSIF
                my ($gene) = $feat_object->get_tag_values('gene');
                control_double($gene, 'gene'); # subroutine
                locus_tag_missing($gene, 'gene'); # subroutine
                $No_Locus_Tag++;
                if ($Opt_Nucleotide && ($start eq 'not_circular' || ($location_stop && $location_stop eq 'not_circular'))) {
                    linear_overlap($gene); # subroutine
                    next;
                }
                print_fasta_ID($feat_object, ''); # subroutine; give empty string for $tag to sub, to be able to determine that the sub call is from GENE-ELSIF and 'g=' should only be included if option '-f' is set

            } else { # if none of the above tags are existent use the internal CDS counter; CDS_COUNT-ELSE
                my $cds_prefix_count = $Prefix.'_'.$Cds_Count;
                locus_tag_missing($cds_prefix_count, 'cds-counter'); # subroutine
                $No_Locus_Tag++;
                if ($Opt_Nucleotide && ($start eq 'not_circular' || ($location_stop && $location_stop eq 'not_circular'))) {
                    linear_overlap($cds_prefix_count); # subroutine
                    next;
                }
                print_fasta_ID($feat_object, $cds_prefix_count); # subroutine
            }

            ### Print location, organism, EC_number(s), and sequence for all $feat_object if '$Locustag_List' not defined
            print_location($start, $stop, $strand, $seq_stop, $location_stop); # subroutine
            print_org_ec($feat_object); # subroutine
            print_seq($feat_object, $start, $stop, $strand, $seq_stop, $seq_object); # subroutine
        }
    }

    if ($Locustag_List) { # skip ALSO residual SEQUENCE objects if all locus_tags are found
        if (grep ($Locus_Tags{$_} == 1, keys %Locus_Tags) == keys %Locus_Tags) {
            last;
        }
    }
}

close $Out_Fh;
if ($No_Locus_Tag > 0) {
    warn "### $No_Locus_Tag CDS feature(s) don't have a locus tag in file '$Ori_Filename'. The respective ID(s) are written to the error file '$Locus_Tag_Miss_Err'!\n";
    close $Locus_Tag_Miss_Err_Fh;
}
if ($Linear_Seq_Overlap > 0) {
    warn "### $Linear_Seq_Overlap CDS feature(s) overlap a sequence border of a linear molecule in '$Ori_Filename' because of your settings for '-u' and/or '-d'. These CDS feature(s) are not included in the output file '$Seq_File'! The respective ID(s) are written to the error file '$Linear_Overlap_Err'!\n";
    close $Linear_Seq_Overlap_Fh;
}



### Exit with error and remove empty result file if no CDS features/annotation found in (multi-)seq file
if (!$Anno_Present) {
    my $anno_err = 'no_annotation_err.txt';
    err_file_exist($anno_err); # subroutine to test for file existence and give warning to STDERR
    open (my $anno_err_fh, ">>", "$anno_err");
    print $anno_err_fh "$Ori_Filename\n";
    close $anno_err_fh;
    unlink $Seq_File;
    die "\n###Fatal error\nNo CDS annotation in '$Ori_Filename', the respective filename was written to '$anno_err'!\nExiting program!\n\n";
}



### Print locus tags that were not found in seq_file
if ($Locustag_List) {
    my @missed = grep ($Locus_Tags{$_} == 0, keys %Locus_Tags);
    if (@missed) {
        print "### The following locus tags were not found in '$Ori_Filename':\n";
        foreach (sort @missed) {
            print "$_\t";
        }
        print "\n";
    }
}


exit;


###############
# Subroutines #
###############

### Control if an identifier for the FASTA header is ambiguous in the $Seq_File and die (they should be unique/unambiguous)
sub control_double {
    my ($tag, $type) = @_;
    if ($Double_Id{$tag}) {
        my $double_id_err = 'double_id_err.txt';
        my $file_exist = err_file_exist($double_id_err); # subroutine
        open (my $double_id_err_fh, ">>", "$double_id_err");
        print $double_id_err_fh "Original_filename\tID-type\tDouble_ID\n" if (!$file_exist); # print header in file only if file doesn't exist
        print $double_id_err_fh "$Ori_Filename\t$type\t$tag\n";
        close $double_id_err_fh;
        unlink $Seq_File;
        die "\n###Fatal error!\n'$tag' of type '$type' exists at least two times in organism '$Organism' of file '$Ori_Filename', but should be unambiguous for downstream analyses! The error is written to the error file '$double_id_err'. Please modify all repetitive occurences.\nExiting program!\n\n";
    } else {
        $Double_Id{$tag} = 1;
    }
    return 1;
}



### Test for error file existence and give warning to STDERR
sub err_file_exist {
    my $file = shift;
    if (-e $file) {
        warn "\nThe error file '$file' exists already, the current errors will be appended to the existing file!\n";
        return 1;
    }
    return 0;
}



### Return value of a tag (replace whitespaces with '_') or catch error if non-existent to return empty value
sub feature_value_eval {
    my ($feat_object, $tag) = @_;
    my $value = '';
    eval {$value = join(';', $feat_object->get_tag_values($tag));}; # values always returned as ARRAYS; catch error if tag doesn't exist. To seperate /EC_numbers ';' is used, the only feature tag needed with more than one occurrence in a single CDS
    $value =~ s/\s/_/g; # replace spaces with '_' (guess only needed for /organism and /product values)
    return $value;
}



### Get the position/location of a feature
sub get_location {
    my ($feat_object, $seq_object) = @_;
    my ($start, $stop, $strand, $seq_stop, $location_stop); # $seq_stop needed for sub 'print_seq' and with $location_stop needed for sub 'print_location'

    $strand = $feat_object->strand; # 1 = feature on leading strand, -1 = feature on lagging strand
    if ($strand == 1) {
        $start = $feat_object->start - $Upstream;
        $stop =  $feat_object->end + $Downstream;
    } elsif ($strand == -1) {
        $start = $feat_object->start - $Downstream;
        $stop = $feat_object->end + $Upstream;
    }

    ### Adjust positions if they overlap replicon sequence borders
    if ($start < 0) { # '$seq_obj->subseq|->trunc' (in subroutine 'print_seq') don't work with negative numbers, but with positions > '$seq_obj->length'
        if ($seq_object->is_circular) { # true if molecule is circular
            $start = $seq_object->length - abs($start);
        } else { # if sequence linear wrong to print overlapping sequences (see sub 'linear_overlap')
            $start = 'not_circular';
        }
        $seq_stop = $seq_object->length + $stop;
    }
    if ($stop > $seq_object->length) { # for sub 'print_location' a stop > '$seq_obj->length' is no use and should start again from the seq start
        if ($seq_object->is_circular) {
            $location_stop = $stop - $seq_object->length;
        } else {
            $location_stop = 'not_circular';
        }
    }

    return ($start, $stop, $strand, $seq_stop, $location_stop);
}



### Inform CDS locations with options '-u' and/or '-d' are overlapping a sequence border of a linear molecule
sub linear_overlap {
    my $tag = shift;
    if ($Linear_Seq_Overlap == 0) { # open error file to append errors only for first overlap error
        my $file_exist = err_file_exist($Linear_Overlap_Err); # subroutine
        open ($Linear_Seq_Overlap_Fh, ">>", "$Linear_Overlap_Err");
        print $Linear_Seq_Overlap_Fh "Original_filename\tOrganism\tID\n" if (!$file_exist); # header of error file
    }
    $Linear_Seq_Overlap++;
    print $Linear_Seq_Overlap_Fh "$Ori_Filename\t$Organism\t$tag\n";
    return 1;
}



### Inform if locus_tags are missing
sub locus_tag_missing {
    my ($tag, $type) = @_;
    if ($No_Locus_Tag == 0) { # open error file to append errors only for first missing locus_tag occurrence
        my $file_exist = err_file_exist($Locus_Tag_Miss_Err); # subroutine
        open ($Locus_Tag_Miss_Err_Fh, ">>", "$Locus_Tag_Miss_Err");
        print $Locus_Tag_Miss_Err_Fh "Original_filename\tOrganism\tAlternative_ID_type\tAlternative_ID\n" if (!$file_exist); # header of error file
    }
    print $Locus_Tag_Miss_Err_Fh "$Ori_Filename\t$Organism\t$type\t$tag\n";
    return 1;
}



### Print identifier, /gene (g=) and /product (p=) (both if existent) of a CDS to FASTA header
sub print_fasta_ID {
    my ($feat_object, $tag) = @_;
    my $gene = feature_value_eval($feat_object, 'gene'); # subroutine
    my $product = feature_value_eval($feat_object, 'product'); # subroutine
    if ($tag) { # $tag is defined for LOCUS_TAG-ELSIF, PROTEIN_ID-ELSIF, and CDS_COUNT-ELSE
        print $Out_Fh ">$tag";
    } elsif (!$tag) { # subroutine call from GENE-ELSIF with empty string for $tag (see below 'g=')
        print $Out_Fh ">$gene";
    }
    print $Out_Fh " g=$gene" if (($gene && $tag) || ($Opt_Full_Header && !$tag)); # print 'g='-gene only if /gene tag present; additionally, only print 'g=' for GENE-ELSIF if '-f' is set, otherwise can't determine if the identifier is a gene in downstream tools (e.g. 'prot_finder.pl')
    print $Out_Fh " p=$product" if ($product);
    return 1;
}



### Print position/location (l=) of a CDS to FASTA header
sub print_location {
    my ($start, $stop, $strand, $seq_stop, $location_stop) = @_;
    $stop = $location_stop if ($location_stop); # see sub 'get_location'
    print $Out_Fh " l=";
    if ($strand == 1) {
        print $Out_Fh ">" if ($seq_stop); # defined if feature start position overlaps replicon's seq start (with '-u' or '-d'), sub 'get_location'
        print $Out_Fh "$start..$stop";
        print $Out_Fh "<" if ($location_stop); # defined if feature stop position overlaps replicon's seq end (with '-u' or '-d')
        return 1;
    } elsif ($strand == -1) {
        print $Out_Fh "<" if ($location_stop); # switched on lagging strand, because $stop and $start switched
        print $Out_Fh "$stop..$start"; # switch $stop and $start to indicate position on lagging strand
        print $Out_Fh ">" if ($seq_stop);
        return -1;
    }
    return 0;
}



### Print organism (o=) and EC numbers (ec=) (both if existent) of a CDS to FASTA header
sub print_org_ec {
    my $feat_object = shift;
    print $Out_Fh " o=$Organism" if ($Organism);
    my $ec_numbers = feature_value_eval($feat_object, 'EC_number'); # subroutine
    print $Out_Fh " ec=$ec_numbers" if ($ec_numbers);
    print $Out_Fh "\n";
    return 1;
}



### Print the protein or nucleic sequence to the result file
sub print_seq {
    my ($feat_object, $start, $stop, $strand, $seq_stop, $seq_object) = @_;
    if ($Opt_Protein) {
        print $Out_Fh $feat_object->get_tag_values('translation'), "\n";

    } elsif ($Opt_Nucleotide) {
        $stop = $seq_stop if ($seq_stop); # see sub 'get_location'
        if ($strand == 1) {
            my $subseq = $seq_object->subseq($start, $stop);
            print $Out_Fh substr($subseq, 0, $Upstream); # upstream bases, not possible to include in '$feat_obj->spliced_seq'
            print $Out_Fh $feat_object->spliced_seq->seq; # to include also features with "join"/fuzzy locations; '$feat_object->spliced_seq' returns a Bio::Seq object thus call funtion '->seq' for the seq string
            print $Out_Fh substr($subseq, length($subseq) - $Downstream, $Downstream), "\n"; # downstream bases

        } elsif ($strand == -1) {
            my $trunc_obj = $seq_object->trunc($start, $stop); # a Bio::Seq object, needed for revcom
            my $rev_obj = $trunc_obj->revcom; # a Bio::Seq object
            print $Out_Fh substr($rev_obj->seq, 0, $Upstream);
            print $Out_Fh $feat_object->spliced_seq->seq;
            print $Out_Fh substr($rev_obj->seq, $rev_obj->length - $Downstream, $Downstream), "\n";
        }
    }
    return 1;
}
