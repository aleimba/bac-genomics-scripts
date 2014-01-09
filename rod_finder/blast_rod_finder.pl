#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

blast_rod_finder.pl                                       07-11-2011

=head1 SYNOPSIS

C<perl blast_rod_finder.pl -q query.embl -r blastn.out -m 2000>

=head1 DESCRIPTION

This script is intended to identify region of differences (RODs)
between a nucleotide query and a nucleotide subject/reference sequence.
In order to do so, a BLASTN (L<http://blast.ncbi.nlm.nih.gov/Blast.cgi>)
needs to be performed beforehand with the query and the subject sequences.
I<blast_rod_finder.pl> is mainly designed to work with bacterial genomes,
while a query genome can be blasted against several subject sequences to
detect RODs over a number of references. Although the results are
optimized towards a complete query genome, both the reference(s) as well
as the query can be used in draft form. To create artificial genomes use
I<cat_seq.pl> or the EMBOSS application union (L<http://emboss.sourceforge.net/>).

The BLASTN report file, the query sequence file (preferably in embl or
genbank format with annotation, see below) and a minimum size for ROD
detection have to be provided for I<blast_rod_finder.pl> to run.
Subsequently, RODs are summarized in a tab-separated result file, a gff3
(usable e.g. in Artemis/DNAPlotter,
L<http://www.sanger.ac.uk/resources/software/artemis/>) and a BRIG
(BLAST Ring Image Generator, L<http://brig.sourceforge.net/>) output file.
Nucleotide sequences of each ROD are written to a multi-fasta file.

The query sequence can be provided in embl or genbank format, but has to
correspond to the fasta file used in querying the BLAST database (the
accession numbers have to correspond to the fasta headers). Use
I<seq_format-converter.pl> to create a corresponding fasta file from 
embl|genbank files for BLASTN if needed. With annotated query files
additional info is given in the result summary and the amino acid
sequences of all non-pseudo CDSs, which are contained or overlap a ROD,
are written to a result file. Furthermore, all detected RODs are saved in
individual sequence files in the corresponding query sequence format.

Run BLASTN and the script I<blast_rod_finder.pl> manually or use the bash
shell scripts I<blast_rod_finder*.sh> (see examples below) to perform the
pipeline consecutively in one folder.

The Perl script runs on BioPerl (L<http://www.bioperl.org>).

=head1 OPTIONS

=head2 Mandatory options

=over 21

=item B<-m>=I<int>, B<-min>=I<int>     Minimum size of RODs that are reported

=item B<-q>=I<str>, B<-query>=I<str>

Query sequence file [fasta, embl, or genbank format]

=item B<-r>=I<str>, B<-report>=I<str>  BLASTN report/output file

=back

=head2 Optional options

=over

=item B<-h>, B<-help>            Help (perldoc POD)

=back

=head1 OUTPUT

=over 17

=item F<./results>

All output files are stored in this result folder

=item F<rod_summary.txt>

Summary of detected ROD regions (for embl/genbank queries includes annotation), tab-separated

=item F<rod.gff>

GFF3 file with ROD coordinates to use in Artemis/DNAPlotter etc.

=item F<rod_BRIG.txt>

ROD coordinates to use in BRIG (BLAST Ring Image Generator), tab-separated

=item F<rod_seq.fasta> 

Nucleotide sequences of ROD regions (>ROD# size start..stop), multi-fasta

=item (F<rod_aa_fasta.txt>)

Only present if query is given with annotation, i.e. embl|genbank format

Amino acid sequences of all CDSs that are contained in or overlap a ROD region
in multi-fasta format (>locus_tag gene product). RODs are seperated in the file via
'~~' (~~ROD# size start..stop).

=item (F<ROD#.embl|gbk>)

Only present if query is given with annotation, i.e. embl|genbank format.

Each identified ROD is written to an individual sequence file (in the same format
as the query).

=back

=head1 EXAMPLES

=head2 Legacy BLASTN

=over

=item C<formatdb -p F -i subject.fasta -n blast_db>

=item C<blastall -p blastn -d blast_db -i query.fasta -o blastn.out -e 2e-11 -F F>

=back

=head2 BLASTN+

=over

=item C<makeblastdb -in subject.fasta -input_type fasta -dbtype nucl -out blast_db>

=item C<blastn -db blast_db -query query.fasta -out blastn.out -evalue 2e-11 -soft_masking false>

=back

=head2 blast_rod_finder.pl

=over

=item C<perl blast_rod_finder.pl -q query.embl|gbk|fasta -r blastn.out -m 5000>

=back

=head2 All-in-one with unix bash-shell scripts

=over

=item C<./blast_rod_finder_legacy.sh subject.fasta query.fasta query.embl|gbk|fasta 5000>

=item C<./blast_rod_finder_plus.sh subject.fasta query.fasta query.embl|gbk|fasta 5000>

=back

=head1 VERSION

0.4                                                update: 13-02-2013

=head1 AUTHORS

Andreas Leimbach        aleimba[at]gmx[dot]de
David Studholme         D[dot]J[dot]Studholme[at]exeter[dot]ac[dot]uk

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
use Getopt::Long; # module to get options from the command line
use Bio::SearchIO; # bioperl module to handle blast reports
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqFeatureI; # bioperl module to handle features in a sequence
use Bio::SeqUtils; # bioperl module with additional methods (including features) for Bio::Seq objects (e.g. revcom, truncate, concatenate)


### Get the options with Getopt::Long, works also abbreviated and with two "--": -r, --r, -report ...
my $usage = "perldoc $0";
my $blast_report = ''; # blastn report/output file
my $q_file = ''; # query sequence file to get the RODs and possible including annotations
my $minimum_size = ''; # minimum size of uncovered regions (RODs) that are reported
my $help = ''; # run perldoc on the POD
GetOptions ('report=s' => \$blast_report, 'query=s' => \$q_file, 'min=s' => \$minimum_size, 'help' => \$help);


### Run perldoc on POD if -h|--h|--help is given as argument or arguments are not given
if (!$blast_report || !$q_file || !$minimum_size) {
    die system($usage);
} elsif ($help) {
    die system($usage);
}


### Parse the blast report/output file
my %q_covered; # stores for each seq position (key) of a query accession if it falls within a hsp hit (value set to 1)
my @blast_q_acc; # store query accessions to test if identical in blast-parse and query-seq (seqio below)
my $blast_db; # add to GFF3 print out
my $parser = new Bio::SearchIO(-file => "<$blast_report", -format => 'blast'); # Bio::SearchIO object
$|++; # turn on autoflush, forces STDOUT to flush right away
print "\nParsing blast report file."; # status display (with autoflush)
while(my $result = $parser->next_result) { # several query sequences possible (result = entire analysis for a single query seq)
    my $query_acc = $result->query_accession;
    $query_acc =~ s/\.$//; # rm a '.' if present at the end of the string (for non-NCBI fasta headers)
    $query_acc =~ s/\.\d$//; # rm version nr. from NCBI query acc (to fit it to the acc.-nr. from Bio::SeqIO below)
    push(@blast_q_acc, $query_acc);
    $blast_db = $result->database_name;
    chop $blast_db; # always an empty character after the name?
    while (my $hit = $result->next_hit) { # several subject sequences in the database might have hits to a single query
	print '.'; # status display
	while(my $hsp = $hit->next_hsp) { # each hit might have one or more hsps (the alignments shown in a blast report)
	    my ($query_start, $query_end) = $hsp->range('query'); # query hsp start/stop coords
	    foreach my $i ($query_start .. $query_end) {
		$q_covered{$query_acc}{$i}++; # changes for each hsp hit the value from 'undef' to '1'
	    }
	}
    }
}


### Read the query sequence file into RAM, for Bio::Seq::RichSeq files (e.g. embl/genbank) also read features
my %q_seqobj; # combine all seqobj (values; for multi-fasta/embl/gbk queries) to acc.-nr.s (keys )
my %features; # stores all features of a query RichSeq file
my $multi = 0; # test if multi-seq query (multi-fasta/embl/genbank) >= 2; and counter
my $seqio_obj = Bio::SeqIO->new(-file => "<$q_file"); # no '-format' to leave it to bioperl guessing
print "\nParsing query sequence file and potential annotation."; # status display (with autoflush)
while (my $seq_obj = $seqio_obj->next_seq) { # Bio::Seq object, query might be multi-fasta/embl/genbank
    print '.'; # status display
    my $query_acc;
    if (ref($seqio_obj) =~ m/\:\:fasta$/i) { # $seqio_obj blessed, thus ref returns package name (here= Bio::SeqIO::fasta) from bioperl guessing
	$query_acc = $seq_obj->display_id; # method '->accession_number' doesn't work with fasta files
	$query_acc =~ s/gi\|\d+\|(emb|gb|dbj|ref)\|(.+)\|/$2/; # get acc. nr. from NCBI fasta ID lines
        $query_acc =~ s/\.\d$//; # rm version nr., to fit it to blast acc. nr. above
    } else { # Bio::Seq::RichSeq file (embl, genbank) inherited from Bio::SeqIO::embl/genbank
	$query_acc = $seq_obj->accession_number;
	@{$features{$query_acc}} = $seq_obj->get_all_SeqFeatures; # store all Bio::SeqFeatureI objects in anonymous array
    }
    if ($query_acc ne $blast_q_acc[$multi]) { # check if acc.s fit; SAME order of acc.s in blast-report and seq-obj
	die "\n\n###Fatal error:\nBlast query accession \'$blast_q_acc[$multi]\' doesn't fit to query sequence accession \'$query_acc\'!\n\n";
    }
    $q_seqobj{$query_acc} = $seq_obj;
    $multi++;
}


### Get the contiguous unaligned regions of the query seq(s)
my %uncovered_regions; # stores coords of uncovered regions (RODs)
print "\nLooking for RODs."; # status display (with autoflush)
foreach my $acc (keys %q_seqobj) {
    print "."; # status display
    my $start;
    my $seq = $q_seqobj{$acc}->seq; # Bio::SeqIO object method 'seq'
    my $previous_state = 1; # begin at base #1 of query seq
    foreach my $i (1 .. length($seq)) { # loop through the whole query seq and look for undef positions
       my $current_state = $q_covered{$acc}{$i};
	if (!defined $current_state and defined $previous_state) {
#	    warn "We have just entered an unaligned region $acc: $i\n";
	    $start = $i;
	} elsif (defined $current_state and !defined $previous_state) {
#	    warn "We have just exited an unaligned region $acc: ", $i-1, "\n";
	    my $end = $i - 1; # position before current_state was undef
	    my $length = $end - $start + 1;
	    push @{$uncovered_regions{$length}{$acc}}, $start; # anonymous array in hash-in-hash data structure
	    undef $start;
	} elsif ($i == length($seq) and !defined $previous_state) {
#	    warn "We have just exited an unaligned region at end of contig $acc: $i\n";
	    my $end = $i;
	    my $length = $end - $start + 1;
	    push @{$uncovered_regions{$length}{$acc}}, $start;
	    undef $start;
	}
	$previous_state = $current_state;
    }
}


### Create results directory, where output files are written to
my $out_dir = './results/';
if (-e $out_dir) {
    print "\n\n###Directory \'$out_dir\' already exists! Replace the directory and all its contents? [y|n] ";
    my $user_ask = <STDIN>;
    if ($user_ask =~ /y/i) {
	unlink glob $out_dir."*"; # remove all files in results directory
	rmdir $out_dir; # remove the empty directory
    } else {
	die "Script abborted!\n";
    }
}
mkdir $out_dir or die "Can't create directory \"$out_dir\": $!\n";


### List the longest uncovered regions and print results in output files
my $rod_result = 'rod_summary.txt'; # overview of parsed ROD results and possible ROD annotations
open(ROD, ">$out_dir"."$rod_result") or die "Failed to create $rod_result file: $!\n";
my $seq_out = 'rod_seq.fasta'; # multi-fasta nucleotide file of ROD seqs
open(SEQ, ">$out_dir"."$seq_out") or die "Failed to create $seq_out file: $!\n";
my $gff_out = 'rod.gff'; # GFF3 file of ROD regions for use in artemis/DNAPlotter etc.
open(GFF, ">$out_dir"."$gff_out") or die "Failed to create $gff_out file: $!\n";
my $brig_out = 'rod_BRIG.txt'; # tab-separated ROD file to load into BRIG
open (BRIG, ">$out_dir"."$brig_out") or die "Failed to create $brig_out file: $!\n";
print ROD "$q_file regions with no blast hits (RODs) in $blast_report\n";
print ROD "Rank\t";
if ($multi >=2) { # for multi-seq query files include acc. nr. column
    print ROD "Accession\t";
}
print ROD "Length\tROD position";
my $cds_out; # only for RichSeq query files (embl/genbank)
if (ref($seqio_obj) !~ m/\:\:fasta$/i) { # query NOT a fasta file
    print ROD "\tFeature primary_tag\tFeature locus_tag\tgene\tFeature product\tFeature position";
    $cds_out = 'rod_aa_fasta.txt'; # CDS amino acid sequence output file
    open(CDS, ">$out_dir"."$cds_out") or die "Failed to create $cds_out file: $!\n";
}
print ROD "\n";
print GFF "##gff-version 3\n"; # header, shows it's a gff version 3 file
print GFF "#$q_file regions with no blast hits (RODs) in $blast_report\n"; # comment line for description
print BRIG "#Start\tStop\tLabel\n"; # column headers (commented)

print "\nPreparing results."; # status display (with autoflush)
my $rank = 0; # number of ROD regions
foreach my $length (sort {$b<=>$a} keys %uncovered_regions) { # sort RODs from large length to small length
    if ($length >= $minimum_size) {
	print "."; # status display
	foreach my $acc (sort keys %{$uncovered_regions{$length}}) { # de-reference hash-in-hash, sort by acc. nr.
	    foreach my $start (sort @{$uncovered_regions{$length}{$acc}}) { # de-reference anonymous array
		my $end = $start + $length - 1;
		my $pos = "$start..$end";
		$rank++;
		print GFF "$acc\tBLASTN\tsequence_difference\t$start\t$end\t.\t+\t.\tName=ROD$rank;Target=$blast_db;color=2\n";
		print BRIG "$start\t$end\tROD$rank\n";
		print ROD "$rank\t";
		if ($multi >=2) { # for multi-seq query files include acc. nr.
		    print ROD "$acc\t";
		}
		print ROD "$length\t$pos";
 		print SEQ ">ROD$rank $length $pos\n";
		print SEQ $q_seqobj{$acc}->subseq($start,$end), "\n\n"; # subseq method of Bio::Seq object

                ### fill ROD feature columns in 'rod_result_summary.txt', print CDS aa seqs to 'rod_aa_fasta.txt' for RichSeq files,
                ### and print each ROD as a separate sequence file (with Bio::Sequtils->'trunc_with_features')
		if (ref($seqio_obj) !~ m/\:\:fasta$/i) { # RichSeq query file
		    if ($rank >1) { # blank line in front of the next ROD aa seq block, if it's not the first block
			print CDS "\n";
		    }
		    my $truncseq = Bio::SeqUtils->trunc_with_features($q_seqobj{$acc}, $start, $end); # truncate current Bio::Seq object to coordinates
		    my $seqout; # variable to hold Bio::SeqIO object
		    if (ref($seqio_obj) =~ /genbank/) { # print the truncated ROD sequences in the same '-format' as $q_file; Bio::SeqIO::genbank object
			$seqout = Bio::SeqIO->new(-file => ">$out_dir"."ROD$rank".".gbk", -format => 'genbank');
		    } elsif (ref($seqio_obj) =~ /embl/) { # Bio::SeqIO::embl object
			$seqout = Bio::SeqIO->new(-file => ">$out_dir"."ROD$rank".".embl", -format => 'embl');
		    }
		    $seqout->write_seq($truncseq); # write truncated sequence object
		    my $loop = 0; # indicates that more than one feature in the ROD range
		    print CDS "~~ROD$rank $length $pos\n";
		    foreach my $feature (@{$features{$acc}}) { # de-reference anonymous array in hash
			if ($feature->location->start >= $start && $feature->location->end <= $end) { # features that are fully within ROD region
			    $loop = print_CDSs($feature, $loop); # subroutine 'print_CDSs'
			} elsif (($feature->location->start <= $start && ($feature->location->end > $start && $feature->location->end <= $end)) || (($feature->location->start >= $start && $feature->location->start < $end) && $feature->location->end > $end)) { # features that overlap the ROD region
			    $loop = print_CDSs($feature, $loop, 1); # 1 defines $overlap in subroutine
			}
		    }
		    if ($loop == 0) { # print "\n" if no features exist in the ROD
			print ROD "\n";
		    }
		} else { # fasta files don't contain feature information
		    print ROD "\n";
		}
	    }
	}
    }
}
close ROD;
close SEQ;
close GFF;
close BRIG;
if (ref($seqio_obj) !~ m/\:\:fasta$/i) {
    close CDS;
}
$| = 0; # turn off autoflush


### Print which files have been created!
print "\n\nThe following files were created in the \"./results\" directory:\n";
print "- $rod_result: Tab-separated summary of the found regions of difference (RODs)\n";
print "- $gff_out: GFF3 file with ROD coordinates to use in Artemis/DNAPlotter ...\n";
print "- $brig_out: Tab-separated file for BRIG (BLAST Ring Image Generator)\n";
print "- $seq_out: Nucleotide sequences of ROD regions in multi-fasta format\n";
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) {
    print "- $cds_out: Amino acid sequences of CDSs in ROD regions\n";
    print "- ROD#.format: Each ROD as a separate sequence file in the query file format\n";
}
print "\n";

exit;


###############
# Subroutines #
###############

### Subroutine to print the position of features according to leading and lagging strand
sub position {
    my $feature = shift;
    if ($feature->strand == 1) { # leading strand
	print ROD "\t", $feature->location->start, "..", $feature->location->end;
    } elsif ($feature->strand == -1) { # lagging strand
	print ROD "\t", $feature->location->end, "..", $feature->location->start;
    }
    return 1;
}


### Subroutine that prints feature information from RichSeq files
sub print_CDSs {
    my $feature = shift;
    my $loop = shift;
    my $overlap = shift; # defined if feature overlaps ROD edge
    my $primary_tag = $feature->primary_tag;
    if (($primary_tag eq 'CDS') || ($primary_tag eq 'misc_RNA') || ($primary_tag eq 'ncRNA') || ($primary_tag eq 'rRNA') || ($primary_tag eq 'tRNA') || ($primary_tag eq 'tmRNA') || ($primary_tag eq 'misc_feature') || ($primary_tag eq 'repeat_region')) {
	skip_columns($loop); # sub to skip to correct column in print out
	print ROD "\t$primary_tag";
	if ($feature->has_tag('locus_tag')) { # TRUE if tag 'locus_tag' exists for this feature
	    print ROD "\t", $feature->get_tag_values('locus_tag');
	    print CDS ">", $feature->get_tag_values('locus_tag') if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo')); # include only CDSs with "/translation" in file 'rod_aa_fasta.txt' (exclude pseudo-genes, RNAs ...)
	} else {
	    print ROD "\t";
	    print CDS ">" if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	}
	if ($feature->has_tag('gene')) {
	    print ROD "\t", $feature->get_tag_values('gene');
	    print CDS " ", $feature->get_tag_values('gene') if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	} else {
	    print ROD "\t";
	    print CDS " " if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	}
	if ($feature->has_tag('product')) {
	    my ($product) = $feature->get_tag_values('product');
	    print ROD "\t", $product;
	    print CDS " ", replace_white($product) if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo')); # subroutine 'replace_white'
	} elsif ($feature->has_tag('rpt_family')) { # don't include feature 'rpt_family' in the aa output
	    print ROD "\t", $feature->get_tag_values('rpt_family');
	} elsif ($feature->has_tag('note')) { # only include feature 'note' if 'product' doesn't exist
	    my ($note) = $feature->get_tag_values('note');
	    print ROD "\t", $note;
	    print CDS " ", replace_white($note) if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	}
	position($feature); # subroutine 'position' to print coordinates of features
	if (defined $overlap) { # include '(overlap)' if feature overlaps edge of ROD
	    print ROD " (overlap)\n";
	    print CDS " (overlap)\n" if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	} else {
	    print ROD "\n";
	    print CDS "\n" if ($primary_tag eq 'CDS' && !$feature->has_tag('pseudo'));
	}
	if ($primary_tag eq 'CDS' && $feature->has_tag('translation')) { # print translations of CDSs
	    print CDS $feature->get_tag_values('translation'), "\n";
	}
	$loop++; # Indicates that more than one feature in the ROD range
    } elsif ($primary_tag eq 'mobile_element') { # include 'mobile_element' features only in ROD
	skip_columns($loop);
	print ROD "\t$primary_tag";
	if ($feature->has_tag('mobile_element_type')) {
	    print ROD "\t\t\t", $feature->get_tag_values('mobile_element_type');
	} else {
	    print ROD "\t\t\t";
	}
	position($feature); # sub to print positions of features
	if (defined $overlap) {
	    print ROD " (overlap)\n";
	} else {
	    print ROD "\n";
	}
	$loop++;
    }
    return $loop;
}


### Subroutine to replace whitespaces in a string
sub replace_white {
    my $string = shift;
    $string =~ s/\s/_/g;
    return $string;
}


### Subroutine to skip to the correct column for features in tab-separated result file 'rod_result_summary.txt'
sub skip_columns {
    my $loop = shift; # $multi is in the same scope ... (actually also $loop)
    if ($loop >= 1) { # for the features after the first ROD line skip to the 'feature primary_tag' column in the tab-separated summary output
	print ROD "\t\t";
	if ($multi >=2) { # multi-seq query files include acc. nr. column
	    print ROD "\t";
	}
    }
    return 1;
}
