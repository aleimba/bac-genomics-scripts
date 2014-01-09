#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SearchIO; # bioperl module to handle blast reports
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqFeatureI; # bioperl module to handle features in a sequence (only possible for rich sequence formats, like embl/genbank, not fasta)

my $usage = "\n".
   "\t#####################################################################\n".
   "\t# $0 blastn_report query-seqfile min_ROD_size      #\n". #$0 = program name
   "\t#                                                                   #\n".
   "\t# Searches for regions of difference (RODs) in a query sequence     #\n".
   "\t# according to a blastn run. Give the BLAST query sequence file in  #\n".
   "\t# either fileformat (multi-fasta, -embl or -genbank) to the script. #\n".
   "\t# embl/genbank queries also return CDSs in RODs.                    #\n".
   "\t# minimum_ROD_size = minimum size of detected RODs (in bp)          #\n".
   "\t# The perl script uses bioperl (www.bioperl.org).                   #\n".
   "\t#                                                                   #\n".
   "\t# version 0.3, update: 23.11.2011  Andreas Leimbach/David Studholme #\n".
   "\t# 07.11.2011                                  aleimba[at]gmx[dot]de #\n".
   "\t#####################################################################\n\n";

### ToDo:
# - 'blasttable' (more efficient runtime) and 'blastxml' (accepts standard blast report format changes) blastreport formats don't work. Might depend on bioperl version. --> possible to check blast report format?


### Print usage if -h|--h|--help is given as argument or arguments are not given
my ($blast_report, $q_file, $minimum_size) = @ARGV; # minimum size of uncovered regions that are reported
if (!defined($blast_report) || !defined($q_file)) {
    die $usage;
} elsif ($blast_report =~ m/-h/) {
    die $usage;
}



### Parse the blast report/output file
my %q_covered; # hash which stores for each position (=hash-key) of the query sequence if it falls within a hsp hit
my @blast_q_acc; # test if identical query accessions in blast parse and seqio (below)
my $blast_db;
my $parser = new Bio::SearchIO(-file => "<$blast_report", -format => 'blast'); # Bio::SearchIO object
$|++; # turn on autoflush, forces STDOUT flush right away (STDOUT normally line-/block-buffered)
print "\nParsing blast report file."; # status display (with autoflush)
while(my $result = $parser->next_result) { # several query sequences possible (result = entire analysis for a single query seq!) -> e.g. usable with multi-fasta query file
    my $query_acc = $result->query_accession;
    $query_acc =~ s/\.$//; # somehow with contig ID-lines '->query_accession' appends a '.'
    $query_acc =~ s/\.\d$//; # rm version number from query accession (to fit it to the acc.-nr. from Bio::SeqIO below)
    push(@blast_q_acc, $query_acc);
    $blast_db = $result->database_name; # For use in the GFF3 print out
    chop $blast_db; # always an empty character after the name?
    while(my $hit = $result->next_hit) { # several subject sequences in the database might have hits
	print '.'; # status display
	while(my $hsp = $hit->next_hsp) { # each hit might have one or more hsps (the alignments shown in a blast report)
	    my $hit_start = $hsp->start('hit'); # reference/subject hsp start coordinate
	    my $hit_end = $hsp->end('hit');
	    my ($query_start, $query_end) = $hsp->range('query'); #equal to 'start/end' above
	    foreach my $i ($query_start .. $query_end) {
		$q_covered{$query_acc}{$i}++; # changes for each hsp hit the value from 'undef' to '1'
	    }
	}
    }
}


### Read the query sequence file into RAM, for Bio::Seq::RichSeq query files (e.g. embl/genbank) also read features
my %q_seqobj; # hash that stores all seqobj as values in reference to acc.-nr.s as keys 
my %features; # hash which stores all features for RichSeq files
my $multi = 0; # test if it is a multi sequence query file (multi-fasta/embl/genbank ...) -> $multi >= 2
my $i = 0; # needed to test for identical acc.-nr.s
my $seqio_obj = Bio::SeqIO->new(-file => "<$q_file"); # didn't use '-format' to leave it to bioperl guessing
print "\nParsing query sequence file."; # status display (with autoflush)
while (my $seq_obj = $seqio_obj->next_seq) { # Bio::Seq object, query might be multi-fasta/embl/genbank
    print '.'; # status display
    my $query_acc; # fits to $query_acc above in blast-result parse
    if (ref($seqio_obj) =~ m/\:\:fasta$/i) { # $seqio_obj blessed, thus ref returns package name (here: Bio::SeqIO::fasta)
	$query_acc = $seq_obj->display_id; # Bio::Seq method, '->accession_number' doesn't work with fasta files
	$query_acc =~ s/gi\|\d+\|(emb|gb|dbj|ref)\|(.+)\|/$2/; # get the accession number from the fasta ID line
        $query_acc =~ s/\.\d$//; # rm version nr., to fit it to the blast acc. number above
    } else { # Bio::Seq::RichSeq file (embl, genbank ...) inherited from Bio::SeqIO::embl/genbank
	$query_acc = $seq_obj->accession_number;
	@{$features{$query_acc}} = $seq_obj->get_all_SeqFeatures; # store all Bio::SeqFeatureI objects in anonymous array
    }
    if ($query_acc ne $blast_q_acc[$i]) { # Script stops if the accessions don't fit together!
	die "\nFatal error:\nBlast query accession \'$blast_q_acc[$i]\' doesn't fit to query sequence accession \'$query_acc\'!\n\n";
    }
    $q_seqobj{$query_acc} = $seq_obj; # hash with the sequence objects, multi-fasta/embl/genbank query possible
    $multi++;
    $i++;
}


### Get the contiguous unaligned regions of query seqs
my %uncovered_regions; # hash to store the uncovered regions
print "\nLooking for RODs."; # status display (with autoflush)
foreach my $acc (keys %q_seqobj) {
    print "."; # status display
    my $start;
    my $seq = $q_seqobj{$acc}->seq; # Bio::SeqIO object method 'seq'
    my $previous_state = 1;
    foreach my $i (1 .. length($seq)) {
       my $current_state = $q_covered{$acc}{$i};
	if (!defined $current_state and defined $previous_state) {
#	    warn "We have just entered an unaligned region $acc: $i\n";
	    $start = $i;
	} elsif (defined $current_state and !defined $previous_state) {
#	    warn "We have just exited an unaligned region $acc: ", $i-1, "\n";
	    my $end = $i - 1; # because at position before current_state was undef
	    my $length = $end - $start + 1;
	    push @{$uncovered_regions{$length}{$acc}}, $start; # anonymous array in hash data structure
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


### List the longest uncovered regions and print results in output files
my $rod_result = 'rod_result_summary.txt'; # Overview of the parsed ROD results
open(ROD, ">$rod_result") or die "Failed to create $rod_result file: $!\n";
my $seq_out = 'rod_seq.fasta'; # Multi-fasta file of the ROD sequences
open(SEQ, ">$seq_out") or die "Failed to create $seq_out file: $!\n";
my $gff_out = 'rod_result.gff'; # GFF3 file of ROD regions for use in artemis/DNAPlotter etc.
open(GFF, ">$gff_out") or die "Failed to create $gff_out file: $!\n";
my $brig_out = 'rod_BRIG.txt'; # tab-seperated file to load into BRIG
open (BRIG, ">$brig_out") or die "Failed to create $brig_out file: $!\n";
print ROD "$q_file regions with no blast hits (RODs) in $blast_report\n";
print ROD "Rank\tLength\tROD position";
my $cds_out;
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # if not a fasta file as query (see above)
    print ROD "\tROD CDSs\tCDSs position";
    $cds_out = 'rod_aa_fasta.txt'; # CDS amino acid sequence output file only makes sense if query file not a fasta file
    open(CDS, ">$cds_out") or die "Failed to create $cds_out file: $!\n";
}
print ROD "\n";
print GFF "##gff-version 3\n"; # header, shows it's a gff version 3 file
print GFF "#$q_file regions with no blast hits (RODs) in $blast_report\n"; # comment line for description
print BRIG "#Start\tStop\tLabel\n"; # comment line as column headers

print "\nPreparing results."; # status display (with autoflush)
my $rank = 0; # number of ROD regions
#my $minimum_size = 100; # minimum size of uncovered regions that are reported
foreach my $length (sort {$b<=>$a} keys %uncovered_regions) { # sort from large length to small length
    if ($length >= $minimum_size) {
	print "."; # status display
	foreach my $acc (sort keys %{$uncovered_regions{$length}}) {
	    foreach my $start (sort @{$uncovered_regions{$length}{$acc}}) {
		my $end = $start + $length - 1;
		my $pos;
		$rank++;
		if ($multi == 1) { # acc.-nr. only included in ROD CDS column of result_summary if query multi-seq file
		    $pos = "$start..$end";
		} elsif ($multi >= 2) {
		    $pos = "$acc: $start..$end";
		}
		print GFF "$acc\tBLASTN\tsequence_difference\t$start\t$end\t.\t+\t.\tName=ROD$rank;Target=$blast_db;color=2\n"; ### Include ...$rank;Target= $hit_acc $hit_start $hit_stop;color ...?
		print BRIG "$start\t$end\tROD$rank\n";
		print ROD "$rank\t$length\t$pos";
 		print SEQ ">ROD$rank\_$length\_$pos\n";
		print SEQ $q_seqobj{$acc}->subseq($start,$end), "\n\n"; # subseq method of Bio::Seq object

                ### Information about the CDS, only makes sense if query file is NOT a fasta file, but a RichSeq file
		if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) {
		    if ($rank >1) { # another \n in front of the next ROD aa seq block, if it's not the first block
			print CDS "\n";
		    }
		    my $loop = 0;
		    my $features_exist = 0; # Are features in the ROD? If not print "\n"!
		    print CDS "--ROD$rank\_$length\_$pos\n";
		    foreach my $feature (@{$features{$acc}}) {
			if ($feature->location->start >= $start && $feature->location->end <= $end) { # features that are fully within ROD region
			    $features_exist = 1;
			    $loop = print_CDSs($feature, $loop); # subroutine
			} elsif (($feature->location->start <= $start && ($feature->location->end > $start && $feature->location->end <= $end)) || (($feature->location->start >= $start && $feature->location->start < $end) && $feature->location->end > $end)) { # features that overlap the ROD region, either only with their end or start coordinates
			    $features_exist = 1;
			    $loop = print_CDSs($feature, $loop, 1);
			}
		    }
		    if ($features_exist == 0) { # Print "\n" if no features exist in the ROD
			print ROD "\n";
		    }
                    ### Print out each ROD in the query sequence fileformat ==> but didn't work well!
#		    my $outfile = "ROD$rank.txt";
#		    my $out = ref($seqio_obj)->new(-file => ">$outfile"); # no format would be needed in the argument, because already decided by the module!
#		    $out->write_seq($q_seqobj{$acc}->trunc($start,$end));
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
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) {
    close CDS;
}
$| = 0; # turn of autoflush

### Print which files have been created!
print "\n\nThe following files were created:\n";
print "- $rod_result: Summary of the found regions of difference (RODs)\n";
print "- $gff_out: GFF3 file with ROD coordinates to use in artemis/DNAPlotter ...\n";
print "- $brig_out: Tab-seperated file for BRIG (blast ring image generator)\n";
print "- $seq_out: Fasta sequences of the ROD regions\n";
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) {
    print "- $cds_out: Fasta amino acid sequences of CDSs in ROD regions\n";
}
print "\n";


exit;


#############
#Subroutines#
#############

### Subroutine that prints feature information from RichSeq files
sub print_CDSs {
    my $feature = shift;
    my $loop = shift;
    my $overlap = shift;
    my $primary_tag = $feature->primary_tag; # variable just to make the code shorter
    if (($primary_tag eq 'CDS') || ($primary_tag eq 'tRNA') || ($primary_tag eq 'rRNA') || ($primary_tag eq 'tmRNA') || ($primary_tag eq 'ncRNA') || ($primary_tag eq 'repeat_region') || ($primary_tag eq 'misc_RNA') || ($primary_tag eq 'misc_binding')) {
	if ($loop >= 1) { # for the features after the first line skip to the 'ROD CDSs' column
	    print ROD "\t\t";
	}
	if ($feature->has_tag('product')) { # TRUE if tag 'product' exists for this feature, otherwise error is thrown
	    print ROD "\t", $feature->get_tag_values('product');
	    print CDS ">", $feature->get_tag_values('product');
	} elsif ($feature->has_tag('note')) {
	    print ROD "\t", $feature->get_tag_values('note');
	    print CDS ">", $feature->get_tag_values('note');
	} elsif ($feature->has_tag('rpt_family')) {
	    print ROD "\t", $feature->get_tag_values('rpt_family');
	    print CDS ">", $feature->get_tag_values('rpt_family');
	}
	print ROD "\t", $feature->location->start, "..", $feature->location->end;
	if (defined $overlap) {
	    print ROD " (overlap)\n";
	    print CDS " (overlap)\n";
	} else {
	    print ROD "\n";
	    print CDS "\n";
	}
	if (($primary_tag eq 'CDS') && ($feature->has_tag('translation'))) {
	    print CDS $feature->get_tag_values('translation'), "\n"; # translation only defined for CDSs not RNAs ...
	} else {
	    print CDS "No translation available\n";
	}
	$loop++; # Indicates that more than one feature in the ROD range
	return $loop;
    } elsif ($primary_tag eq 'mobile_element') {
	print ROD "\t", $feature->get_tag_values('mobile_element_type') if $feature->has_tag('mobile_element_type');
	print CDS ">", $feature->get_tag_values('mobile_element_type') if $feature->has_tag('mobile_element_type');
    }
    return $loop;
}
