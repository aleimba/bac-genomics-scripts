#!/usr/bin/perl

use strict;
use warnings;
# use CGI; # relic from the html print out on the bottom
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SearchIO; # bioperl module to handle blast reports
use Bio::Seq; # bioperl module to handle sequences with features
use Bio::SeqFeatureI; # bioperl module to handle features in a sequence (only possible for rich sequence formats, like embl/genbank, not fasta)

my $usage = "\n".
   "\t###################################################################\n".
   "\t# $0 blastn.out query-seq_file                     #\n". #$0 = program name
   "\t#                                               #\n".
   "\t# blastformat? ('blast', 'blasttable', 'blastxml')          #\n".
   "\t# Searches for uncovered regions in a blastn output file in      #\n".
   "\t# correspondence to a given query sequence file (either format: fasta, embl, genbank) #\n".
   "\t# (which has to correspond     #\n".
   "\t# to the query.fasta file in the blastn run) (use seq_format-converte.pl!)   #\n".
#   "\t# If embl/genbank files are given as a query also returns CDSs in ROD regions (overlap marks!)!           #\n".
   "\t#  ROD = region of difference                                     #\n".
   "\t#  for the blast used an evalue cutoff of '-e 0.1'                #\n".
   "\t#                                               #\n".
   "\t# The perl script uses bioperl!                                  #\n".
   "\t#                                                                 #\n".
   "\t# blastn.out = report/output file from a blast analysis          #\n".
   "\t# query.fasta = query fasta file                                 #\n".
   "\t#                                                                 #\n".
   "\t# version 0.1                   Andreas Leimbach/David Studholme  #\n".
   "\t# 07.11.2011                               aleimba[at]gmx[dot]de  #\n".
   "\t###################################################################\n\n";

#print usage if -h|--h|--help is given as argument or arguments are not given
my ($blast_report, $q_file) = @ARGV;
if (!defined($blast_report) || !defined($q_file)) {
    die $usage;
} elsif ($blast_report =~ m/-h/) {
    die $usage;
}

#############################################################################################
# Old code from David Studholme is labeled with '#####'!
# Studholme included in the script also the creation of a html file with perl package CGI.
# In this html file the unique regions a directly linked to an installation
# to GBrowser, probably to directly have a possibility to look at each region in GBrowse!
#############################################################################################

my $minimum_size = 100; # minimum size of uncovered regions, that should be shown######give in arguments?

##### my $sequence_file = $q_file;
##### my $gbrowse_name = 'Erwinia_toletana';
##### $gbrowse_name = 'Xsp1132' if $q_file =~ m/Xsp1132/;
##### $gbrowse_name = 'Xsp1131' if $q_file =~ m/Xsp1131/;
##### $gbrowse_name = 'Xsp4393' if $q_file =~ m/Xsp4393/;
##### warn "Using Gbrowse $gbrowse_name\n";


my %q_covered; # hash which stores for each position (=key) of the query genome if it falls within a hsp hit, by setting the corresponding hash value to 1!
# When finished parsing the blast output any positions in the hash that remains undefined must be non matching regions.

#########
# Filter for good hits and/or good hsps (e-value, length, percent id ...) -> make an e-value cutoff with the blast program call (write in program description!)? What about repeats, alignment might be to only one repeat in the query (e.g. use method amgibuous_aln? See HOWTO:SearchIO)? But maybe not important if only unique regions in the query are needed?
#########

### Parse the blast report/output file
my $parser = new Bio::SearchIO(-file => "<$blast_report", -format => 'blast'); # Bio::SearchIO object ### better for perfomance to use blast -m8/9 (-format => 'blasttable') (or at least to inhibit failures via possible formatchange use 'blastxml') --> possible to check blast report format?
while(my $result = $parser->next_result) { # several query sequences possible (result = entire analysis for a single query seq!) -> e.g. usable with multi-fasta query file!
    my $query_acc = $result->query_accession;
    $query_acc =~ s/\.\d$//; # rm version number from query accession (in order to fit it to the accession number recieved from Bio::SeqIO below!)
    warn "query_acc blast: $query_acc\n";
    while(my $hit = $result->next_hit) { # several subject sequences in the database might have hits!
#	my $hit_acc = $hit->accession; #### needed? $hit_acc doesn't come up again!

#	my $ambiguous_aln = $hit->ambiguous_aln; ### possibility to control if a hit has overlapping hsp hits (see HOWTO:SearchIO)!
#	print "ambiguous align = $ambiguous_aln\n";

	while(my $hsp = $hit->next_hsp) { # each hit might have one or more hsps (the alignemnts shown in a blast report)!
	    my $hit_start = $hsp->start('hit'); # reference/subject hsp start coordinate; also here range possible (see below)
	    my $hit_end = $hsp->end('hit');
	    my ($query_start, $query_end) = $hsp->range('query'); ### should be eq to $hsp->length('query'); -> Control?!
#	    my $query_start = $hsp->start('query');
#	    my $query_end = $hsp->end('query');
	    warn "hitlength: ", $hit_end-$hit_start, "\tquerylength: ", $query_end-$query_start, "\n";
	    my $query_string = $hsp->query_string; # not used afterwards!
	    my $hit_string = $hsp->hit_string; # not used afterwards!
	    my $frac_identical = $hsp->frac_identical; # not used afterwards!

#####	    my ($tag_r, $start_r, $start_q, $length) = ($1, $2, $3, $4); # must be old code from David, makes no sense here ...

#if ($hsp->evalue > 0.0001) {} # filter for low e-values (length ...) here? Or do this directly in the blast call?

	    foreach my $i ($query_start .. $query_end) {
		$q_covered{$query_acc}{$i}++; # changes for each query sequence position the value from 'undefined' to '1'!
		# warn "\$q_covered{$query_acc}{$i}++\t";
	    }
	}
    }
}


### Read the query sequence file into RAM, for RichSeq query files (e.g. embl/genank) also read features into RAM
my %q_seqobj; # hash that stores all seqobj, so '->seq',  '->subseq/trunc' (for ROD seqs) can be used later on
my %features; # hash which stores all features for RichSeq files
my $seqio_obj = Bio::SeqIO->new(-file => "<$q_file"); # Bio::SeqIO object; didn't use '-format' to leave it to bioperl guessing
my $multi; # test if it is a multi sequence query file (e.g. multi-fasta/embl/genbank ...)
while (my $seq_obj = $seqio_obj->next_seq) { # a Bio::Seq object, query might be multi-fasta/embl/genbank
    my $query_acc; # should fit to $query_acc above in blast-result parse!
    warn "ref-seqio_obj: ", ref($seqio_obj), "\n";
    if (ref($seqio_obj) =~ m/\:\:fasta$/i) { # ref returns value if $seq_obj is a reference, and singe object is blessed into package it returns package name (here: Bio::SeqIO::fasta)
	$query_acc = $seq_obj->display_id; # Bio::Seq methods, '->accession_number' doesn't work with fasta files
	$query_acc =~ s/gi\|\d+\|(emb|gb|dbj|ref)\|(.+)\|/$2/; # get the accession number from the fasta ID line#########what if multi-fasta is not E. coli genome IDs but contig names ...?!
        $query_acc =~ s/\.\d$//; # rm version number from query accession (see above) (in order to fit it to the accession number recieved from Bio::SeqIO below!)
    } else { # RichSeq file (embl, genbank ...)
	$query_acc = $seq_obj->accession_number; # doesn't work with fasta files
	@{$features{$query_acc}} = $seq_obj->get_all_SeqFeatures; # store all features in anonymous array of Bio::SeqFeatureI objects
    }
    $q_seqobj{$query_acc} = $seq_obj; # anonymous array with the sequence objects, multi-fasta/embl/genbank query possible
    $multi++;
}


### Get the contiguous unaligned regions of query seqs
my %uncovered_regions; # hash to store the uncovered regions
foreach my $acc (sort keys %q_seqobj) {
    my $start;
    my $seq = $q_seqobj{$acc}->seq; # Bio::SeqIO object with method 'seq'
    my $previous_state = 1;
    foreach my $i (1 .. length($seq)) { ####there was a '-1' here? why????
       my $current_state = $q_covered{$acc}{$i};
	if (defined $current_state and !defined $previous_state) {
	    warn "We have just exited an unaligned region $acc: ", $i-1, "\n";
	    my $end = $i - 1; # because $end at the position before current_state defined!
	    my $length = $end - $start + 1;
	    push @{$uncovered_regions{$length}{$acc}}, $start; # anonymous array in hash data structure!
	    undef $start;
#####	    $start = undef;
	}
	if (!defined $current_state and defined $previous_state) {
	    warn "We have just entered an unaligned region $acc: $i\n";
	    $start = $i;
	}
	if ($i == length($seq) and !defined $previous_state) {
	    warn "We have just exited an unaligned region at end of contig $acc: $i\n";
#	    my $end = $i - 1; #### needed here?
	    my $end = $i;
	    my $length = $end - $start + 1;
	    push @{$uncovered_regions{$length}{$acc}}, $start;
	    undef $start;
#####	    $start = undef;
	}
	$previous_state = $current_state;
    }
}


#List the longest uncovered regions and print in text files
my $rod_result = 'rod_result.txt'; # Overview of the ROD results
open(ROD, ">$rod_result") or die "Failed to create $rod_result file: $!\n";
my $seq_out = 'rod_seq.fasta';
open(SEQ, ">$seq_out") or die "Failed to create $seq_out file: $!\n";
my $gff_out = 'rod_result.gff';
open(GFF, ">$gff_out") or die "Failed to create $gff_out file: $!\n";
print ROD "$q_file regions with no blast hits (RODs) in $blast_report\n";
print ROD "Rank";
print ROD "\tLength";
print ROD "\tROD Position";
print GFF "##gff-version 3\n"; # shows it is a gff version 3 file
print GFF "#$q_file regions with no blast hits (RODs) in $blast_report\n"; # comment line for description
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # if not a fasta file as query (see above)
    print ROD "\tROD CDSs";
    print ROD "\tCDSs start:stop";
    my $cds_out = 'aa_seq_fasta.txt'; # CDS amino acid sequence output file only makes sense if query file not a fasta file
    open(CDS, ">$cds_out") or die "Failed to create $cds_out file: $!\n";
}
print ROD "\n";
my $rank = 0;
foreach my $length (sort {$b<=>$a} keys %uncovered_regions) { # sort from large length to small length
    if ($length >= $minimum_size) {
	foreach my $acc (sort keys %{$uncovered_regions{$length}}) {
	    foreach my $start (@{$uncovered_regions{$length}{$acc}}) { #### Sort also the start coordinates?!
		my $end = $start + $length - 1;
		my $pos;
		if ($multi <= 1) { # accession-nr. only included if the query is a multi-sequence file
		    $pos = "$start..$end";
		} elsif ($multi >= 2) {
		    $pos = "$acc: $start..$end";
		}
		$rank++;
                ### Print gff3 for use in artemis ...
		print GFF "$acc\tBLASTN\tsequence_difference\t$start\t$end\t.\t+\t.\tName=ROD$rank;Target=;color=2\n";

		print ROD "$rank";
		print ROD "\t$length";
		print ROD "\t$pos";
 		print SEQ ">ROD$rank\_$length\_$pos\n";
		print SEQ $q_seqobj{$acc}->subseq($start,$end), "\n\n"; # Bio::SeqIO object, multi-fasta file of 
                ### Information about the CDS, only makes sense if query file is NOT a fasta file, but a RichSeq file
		if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # see above
		    if ($rank >1) { # another \n in front of the next ROD CDS block, if it's not the first block
			print CDS "\n";
		    }

                    ### Print out each ROD in the query sequence fileformat ==> but didn't work well!
#		    my $outfile = "ROD$rank.txt";
#		    my $out = ref($seqio_obj)->new(-file => ">$outfile"); # no format would be needed in the argument, because already decided by the module!
#		    $out->write_seq($q_seqobj{$acc}->trunc($start,$end));

		    my $loop = 0;
		    my $counter;
		    print CDS "--ROD$rank\_$length\_$pos\n";
		    foreach my $feature (@{$features{$acc}}) {
			if ($feature->location->start >= $start && $feature->location->end <= $end) {
			    if ($feature->primary_tag eq 'CDS') {
				if ($loop >= 1) { # for the features after the first line skip to the 'ROD CDSs' column
				    print ROD "\t\t";
				}
				print ROD "\t", $feature->get_tag_values('product');
				print CDS ">", $feature->get_tag_values('product'), "\n";
				print ROD "\t", $feature->location->start, "..", $feature->location->end, "\n";
				print CDS $feature->get_tag_values('translation'), "\n";
				print "Counter inside: ", $counter++, " ", $feature->get_tag_values('product'), ", loop: $loop\n";
				$loop = 1;
			    }
			} elsif (($feature->location->start <= $start && ($feature->location->end > $start && $feature->location->end <= $end)) || (($feature->location->start >= $start && $feature->location->start < $end) && $feature->location->end > $end)) {
			    if ($feature->primary_tag eq 'CDS') {
				if ($loop >= 1) {
				    print ROD "\t\t";
				}
				print ROD "\t", $feature->get_tag_values('product'), " (overlap)";
				print CDS ">", $feature->get_tag_values('product'), "\n";
				print ROD "\t", $feature->location->start, ":", $feature->location->end, "\n";
				print CDS $feature->get_tag_values('translation'), "\n";
				print "Counter end-overlap: ", $counter++, " ", $feature->get_tag_values('product'), ", loop: $loop\n";
				$loop = 1;
			    }
			}
		    }
		}
	    }
	}
    }
}
close ROD;
close SEQ;
close GFF;
if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # if not a fasta file as query (see above)
    close CDS;
}








##### Old code for printing out the result in a html file!
# List the longest uncovered regions and print in html-file
# my $rod_result = 'rod_result.html';
# open(ROD, ">$rod_result") or die "Failed to create $rod_result file: $!\n";
# my $cgi = new CGI;
# print ROD $cgi->start_html(-title => "$q_file regions with no blast hits (RODs) in $blast_report");

# ##### print $cgi->start_html(-title => "$gbrowse_name: $q_file regions with no blast hits in $blast_report");
# ##### print "<link rel=\"stylesheet\" href=\"http://www.sainsbury-laboratory.ac.uk/sl.css\">\n";

# print ROD "\n<h3>$q_file regions with no blast hits (RODs) in $blast_report</h3>\n";

# ##### print "\n<h3>$gbrowse_name: $q_file regions with no blast hits in $blast_report</h3>\n";

# print ROD "\n<table border=1>\n";
# print ROD "<tr>";
# print ROD "<th>Rank</th>";
# print ROD "<th>Length</th>";
# print ROD "<th>ROD Position</th>";
# if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # see above
#     print ROD "<th>ROD CDSs</th>";
#     print ROD "<th>CDSs start:stop</th>";
#     print ROD "<th>Feature aa sequence</th>"; # Amino acid sequence
# }
# print ROD "</tr>";
# my $i = 0;
# foreach my $length (sort {$b<=>$a} keys %uncovered_regions) { # sort from large length to small length
#     if ($length >= $minimum_size) {
# 	foreach my $acc (sort keys %{$uncovered_regions{$length}}) {
# 	    foreach my $start (@{$uncovered_regions{$length}{$acc}}) {
# 		my $end = $start + $length - 1; ### why length-1 again?
# 		my $pos;
# 		if ($multi <= 1) { # accession# only included if the query is a multi-sequence file
# 		    $pos = "$start..$end";
# 		} elsif ($multi >= 2) {
# 		    $pos = "$acc: $start..$end";
# 		}

# #####		my $url="http://144.173.144.7/cgi-bin/gb2/gbrowse/$gbrowse_name?name=$pos";
# #####		#my $url="http://jic55119.jic.bbsrc.ac.uk/cgi-bin/gbrowse/$gbrowse_name?name=$pos";
# #####		my @CDS = get_CDSs($pos, $gbrowse_name);

# 		$i++;
# 		print ROD "<tr>";
# 		print ROD "<td>$i</td>";
# 		print ROD "<td>$length</td>";
# 		print ROD "<td>$pos</td>";

# # Include a link here to the sequences of the ROD included genes, i.e. print them out in another html file (instead of printing them out also in rod_result.html'). This should enhance visibility.
# #####		print ROD "<td><a href=\"$url\" target=_blank>$pos</a></td>";
# ### Information about the CDS

# 		if (!(ref($seqio_obj) =~ m/\:\:fasta$/i)) { # see above
# 		    my $loop = 0;
# 		    my $counter;
# 		    foreach my $feature (@{$features{$acc}}) {
# 			if ($feature->location->start >= $start && $feature->location->end <= $end) { #####what about features that overlap start and stop coordinates?!
# 			    if ($feature->primary_tag eq 'CDS') {
# 				if ($loop >= 1) {
# 				    print ROD "<tr>", "<td></td><td></td><td></td>";
# 				}
# 				print ROD "<td><b>", $feature->get_tag_values('product'), "</b></td>";
# 				print ROD "<td>", $feature->location->start, ":", $feature->location->end, "</td>";
# 				print ROD "<td>", $feature->get_tag_values('translation'), "</td>";
# 				print ROD "</tr>\n";
# 				print "Counter inside: ", $counter++, " ", $feature->get_tag_values('product'), ", loop: $loop\n";
# 				$loop = 1;
# 			    }
# 			} elsif ($feature->location->start <= $start && ($feature->location->end > $start && $feature->location->end <= $end)) {
# 			    if ($feature->primary_tag eq 'CDS') {
# 				if ($loop >= 1) {
# 				    print ROD "<tr>", "<td></td><td></td><td></td>";
# 				}
# 				print ROD "<td><b>", $feature->get_tag_values('product'), "</b> (overlap)</td>";
# 				print ROD "<td>", $feature->location->start, ":", $feature->location->end, "</td>";
# 				print ROD "<td>", $feature->get_tag_values('translation'), "</td>";
# 				print ROD "</tr>\n";
# 				print "Counter end-overlap: ", $counter++, " ", $feature->get_tag_values('product'), ", loop: $loop\n";
# 				$loop = 1;
# 			    }
# 			} elsif (($feature->location->start >= $start && $feature->location->start < $end) && $feature->location->end > $end) {
# 			    if ($feature->primary_tag eq 'CDS') {
# 				if ($loop >= 1) {
# 				    print ROD "<tr>", "<td></td><td></td><td></td>";
# 				}
# 				print ROD "<td><b>", $feature->get_tag_values('product'), "</b> (overlap)</td>";
# 				print ROD "<td>", $feature->location->start, ":", $feature->location->end, "</td>";
# 				print ROD "<td>", $feature->get_tag_values('translation'), "</td>";
# 				print ROD "</tr>\n";
# 				print "Counter start-overlap: ", $counter++, " ", $feature->get_tag_values('product'), ", loop: $loop\n";
# 				$loop = 1;
# 			    }
# 			}
# 		    }
# 		}
# 	    }
# 	}
#     }
# }
# print ROD "\n</table>\n";
# print ROD $cgi->end_html;
# close ROD;




### Get the features in the longest uncovered regions
# Print out a GFF3 file with ROD region coordinates for use e.g. in Artemis/DNAPlotter/BRIG ...!


exit;


#############
#Subroutines#
#############

########## not needed, old code from David Studholme!
sub get_CDSs{
    use Bio::DB::GFF; # a SQL database with stored annotation? Used for implementing relational databases when using bioperl-db (The bioperl-db package is intended to enable the easy access and manipulation of sequences in a relational databases via a Perl interface)
    my $pos = shift or die; # get position from @_
    my $gbrowse_name = shift or die; # get gbrowse name from @_
    my @CDS;

    my $gbrowse_db_name = $gbrowse_name;

    if ($pos =~ m/(\S+):(\d+)\.\.(\d+)/) {

	my ($seq_id, $window_start, $window_end) = ($1, $2, $3);

	warn "Using gbrowse db '$gbrowse_db_name' and gbrowser '$gbrowse_name'\n";
	my $db = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt',
				    -dsn     => "dbi:mysql:$gbrowse_db_name",
				    -user    => 'www-data',
	    );
	# get all feature types in the database
	my %ftypes;
	foreach my $type ($db->types) {
	    warn "Type: $type\n" unless $type =~ /ssaha/i;
	    $ftypes{$type} = 1;
	}



	### Get all CDS in this window
	my @CDS_features  = $db->overlapping_features(-refseq => $seq_id,
						      -start => $window_start,
						      -stop => $window_end,
						      -types => [
							   'CDS',
							   'tRNA',
							   'rRNA',
							   'mRNA',

						      ]);
	foreach my $feature (@CDS_features) {

	    #warn "Analysing feature $feature\n";

	    my $feature_name = $feature->name;
	    my $feature_seq = $feature->seq;
	    my $feature_source = $feature->source;
	    my $feature_method = $feature->method;
	    my $feature_type = $feature->type;

	    $feature_name =~ s/^(\S+)/<b>$1<\/b>/;

	    my %attributes = $feature->attributes();
	    my $locus_tag = $attributes{'Parent'};
	    my $description = $attributes{'Note'};

	    $locus_tag =~ s/Erwinia_toletana\.gene\.//g;

	    my $annotation = "<b>$locus_tag</b> $description";

	    push @CDS, $annotation;
	}
	return @CDS;
    }
}
