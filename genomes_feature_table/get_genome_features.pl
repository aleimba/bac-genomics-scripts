#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO; # bioperl module to handle sequence input/output
use Bio::SeqFeatureI;  # bioperl module to handle features in a sequence

my $usage = "\n".
   "\t#######################################################\n".
   "\t# $0 file-extension               #\n". #$0 = program name
   "\t#                                                     #\n".
   "\t# Extracts all primary features from all genome files #\n".
   "\t# ((multi)-embl or -genbank, and drafts) with the     #\n".
   "\t# given extension in the current directory and counts #\n".
   "\t# these. Results are written to *_features.txt.       #\n".
   "\t# CAREFUL: Don't put drafts in the same file, if the  #\n".
   "\t# accession-numbers are the same (except for the last #\n".
   "\t# four digits). Otherwise, features will be added for #\n".
   "\t# both draft WGS wrongly to only one name (and the    #\n".
   "\t# other name omitted)!                                #\n".
   "\t# The perl script uses bioperl (www.bioperl.org).     #\n".
   "\t#                                                     #\n".
   "\t# version 0.2 update: 19.09.2012     Andreas Leimbach #\n".
   "\t# 25.11.2011                    aleimba[at]gmx[dot]de #\n".
   "\t#######################################################\n\n";


### Print usage if -h|--h|--help is given as argument or file-extension is not given
my $extension = shift;
if (!defined $extension) {
    die $usage;
} elsif ($extension =~ m/-h/) {
    die $usage;
}


### Save all primary features from all seq-files (to include all possibilities) and count them for each genome/replicon individually (including draft genomes in 'multi'-format)
my $dirname = '.'; # use current directory
my %prim_features; # store all primary features in ALL the seq-files for the print out
my %strain_features; # store all primary features of each strain, counted
my %desc; # hash to save the description (E. coli ...) of each strain
my %length; # store the sequence length for each strain
my %gc; # store the gc_content for each strain
my %unresolved_n; # store unresolved bases ('n/N's) for each strain
my %coding_percentage; # save the coding percentage for each strain
my $file_count = 0; # drafts have very similar acc#s if the last four digits are removed (see below), but they should be in different files
opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
while (defined(my $file = readdir(DIR))) {
    if ($file =~ m/.+\.$extension$/)  {
	$file_count++;
	my $seqio_obj = Bio::SeqIO->new(-file => "<$file");
	my %draft; # count the contigs/scaffolds for each strain
	my %draft_seq; # store the concatenated seq for all draft contigs/scaffolds for subroutine gc_content for each strain
	my %coding_bases; # save base count to calculate coding percentage for each strain
	while (my $seq_obj = $seqio_obj->next_seq) { # multi-seq files possible (e.g. chromosome plus plasmids, and also drafts with several contigs/scaffolds)
	    my $acc = $seq_obj->accession_number;
	    if ($seq_obj->keywords =~ /WGS/) { # look for draft sequences
		$acc =~ s/^(\S+)(\d{4})$/$1/; # get rid of the last four numbers to include all contigs for drafts (hopefully the contig count doesn't reach the ten thousand mark ...)
		$acc = $acc.$file_count; # append the file_count to the acc# to make it unique for several draft files!
 	    	$draft{$acc}++; # count contigs/scaffolds
		$prim_features{'contigs/scaffolds'} = 1; # include in the header of the result file
		$strain_features{$acc}{'contigs/scaffolds'} = $draft{$acc};
		if ($draft{$acc} == 1) { # give a the range of acc#s for drafts
		    $strain_features{$acc}{'acc_start'} = $seq_obj->accession_number;
		} else {
		    $strain_features{$acc}{'acc_stop'} = $seq_obj->accession_number;
		}
	    }
#	    my $species_obj =  $seq_obj->species; # doesn't work, because '$species_obj->sub_species' gives only 'str.' for MG1655 
	    $desc{$acc} = $seq_obj->desc; # for drafts only the last sequence object
	    $desc{$acc} =~ s/(,* complete genome.*|DNA, complete genome.*|chromosome, complete sequence.*|complete genome, strain.*|, complete sequence.*)//; # shorten strain descriptions
	    if (defined($draft{$acc}) && $draft{$acc} >= 1) { # draft seqs only have the shortened acc now
		if ($draft{$acc} == 1) { # the first contig/scaffold of drafts
		    print "Processing draft genome: ", $desc{$acc}, ", ", $seq_obj->accession_number, "\n"; # print only for the first contig
		    $length{$acc} = $seq_obj->length;
		    ($gc{$acc}, $unresolved_n{$acc}) = gc_content($seq_obj->seq); # subroutine to calculate the GC-content and count 'N's in the sequence
		    $draft_seq{$acc} = $seq_obj->seq; # save the sequence to concatenate to the following contigs/scaffolds
		    $coding_bases{$acc} = 0; # set to zero for the first contig of the draft
		} elsif ($draft{$acc} > 1) { # calculations needed for all further contigs/scaffolds
		    $length{$acc} = $seq_obj->length + $length{$acc}; # add all previous lengths to the current
		    $draft_seq{$acc} = $seq_obj->seq . $draft_seq{$acc}; # concatenate all contig/scaffold sequences for one draft genome
		    my $prev_n = $unresolved_n{$acc};
		    ($gc{$acc}, $unresolved_n{$acc}) = gc_content($draft_seq{$acc}); # use the concatenated sequence for subroutine
		    $unresolved_n{$acc} = $unresolved_n{$acc} + $prev_n;
		}
	    } else {
		print "Processing complete replicon: ", $desc{$acc}, ", $acc\n";
		$length{$acc} = $seq_obj->length;
		($gc{$acc}, $unresolved_n{$acc}) = gc_content($seq_obj->seq);
		$coding_bases{$acc} = 0; # set to zero for non-draft sequences
	    }
	    foreach my $feat_obj ($seq_obj->get_SeqFeatures) {
		if ($feat_obj->primary_tag eq 'source') { # exclude source primary tag
		    next;
		} elsif ($feat_obj->primary_tag eq 'gene') {
		    $prim_features{$feat_obj->primary_tag} = 1; # NOT really needed as always excluded below
		    $strain_features{$acc}{$feat_obj->primary_tag}++;
		    foreach my $tag ($feat_obj->get_all_tags) {
			if ($tag eq 'pseudo') {
			    $strain_features{$acc}{'pseudo_gene'}++;
			}
		    }
		} elsif ($feat_obj->primary_tag eq 'CDS') {
		    $prim_features{$feat_obj->primary_tag} = 1; # NOT really needed as always excluded below
		    $strain_features{$acc}{$feat_obj->primary_tag}++;
		    my $pseudo = 0; # pseudo switch to exclude pseudo locations from the coding percentage
		    foreach my $tag ($feat_obj->get_all_tags) {
			if ($tag eq 'pseudo') {
			    $strain_features{$acc}{'pseudo_CDS'}++;
			    $pseudo = 1;
			}
		    }
		    if ($pseudo == 0) { # Don't include pseudogenes in coding percentage
			$coding_bases{$acc} = ($feat_obj->location->end - $feat_obj->location->start) + 1 + $coding_bases{$acc};
		    }
		} else { # get ALL the other primary tags
		    $prim_features{$feat_obj->primary_tag} = 1;
		    $strain_features{$acc}{$feat_obj->primary_tag}++;
		}
	    }
	    $coding_percentage{$acc} = ($coding_bases{$acc}/$length{$acc})*100;
	    $coding_percentage{$acc} = sprintf("%.2f", $coding_percentage{$acc}); # round percentage to two decimals
	}
    }
}
closedir DIR;


### Print header of output with all possible primary features
my $temp = "temp.txt"; # use a temp file, because can't get rid of the trailing \t's otherwise (see below)
open(TEMP, ">$temp") or die "Failed to create file \'$temp\': $!\n";
print TEMP "Name (last contig/scaffold for drafts)\tSize\tGC-content\tCoding percentage\tTotal CDS (pseudo)\tTotal genes (pseudo)\trRNA\ttRNA\ttmRNA\tncRNA\tAccession (start..stop for drafts)\tcontigs/scaffolds\tUnresolved bases (Ns)\t"; # specific print order for the most common primary features
foreach (sort keys %prim_features) { # print the residual primary features
    if ($_ =~ m/^(CDS|gene|rRNA|tRNA|tmRNA|ncRNA|contigs\/scaffolds)$/) { # exclude the print order ones
	next;
    } else {
	print TEMP "$_\t";
    }
}
print TEMP "\n";
print TEMP "? = tag, key or qualifier not existent\n"; # explanation line


### Print the primary features for each strain, print '?' if key is not used in a seq-file
foreach my $acc (sort keys %strain_features) {
    print TEMP $desc{$acc}, "\t";
    print TEMP $length{$acc}, "\t";
    print TEMP $gc{$acc}, "\t";
    if ($coding_percentage{$acc} > 0) { # Some genbanks have only 'gene' tags, no 'CDS's
	print TEMP $coding_percentage{$acc}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'CDS'}) { # Some genbanks have only 'gene' tags, no 'CDS's
	print TEMP $strain_features{$acc}{'CDS'};
    } else {
	print TEMP "?";
    }
    if (defined $strain_features{$acc}{'pseudo_CDS'}) {
	print TEMP  " (", $strain_features{$acc}{'pseudo_CDS'}, ")\t";
    } else {
	print TEMP " (?)\t";
    }
    if (defined $strain_features{$acc}{'gene'}) {
	print TEMP $strain_features{$acc}{'gene'};
	if (defined $strain_features{$acc}{'pseudo_gene'}) {
	    print TEMP  " (", $strain_features{$acc}{'pseudo_gene'}, ")\t";
	} else {
	    print TEMP " (?)\t";
	}
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'rRNA'}) {
	print TEMP $strain_features{$acc}{'rRNA'}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'tRNA'}) {
	print TEMP $strain_features{$acc}{'tRNA'}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'tmRNA'}) {
	print TEMP $strain_features{$acc}{'tmRNA'}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'ncRNA'}) {
	print TEMP $strain_features{$acc}{'ncRNA'}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $strain_features{$acc}{'acc_start'}) { # for drafts first and last acc#
	print TEMP $strain_features{$acc}{'acc_start'}, "..", $strain_features{$acc}{'acc_stop'}, "\t";
    } else {
	print TEMP $acc, "\t";
    }
    if (defined $strain_features{$acc}{'contigs/scaffolds'}) {
	print TEMP $strain_features{$acc}{'contigs/scaffolds'}, "\t";
    } else {
	print TEMP "?\t";
    }
    if (defined $unresolved_n{$acc}) {
	print TEMP $unresolved_n{$acc}, "\t";
    } else {
	print TEMP "?\t";
    }
    foreach (sort keys %prim_features) { # print the residual primary features
	if ($_ =~ m/^(CDS|gene|rRNA|tRNA|tmRNA|ncRNA|contigs\/scaffolds)$/) { # exclude the above ones
	    next;
	} elsif (defined $strain_features{$acc}{$_}) {
	    print TEMP $strain_features{$acc}{$_}, "\t";
	} else {
	    print TEMP "?\t";
	    next;
	}
    }
    print TEMP "\n";
}
close TEMP;


### Get rid of the trailing '\t's at each line, couldn't find another solution
open (FILE, "<$temp");
my $output = "$extension\_features.txt";
if (-e $output) {
    print "\n###The result file $output exists already and will be replaced!";
}
open(OUT, ">$output") or die "Failed to create file \'$output\': $!\n";
while (<FILE>) {
    chomp;
    $_ =~ s/\t$//;
    print OUT "$_\n";
}
close FILE;
unlink $temp or warn "Could not delete file \'$temp\': $!";
close OUT;
print "\n###Result file \'$extension\_features.txt\' was created!\n\n";


exit;


#############
#Subroutines#
#############

### Subroutine to calculate gc-content
sub gc_content {
    my $seq = shift;
    my $A = ($seq =~ tr/[aA]//); # transliterations don't accept modifiers like case-insensitive 'i'
    my $C = ($seq =~ tr/[cC]//);
    my $G = ($seq =~ tr/[gG]//);
    my $T = ($seq =~ tr/[tT]//);
    my $N = ($seq =~ tr/[nN]//);
    my $gc_content = (($C + $G)/($A + $C + $G + $T))*100;
    $gc_content = sprintf("%.2f", $gc_content); # round percentage to two decimals
    return ($gc_content, $N);
}
