#!/usr/bin/perl

#######
# POD #
#######

=pod

=head1 NAME

C<po2anno.pl> - create an annotation comparison matrix from Proteinortho5 output

=head1 SYNOPSIS

C<perl po2anno.pl -i matrix.proteinortho -g genome_fasta_dir/ -l -a E<gt> annotation_comparison.tab>

=head1 DESCRIPTION

Supplement an ortholog/paralog output matrix from a
L<I<Proteinortho5>|http://www.bioinf.uni-leipzig.de/Software/proteinortho/>
calculation with annotation information. The resulting tab-separated
annotation comparison matrix (ACM) is mainly intended for the
transfer of high quality annotations from reference genomes to
homologs (orthologs and co-orthologs/paralogs) in a query genome
(e.g. in conjunction with L<C<tbl2tab.pl>|/tbl2tab>). But of course
it can also be used to have a quick glance at the annotation of
genes present only in a couple of input genomes in comparison to the
others.

Annotation is retrieved from multi-FASTA files created with
L<C<cds_extractor.pl>|/cds_extractor>. See
L<C<cds_extractor.pl>|/cds_extractor> for a description of the
format. These files are used as input for the PO analysis.

I<Proteinortho5> (PO) has to be run with option B<-singles> to include
also genes without orthologs, so-called singletons/ORFans, for each
genome in the PO matrix (see the
L<PO manual|http://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html>).
Additionally, option B<-selfblast> is recommended to enhance paralog
detection by PO.

Each orthologous group (OG) is listed in a row of the resulting ACM,
the first column holds the OG numbers from the PO input matrix (i.e.
line number minus one). The following columns specify the
orthologous CDS for each input genome. For each CDS the ID,
optionally the length in bp (option B<-l>), gene, EC number(s), and
product are shown depending on their presence in the CDS's
annotation. The ID is in most cases the locus tag (see
L<C<cds_extractor.pl>|/cds_extractor>). If several EC numbers exist
for a single CDS they're separated by ';'. If an OG includes
paralogs, i.e. co-orthologs from a single genome, these will be
printed in the following row(s) B<without> a new OG number in the
first column. The order of paralogous CDSs within an OG is
arbitrarily.

The OGs are sorted numerically via the query ID (see option B<-q>).
If option B<-a> is set, the non-query OGs are appended to the output
after the query OGs, sorted numerically via OG number.

=head1 OPTIONS

=head2 Mandatory options

=over 20

=item B<-i>=I<str>, B<-input>=I<str>

Proteinortho (PO) result matrix (*.proteinortho or *.poff), or piped C<STDIN> (-)

=item B<-g>=I<str>, B<-genome_dir>=I<str>

Path to the directory including the genome multi-FASTA PO input
files, created with C<cds_extractor.pl>

=back

=head2 Optional options

=over 20

=item B<-h>, B<-help>

Help (perldoc POD)

=item B<-q>=I<str>, B<-query>=I<str>

Query genome (has to be identical to the string in the PO matrix)
[default = first one in alphabetical order]

=item B<-l>, B<-length>

Include length of each CDS in bp

=item B<-a>, B<-all>

Append non-query orthologous groups (OGs) to the output

=item B<-v>, B<-version>

Print version number to C<STDERR>

=back

=head1 OUTPUT

=over 20

=item C<STDOUT>

The resulting ACM is printed to C<STDOUT>. Redirect or pipe into
another tool as needed (e.g. C<cut>, C<grep>, C<head>, or C<tail>).

=back

=head1 EXAMPLES

=head2 L<C<cds_extractor.pl>|/cds_extractor>

=over

=item C<for i in *.[gbk|embl]; do perl cds_extractor.pl -i $i [-p|-n]; done>

=item C<rename 's/_cds_[aa|nuc].fasta/.[faa|fna]/' *_cds_[aa|nuc].fasta>

=back

=head2 L<I<Proteinortho5>|http://www.bioinf.uni-leipzig.de/Software/proteinortho/>

=over

=item C<proteinortho5.pl -graph [-synteny] -cpus=# -selfblast -singles -identity=50 -cov=50 -blastParameters='-use_sw_tback' *.[faa|fna]>

=back

=head2 C<po2anno.pl>

=over

=item C<perl po2anno.pl -i matrix.[proteinortho|poff] -g genome_fasta_dir/ -q query.[faa|fna] -l -a E<gt> annotation_comparison.tab>

=back

=head1 VERSION

 0.2                                               update: 15-01-2015
 0.1                                                       18-12-2014

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

### Get the options with Getopt::Long
my $Input_File; # PO result matrix (*.proteinortho or *.poff)
my $Genome_Dir; # directory with genome multi-fastas (PO input)
my $Query; # query genome (first column in output)
my $Opt_Length; # include CDS nucleotide lengths in output
my $Opt_All_OGs; # print also all non-query OGs to output
my $VERSION = 0.2;
my ($Opt_Version, $Opt_Help);
GetOptions ('input=s' => \$Input_File,
            'genome_dir=s' => \$Genome_Dir,
            'query=s' => \$Query,
            'length' => \$Opt_Length,
            'all' => \$Opt_All_OGs,
            'version' => \$Opt_Version,
            'help|?' => \$Opt_Help);



### Run perldoc on POD
pod2usage(-verbose => 2) if ($Opt_Help);
die "$0 $VERSION\n" if ($Opt_Version);
if (!$Input_File || !$Genome_Dir) {
    my $warning = "\n### Fatal error: Mandatory options '-i' or '-g' or their arguments are missing!\n";
    pod2usage(-verbose => 1, -message => $warning, -exitval => 2);
}
die "\n### Fatal error: Directory '$Genome_Dir' does not exist: $!\n" if (!-d $Genome_Dir);



### Pipe input from STDIN or open input file
my $Input_Fh;
if ($Input_File eq '-') { # file input via STDIN
    $Input_Fh = *STDIN; # capture typeglob of STDIN
} else { # input via input file
    open ($Input_Fh, "<", "$Input_File");
}



### Parse OGs in input PO matrix
print STDERR "Parsing Proteinortho input matrix ...\n"; # run status of script
my @Genome_Files; # store genome input files
my %Ortho_Groups; # anonymous array in two-dimensional hash to store OGs
while (<$Input_Fh>) {
    chomp;

    # check PO input file header and get genome file names
    if ($. == 1) { # header of PO matrix file (first line)
        die "\n### Fatal error:\nProteinortho input matrix '$Input_File' does not have the mandatory header line, which starts with the first three tab-separated mandatory columns:\n# Species\tGenes\tAlg.-Conn.\n" if (!/^# Species\tGenes\tAlg\.-Conn\./);

        # get input genomes from PO matrix and check $Query existence
        my $query_present = 0 if ($Query);
        foreach (split(/\t/, $_)) {
            next if (/# Species|Genes|Alg\.-Conn\./); # non-genome header columns
            $query_present = 1 if ($Query && $Query eq $_);
             push (@Genome_Files, $_); # store genome in array
        }
        die "\n### Fatal error:\nGiven query '$Query' not found in the header of the input Proteinortho matrix '$Input_File'! Make sure the string corresponds exactly (also lower/uppercase) to the header of the matrix and the input file in '$Genome_Dir'!\n" if ($Query && !$query_present);
        next; # skip to second line of PO matrix
    }

    # parse PO ortholog matrix
    my ($organisms, $genes, $alg_conn, @og) = split(/\t/); # $alg_conn not used
    my $organism_count = 0; # count number of organisms in OG to compare to '$organisms' (should be ok from PO)
    my $gene_count = 0; # count number of genes in OG to compare to '$genes' (should be ok from PO)
    my $column = -1; # column counter to associate column to correct genome file in OG (@Genome_Files array zero-based, thus start with '-1')
    foreach (@og) {
        $column++;
        next if (/^\*$/); # skip empty columns in PO matrix (indicated by '*')
        $organism_count++;
        my @paralogs = split(',');
        foreach (sort {$a cmp $b} @paralogs) {
            $gene_count++;
            push(@{ $Ortho_Groups{$.-1}->{$Genome_Files[$column]} }, $_); # push into anonymous array in two-dimensional hash ($.-1 = OG number, because of header)
        }
    }
    die "\n### Fatal error:\nThe indicated number of species ($organisms) does not fit to the counted number of organisms ($organism_count) in the orthologous group of line '$.' of the Proteinortho input matrix '$Input_File'!\n" if ($organism_count != $organisms); # to check PO matrix
    die "\n### Fatal error:\nThe indicated number of genes ($genes) does not fit to the counted number of orthologs/paralogs ($gene_count) in the orthologous group of line '$.' of the Proteinortho input matrix '$Input_File'!\n" if ($gene_count != $genes); # to check PO matrix
}



### Sort @Genome_Files and set $Query as first element (not possible above because has to follow $column!)
@Genome_Files = sort {lc($a) cmp lc($b)} @Genome_Files; # sort regardless of upper/lowercase
if ($Query) {
    @Genome_Files = grep ($_ ne $Query, @Genome_Files); # remove $Query from @Genome_Files
    unshift (@Genome_Files, $Query); # set as first element
} elsif (!$Query) {
    $Query = $Genome_Files[0]; # $Query first element in sorted array, when not given with option '-q'
}



### Parse annotations in genome multi-fastas
print STDERR "Parsing annotation from multi-FASTA CDS genome files ...\n";
my %Annotation; # two-dimensional hash to store the annotation of the genome files
my %Anno_Features; # two-dimensional hash to store which annotation features are present overall in each individual genome multi-FASTA (optional are 'gene/g=', 'product/p=' and 'EC_number/ec=' tags)
foreach my $genome (@Genome_Files) {
    my $genome_file_path = "./$Genome_Dir/$genome";
    check_file_exist($genome_file_path); # subroutine
    open (my $genome_fh, "<", "$genome_file_path");

    # parse anno in current genome file
    while (my $anno_line = <$genome_fh>) {
        chomp $anno_line;
        die "\n### Fatal error:\n'$genome_file_path' is not a FASTA input file. First line of the file should be a FASTA ID/header line and start with a '>':\n$anno_line\n" if ($anno_line !~ /^>/ && $. == 1);
        next if ($anno_line !~ /^>/); # skip non-ID FASTA lines
        my @anno = split(/\s/, $anno_line);
        my ($id, $gene, $product, $length, $organism, $ec); # $organism not used
        foreach (@anno) {
            if (/^>(.+)$/) {
                $id = $1;
            } elsif (/g=(.*)$/) {
                $gene = $1;
            } elsif (/p=(.*)$/) {
                warn "\n### Warning:\nNo product annotation (p=) in file '$genome_file_path' on line:\n$anno_line\nProceeding ...\n" if (!$1);
                $product = $1;
            } elsif (/l=(.*)$/) {
                $length = $1;
                my ($start, $stop) = split(/\.\./, $length);
                $length = abs($stop - $start) + 1;
            } elsif (/o=(.*)$/) { # $organism not used
                $organism = $1;
            } elsif (/ec=(.*)$/) {
                $ec = $1;
            } else {
                die "\n### Fatal error:\nAnnotation in file '$genome_file_path' is not recognized on the following line. Please use 'cds_extractor.pl' to create your multi-FASTA protein genome files.\n$anno_line\n";
            }
        }
        die "\n### Fatal error:\nThe following FASTA ID of file '$genome_file_path' is not unique but has to be considering ALL input genome files (actually Proteinortho should have complained already). Please modify all repetitive occurences.\n$id\n" if ($Annotation{$id});

        # fill two-dimensional annotation hash
        $Annotation{$id} = {'genome' => $genome,
                            'gene' => $gene,
                            'product' => $product,
                            'length' => $length,
                            'ec' => $ec};

        # store which 'optional' annotation features are present overall in a genome
        # (used to have lots of greps on %Annotation for this purpose below, but lets save performance)
        $Anno_Features{$genome}{'gene'} = 1 if ($gene && !$Anno_Features{$genome}->{'gene'});
        $Anno_Features{$genome}{'product'} = 1 if ($product && !$Anno_Features{$genome}->{'product'});
        $Anno_Features{$genome}{'ec'} = 1 if ($ec && !$Anno_Features{$genome}->{'ec'});
    }
}



### Check if CDS counts in multi-FASTA genome files and PO matrix are equal and print header of output
print STDERR "Comparing annotation CDS counts to orthologous group CDS counts and printing header for output matrix ...\n";
my $Query_CDS_Count; # store number of query CDS for stat report at end
print "# OG"; # first column of header
foreach my $genome (@Genome_Files) {
    my $ortho_cds_count = map ($Ortho_Groups{$_}->{$genome} ? @{ $Ortho_Groups{$_}->{$genome} } : (), keys %Ortho_Groups); # de-reference anonymous array in two-dimensional hash
    # map evaluates BLOCK or EXPR in list context and returns the LIST value composed of the results of each such evaluation (each element of LIST may produce zero, one, or more elements in the returned value <=> grep returns only one value). In scalar context, returns the total number of elements so generated.

    # the map line above is a fancy way of writing
    #foreach (keys %Ortho_Groups) {
        #$ortho_cds_count += @{ $Ortho_Groups{$_}->{$genome} } if ($Ortho_Groups{$_}->{$genome});
    #}

    my $anno_cds_count = grep ($Annotation{$_}->{'genome'} eq $genome, keys %Annotation);
    $Query_CDS_Count = $anno_cds_count if ($Query eq $genome);
    die "\n### Fatal error:\nThere are more CDSs in file '$genome' than for this genome in the Proteinortho matrix '$Input_File', but counts have to be equal. Please run Proteinortho5 with option '-singles' to include also genes without orthologs, so-called singletons/ORFans (recommended is also option '-selfblast' to enhance paralog detection).\n" if ($ortho_cds_count < $anno_cds_count);
    die "\n### Fatal error:\nThere are less CDSs in file '$genome' than for this genome in the Proteinortho matrix '$Input_File', but counts have to be equal. Please check if the files are correct.\n" if ($ortho_cds_count > $anno_cds_count);

    # print header fields for output
    print "\t$genome";
    print "\tlength [bp]" if ($Opt_Length);
    print "\tgene" if ($Anno_Features{$genome}->{'gene'});
    print "\tec" if ($Anno_Features{$genome}->{'ec'});
    print "\tproduct" if ($Anno_Features{$genome}->{'product'}); # very unlikely that there's a whole genome without product annotation ...
}
print "\n";



### Print annotation comparison matrix
print STDERR "Printing output annotation comparison matrix for query OGs ...\n";
my %Query_ID_Seen; # query IDs already processed for output (because of paralogous CDSs in the same OGs)
my %Query_OGs; # OGs containing query CDSs, to skip for non-query OGs (option '-a') below
my $Query_Specific_OGs = 0; # count number of OGs including only query CDSs (not including CDSs of the other genomes)
my $Query_Singletons = 0; # count query singletons/ORFans (different to '$Query_Singleton_OGs' because of paralogous CDSs)
foreach my $id (sort { $a cmp $b } grep ($Annotation{$_}->{'genome'} eq $Query, keys %Annotation)) { # get IDs only for the query genome, sorted
    next if ($Query_ID_Seen{$id});

    # get OG for current query $id (I'm sure this can also be resolved by a fancy map/grep construction, but I didn't get it to work)
    my $id_og; # OG containing the current query ID
    foreach my $og (keys %Ortho_Groups) {
        foreach (@{ $Ortho_Groups{$og}->{$Query} }) {
            if ($_ eq $id) {
                $id_og = $og;
                last;
            }
        }
        last if ($id_og);
    }
    $Query_OGs{$id_og} = 1;

    print_matrix($id_og, \%Query_ID_Seen); # subroutine to print the annotation comparison matrix for the current OG to STDOUT
}



### Optionally, print also non-query OGs/singletons
if ($Opt_All_OGs) {
    print STDERR "Printing output annotation comparison matrix for non-query OGs (option '-all') ...\n";
    foreach my $og (sort {$a <=> $b} keys %Ortho_Groups) { # sort from smallest to largest OG
        next if ($Query_OGs{$og});
        print_matrix($og); # subroutine
    }
}



### Print some final stats
print STDERR "\nFinished creating Proteinortho annotation comparison matrix with query genome '$Query':\n";
print STDERR "Total genomes: ", scalar @Genome_Files, "\n";
print STDERR "Total CDSs: ", scalar keys %Annotation, "\n";
print STDERR "Total query CDSs: $Query_CDS_Count\n";
print STDERR "Total OGs: ", scalar keys %Ortho_Groups, "\n";
print STDERR "OGs including query CDSs: ", scalar keys %Query_OGs, "\n";
print STDERR "Query-specific OGs (not including CDSs of the other genomes): $Query_Specific_OGs\n";
print STDERR "Total query singletons/ORFans: $Query_Singletons\n";


exit;


#############
#Subroutines#
#############

### Subroutine to test for file existence
sub check_file_exist {
    my $file = shift;
    die "\n### Fatal error:\nGenome file '$file' is listed in the Proteinortho matrix, but does not exist in directory '$Genome_Dir': $!\n" if (!-e $file);
    return 1;
}



### Subroutine to print resulting annotation comparison matrix to STDOUT
sub print_matrix {
    my ($og, $query_id_seen) = @_; # $query_id_seen hash-ref to %Query_ID_Seen ($query_id_seen not given/defined if called above for non-query OGs, option '-a')

    # get max count of paralogs/co-orthologs over all genomes in current OG for row printing below
    my $max_paralog_count = 0;
    foreach my $genome (keys %{$Ortho_Groups{$og}}) {
        my $paralog_count = 0;
        foreach (@{ $Ortho_Groups{$og}->{$genome} }) { # only needed for genomes in the CURRENT $og <=> ALL genomes needed below
            $query_id_seen->{$_} = 1 if ($query_id_seen && $genome eq $Query); # store query IDs that have been processed
            $paralog_count++;
        }
        if ($genome eq $Query && keys %{$Ortho_Groups{$og}} == 1) { # query-specific OGs (not including CDSs of the other genomes)
            $Query_Specific_OGs++;
            $Query_Singletons += $paralog_count;
        }
        $max_paralog_count = $paralog_count if ($paralog_count > $max_paralog_count);
    }

    # print output matrix for current OG
    print "$og";
    for (my $i = 0; $i < $max_paralog_count; $i++) { # '<' because array zero-based
        print "\t"; # tab after OG-nr

        foreach my $genome (@Genome_Files) { # for ALL genomes

            # CDS(s) of genome present in current $og and within $max_paralog_count
            if ($Ortho_Groups{$og}->{$genome}->[$i]) {
                print "$Ortho_Groups{$og}->{$genome}->[$i]\t"; # print CDS ID
                print $Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'length'}, "\t" if ($Opt_Length);
                if ($Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'gene'}) {
                    print $Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'gene'}, "\t";
                } elsif ($Anno_Features{$genome}->{'gene'}) { # tab only if this genome has a 'g=/gene' annotation at all (stored in two-dimensional hash %Anno_Features)
                    print "\t";
                }
                if ($Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'ec'}) {
                    print $Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'ec'}, "\t";
                } elsif ($Anno_Features{$genome}->{'ec'}) {
                    print "\t";
                }
                if ($Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'product'}) {
                    print $Annotation{$Ortho_Groups{$og}->{$genome}->[$i]}->{'product'}, "\t";
                } elsif ($Anno_Features{$genome}->{'product'}) {
                    warn "\n### Warning:\nNo product annotation (p=) in file '$genome' for ID:\n$Ortho_Groups{$og}->{$genome}->[$i]\nProceeding ...\n";
                    print "\t";
                }

            } else { # genome without a CDS in current $og or less paralogous CDSs than $max_paralog_count
                my $times = 5;
                $times -= 1 if (!$Opt_Length);
                $times -= 1 if (!$Anno_Features{$genome}->{'gene'}); # genome doesn't have a 'g=/gene' annotation at all
                $times -= 1 if (!$Anno_Features{$genome}->{'ec'});
                $times -= 1 if (!$Anno_Features{$genome}->{'product'}); # very unlikely that there's a whole genome without product annotation ...
                print "\t" x $times;
            }
        }
        print "\n";
    }

    return 1;
}
