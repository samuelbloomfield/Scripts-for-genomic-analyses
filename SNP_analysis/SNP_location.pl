###
#	SNP_location.pl
#	Version 2.0
#	Script to locate the position and types of SNPs in a .gff file
#	Author: Samuel Bloomfield
###


use warnings;
use strict;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

my $GFF;
my $SNP;

GetOptions ('GFF:s'		=> \$GFF,
			'SNP:s'	=> \$SNP);




#Creates output files from GFF file name or uses 'SNP' as the base name
my $name;

if ($GFF=~m/(.*?)\.gff/) {
    $name = $1;
} else {
    $name = 'SNP';
}    

my $output = ($name . '_SNP_types.txt');
my $gene_output = ($name . '_gene_SNP_counts.txt');
my $SNP_summary = ($name . '_SNP_summary.txt');

my $line;
my $first;
my $second;
my $orientation;
my $protein_counter;

my @count = ();
my @codon=();
my @codon_sequence=();
my @lines = ();

my $x;
my $y;
my $z;
my $position;

my $line_2;
my $size;
my @codon_variant=();
my @sequences = ();
my $sequence;


my @DNA_sequences=();
my $temp_DNA;

my $temp1;
my $temp2;

my $intergenic = 0;
my $single_gene = 0;
my $multiple_gene = 0;

my $synonymous = 0;
my $non_synonymous = 0;

open FILE, $GFF or die "GFF $GFF NOT FOUND - $!\n";

#Reads sequence contained within .gff file and compiles into a single string

foreach $_ (<FILE>) {
	$line_2 = $_;
	if ($line_2=~m/sequence-region\s(.*?)\s\d+\s(\d+)\n/) { 
		print "$1\t$2\n";
		$size = $2;
	} elsif ($line_2=~/ |\t|#|>/) {
	} elsif ($line_2=~/([GATCN]+)/) {
		push @sequences, $1;
	}
}


close FILE;

$sequence = join("", @sequences);

my $SNP_line;

my $loc;
my @locations=();
my $pos;
my $base;
my @bases;
my @SNP_info;
my $bases_raw;

my @codon_variant_line;
my @codon_line;

my @protein_variant;

my $amino_acid;
my $aa_count;
my $aa_type;
my @SNP_type_line;
my @SNP_type=();

my $bio_sequence;
my $protein_sequence;
my $codon_var;
my $base_length;

my $counter = 0;

#Reads information provided on the SNPs in the SNP file

open FILE2, $SNP or die "FILE $SNP NOT FOUND - $!\n";

foreach $_ (<FILE2>) {
	$SNP_line = $_;
	

	@SNP_info = split (/\s/, $SNP_line);
	
	$loc = $SNP_info[0];
	$pos = $SNP_info[1];
	
    push @locations, $loc;
	
	$base_length = scalar @SNP_info - 2;
	
	@bases= splice (@SNP_info, 2, $base_length);
	
	@bases = grep { $_ ne '' } @bases;

	$protein_counter=0;
	
	
	@codon_line=();
	@SNP_type_line=();
	
	open FILE1, $GFF or die "FILE $GFF NOT FOUND - $!\n";
	foreach $_ (<FILE1>) {
		$line = $_;

		
		#Identifies what genes, if any, the SNP is associated with

		if($line=~m/CDS\t(\d+)\t(\d+)\t\.\t([+-])\t/) {
			$first = $1;
			$second = $2;
			$orientation = $3;
			
			if($loc>($first-1) and $loc<($second+1)){
				$protein_counter++;
				
				if ($orientation eq '+') {
					$temp_DNA = substr($sequence, ($first -1), ($second - $first + 1));
					push @DNA_sequences, $temp_DNA;
				} elsif ($orientation eq '-') {
					$temp_DNA = substr($sequence, ($first -1), ($second - $first + 1));
					$temp_DNA =~ tr/ACGTN/TGCAN/;
					$temp_DNA = reverse $temp_DNA;
					push @DNA_sequences, $temp_DNA;
				}

                @protein_variant=();
                @codon_variant_line=();
                
				#Identifed the codon position of the SNP and extracts the codon
				
				foreach(@bases){
					$base = $_;
					
					
				
					if($orientation eq '+') {
						if($pos eq "R") {
							$base =~ tr/ACGTN/TGCAN/;
						}
						
						$x = $loc - $first;
						$y = 0;
						$z = ($second - $first);
						$position = (($x%3)+1);
						push @codon_line, $position;
						if($position==1){
							push @codon_variant_line, ($base . (substr($sequence, ($loc), 2)));
						} elsif ($position==2){
							push @codon_variant_line, ((substr($sequence, ($loc -2), 1)) . $base . (substr($sequence, ($loc), 1)));
						} elsif ($position==3){
							push @codon_variant_line, ((substr($sequence, ($loc -3), 2)) . $base);
						}
					} elsif($orientation eq '-') {
						$base =~ tr/ACGTN/TGCAN/;
						if($pos eq "R") {
							$base =~ tr/ACGTN/TGCAN/;
						}

						$x = $second - $loc;
						$y = $second - $first;
						$z = 0;
						$position = (($x%3)+1);
						push @codon_line, $position;
							if($position==1){
							$temp1= substr($sequence, ($loc -3), 1);
							$temp2= substr($sequence, ($loc -2), 1);
							$temp1=~ tr/ACGTN/TGCAN/;
							$temp2=~ tr/ACGTN/TGCAN/;
							push @codon_variant_line, ($base . $temp2 . $temp1);
						} elsif ($position==2){
							$temp1= substr($sequence, ($loc), 1);
							$temp2= substr($sequence, ($loc -2), 1);
							$temp1=~ tr/ACGTN/TGCAN/;
							$temp2=~ tr/ACGTN/TGCAN/;
							push @codon_variant_line, ($temp1 . $base . $temp2);
						} elsif ($position==3) {
							$temp1= substr($sequence, ($loc), 1);
							$temp2= substr($sequence, ($loc + 1), 1);
							$temp1=~ tr/ACGTN/TGCAN/;
							$temp2=~ tr/ACGTN/TGCAN/;
							push @codon_variant_line, ($temp2 . $temp1 . $base);
						}
					}
      
				}
			
			    #Converts the codon sequence into an amino acid
	
                foreach(@codon_variant_line) {
	                $codon_var = $_;

	                $bio_sequence = Bio::Seq->new(-seq => $codon_var, -alphabet =>'dna');

	                $protein_sequence = $bio_sequence->translate;
		
	                push @protein_variant, $protein_sequence->seq;
			
                }
		

                #Determines whether the SNP is synonymous or non-synonymous
		
                if(scalar(@protein_variant) < 2){
	                push @SNP_type_line, "Synonymous";
                } else {
	                $amino_acid = $protein_variant[0];
	                $aa_count=0;
	                for(my $i=0; $i < scalar(@protein_variant); $i++) {
		                if($protein_variant[$i] eq $amino_acid){
		                } else {
			                $aa_count++
		                }
	                }
			
	                if ($aa_count > 0){
		                push @SNP_type_line, "Non-synonymous";
	                } else {
		                push @SNP_type_line, "Synonymous";
	                }
                }
            }
    	}
	}
	
	push @count, $protein_counter;
	
	if ($protein_counter == 0) {
		push @codon_variant, "NA";
		push @codon, "NA";
		push @lines, "NA";
		push @SNP_type, "Synonymous";
		$intergenic++;
		$synonymous++;
	} elsif($protein_counter == 1) {
		push @codon_variant, $codon_variant_line[0];
		push @codon, $codon_line[0];
		push @lines, $lines[0];
		push @SNP_type, $SNP_type_line[0];
		if ($SNP_type_line[0] eq 'Non-synonymous'){
		    $non_synonymous++;
		} elsif ($SNP_type_line[0] eq 'Synonymous') {
		    $synonymous++;
		}
		$single_gene++;
	} elsif($protein_counter > 1) {
		push @codon_variant, "NA";
		push @codon, "NA";
		push @lines, "NA";
		$multiple_gene++;
		$aa_type = 0;
		foreach(@SNP_type_line){
			if($_ eq 'Non-synonymous'){
				$aa_type++;
			}
		}
		if($aa_type > 0){
			push @SNP_type, 'Non-synonymous';
			$non_synonymous++;
		} else {
			push @SNP_type, 'Synonymous';
			$synonymous++;
		}
	}
	
	$counter++;
	
	print ("Analysed " . $counter . " SNPs\n");
}


foreach(@SNP_type_line){
    print ($_ . "\n");
}


close FILE1;

open(OUT, ">$output") or die "Couldn't open OUT $output $!\n";
print OUT "Locations\t";
print OUT join ("\t", @locations);
print OUT "\n";
print OUT "Gene_count\t";
print OUT join("\t", @count);
print OUT "\n";
print OUT "Codon_position\t";
print OUT join("\t", @codon);
print OUT "\n";
#print OUT "Codon_variant\t";
#print OUT join("\t", @codon_variant);
#print OUT "\n";
print OUT "SNP_type\t";
print OUT join("\t", @SNP_type);
print OUT "\n";

close OUT;


open(OUT2, ">$SNP_summary") or die "Couldn't open OUT $SNP_summary $!\n";
print OUT2 ("Total number of SNPs analysed: " . scalar(@locations) . "\n\n");
print OUT2 ("SNP types:\n");
print OUT2 ("Intergenic\tSingle_genes\tMultiple_genes\tSynonymous\tNon-synonymous\n");
print OUT2 ($intergenic . "\t" . $single_gene . "\t" .$multiple_gene . "\t" . $synonymous . "\t" . $non_synonymous);

close OUT2;

print ("\nIdentifying the number of SNPs contained in each gene\n");

my @unique_DNA_sequences = uniq(@DNA_sequences);
my $unique_gene;
my $gene_SNP_count;

open(OUT3, ">$gene_output") or die "Couldn't open OUT $gene_output $!\n";

print OUT3 ("Number_of_SNPs\tGene_sequence\n");

foreach(@unique_DNA_sequences){
    $unique_gene = $_;
    $gene_SNP_count = grep { $_ eq $unique_gene  } @DNA_sequences;
    print OUT3 ($gene_SNP_count . "\t" . $unique_gene . "\n");
}

close OUT3;
