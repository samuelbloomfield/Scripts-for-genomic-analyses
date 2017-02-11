###
#	SNP_location.pl
#	Script to locate the position and types of SNPs in a .gff file
#	Author: Samuel Bloomfield
###


use warnings;
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

#Location of working directory
my $wkdir="/working/directory/location/";

#.gff file name
my $input_gff="Isolate.gff";

#Location of SNPs within genome
my @locations = (1001, 1234, 5002);

#Direction of SNPs
#Note that for the snippy all the directions will be forward (F), whilst for kSNP they can be forward(F) or reverse (R)
my @position = ("F",	"F",	"R");

#First SNP variant
my @SNP_variant_1 =("G",	"C",	"A");

#Second SNP variant
my @SNP_variant_2 =("T",	"T",	"G");

my @codon_variant_1=();
my @codon_variant_2=();

my $file1=($wkdir . $input_gff);
my $output=($file1 . '_SNP_locations');

my $line_no;

my $loc;
my $line;
my $first;
my $second;
my $orientation;
my $counter;
my $protein_counter;

my @count = ();
my @loc_found = ();
my @codon=();
my @codon_sequence=();
my @lines = ();

my $x;
my $y;
my $z;
my $position;

my $line_2;
my $size;

my @sequences = ();
my $sequence;

my $temp_sequence_1;
my $temp_sequence_2;

my @DNA_sequences=();
my $temp_DNA;

my $temp1;
my $temp2;

open FILE, $file1 or die "FILE $file1 NOT FOUND - $!\n";

#Reads sequence data from a .gff file

foreach $_ (<FILE>) {
	$line_2 = $_;
	if ($line_2=~m/sequence-region\s(.*?)\s\d+\s(\d+)\n/) { 
		print "$1\t$2";
		$size = $2;
	} elsif ($line_2=~/ |\t|#|>/) {
	} elsif ($line_2=~/([GATCN]+)/) {
		push @sequences, $1;
	}
}


close FILE;

$sequence = join("", @sequences);

$counter=0;

foreach (@locations) {
	$loc = $_;
	$protein_counter=0;
	
	$line_no = 0;
	
	open FILE1, $file1 or die "FILE $file1 NOT FOUND - $!\n";
	foreach $_ (<FILE1>) {
		$line = $_;
		
		$line_no++;
		
		#Identifies what genes, if any, the SNP is associated with

		if($line=~m/CDS\t(\d+)\t(\d+)\t\.\t([+-])\t\d+\tID/) {
			$first = $1;
			$second = $2;
			$orientation = $3;
			
			if($loc>($first-1) and $loc<($second+1)){
				push @loc_found, $loc;
				$protein_counter++;
				
				push @lines, $line_no;
				
				if ($orientation eq '+') {
					$temp_DNA = substr($sequence, ($first -1), ($second - $first + 1));
					push @DNA_sequences, $temp_DNA;
				} elsif ($orientation eq '-') {
					$temp_DNA = substr($sequence, ($first -1), ($second - $first + 1));
					$temp_DNA =~ tr/ACGTN/TGCAN/;
					$temp_DNA = reverse $temp_DNA;
					push @DNA_sequences, $temp_DNA;
				}

				#Identifed the codon position of the SNP and extracts the codon
				
				if($orientation eq '+') {
					$temp_sequence_1 = $SNP_variant_1[$counter];
					$temp_sequence_2 = $SNP_variant_2[$counter];
					if($position[$counter] eq "R") {
						$temp_sequence_1 =~ tr/ACGTN/TGCAN/;
						$temp_sequence_2 =~ tr/ACGTN/TGCAN/;
					}
					$x = $loc - $first;
					$y = 0;
					$z = ($second - $first);
					$position = (($x%3)+1);
					push @codon, $position;
					if($position==1){
						push @codon_variant_1, ($temp_sequence_1 . (substr($sequence, ($loc), 2)));
						push @codon_variant_2, ($temp_sequence_2 . (substr($sequence, ($loc), 2)));
					} elsif ($position==2){
						push @codon_variant_1, ((substr($sequence, ($loc -2), 1)) . $temp_sequence_1 . (substr($sequence, ($loc), 1)));
						push @codon_variant_2, ((substr($sequence, ($loc -2), 1)) . $temp_sequence_2 . (substr($sequence, ($loc), 1)));
					} elsif ($position==3){
						push @codon_variant_1, ((substr($sequence, ($loc -3), 2)) . $temp_sequence_1);
						push @codon_variant_2, ((substr($sequence, ($loc -3), 2)) . $temp_sequence_2);
					}
				} elsif($orientation eq '-') {
					$temp_sequence_1 = $SNP_variant_1[$counter];
					$temp_sequence_2 = $SNP_variant_2[$counter];
					if($position[$counter] eq "R") {
						$temp_sequence_1 =~ tr/ACGTN/TGCAN/;
						$temp_sequence_2 =~ tr/ACGTN/TGCAN/;
					}
					$temp_sequence_1 =~ tr/ACGTN/TGCAN/;
					$temp_sequence_2 =~ tr/ACGTN/TGCAN/;
					$x = $second - $loc;
					$y = $second - $first;
					$z = 0;
					$position = (($x%3)+1);
					push @codon, $position;
					if($position==1){
						$temp1= substr($sequence, ($loc -3), 1);
						$temp2= substr($sequence, ($loc -2), 1);
						$temp1=~ tr/ACGTN/TGCAN/;
						$temp2=~ tr/ACGTN/TGCAN/;
						push @codon_variant_1, ($temp_sequence_1 . $temp2 . $temp1);
						push @codon_variant_2, ($temp_sequence_2 . $temp2 . $temp1);
					} elsif ($position==2){
						$temp1= substr($sequence, ($loc), 1);
						$temp2= substr($sequence, ($loc -2), 1);
						$temp1=~ tr/ACGTN/TGCAN/;
						$temp2=~ tr/ACGTN/TGCAN/;
						push @codon_variant_1, ($temp1 . $temp_sequence_1 . $temp2);
						push @codon_variant_2, ($temp1 . $temp_sequence_2 . $temp2);
					} elsif ($position==3) {
						$temp1= substr($sequence, ($loc), 1);
						$temp2= substr($sequence, ($loc + 1), 1);
						$temp1=~ tr/ACGTN/TGCAN/;
						$temp2=~ tr/ACGTN/TGCAN/;
						push @codon_variant_1, ($temp2 . $temp1 . $temp_sequence_1);
						push @codon_variant_2, ($temp2 . $temp1 . $temp_sequence_2);
					}
				}
			}
		}
	}
	push @count, $protein_counter;
	
	if ($protein_counter == 0) {
		push @codon_variant_1, "NA";
		push @codon_variant_2, "NA";
		push @codon, "NA";
		push @loc_found, "NA";
		push @DNA_sequences, "NA";
		push @lines, "NA";
	}
	$counter++;
}

my @protein_variant_1=();
my @protein_variant_2=();

my @SNP_type=();

my $bio_sequence;
my $protein_sequence;

#Converts the codon sequence into an amino acid

for(my $i=0; $i < scalar(@codon_variant_1); $i++) {

	if($codon_variant_1[$i] eq 'NA') {
		push @protein_variant_1, 'NA';
	} else {
		$bio_sequence = Bio::Seq->new(-seq => $codon_variant_1[$i], -alphabet =>'dna');

		$protein_sequence = $bio_sequence->translate;
	
		push @protein_variant_1, $protein_sequence->seq;
	
	}

	if($codon_variant_2[$i] eq 'NA') {
		push @protein_variant_2, 'NA';
	} else {
		$bio_sequence = Bio::Seq->new(-seq => $codon_variant_2[$i], -alphabet =>'dna');

		$protein_sequence = $bio_sequence->translate;
	
		push @protein_variant_2, $protein_sequence->seq;
	
	}
}


#Determines whether the SNP is synonymous or non-synonymous

for(my $i=0; $i < scalar(@protein_variant_1); $i++) {
	if($protein_variant_1[$i] eq $protein_variant_2[$i]){
		push @SNP_type, "Synonymous";
	} else {
		push @SNP_type, "Non-synonymous";
	}
}

close FILE1;

open(OUT, ">$output") or die "Couldn't open OUT $output $!\n";
print OUT "Locations\t";
print OUT join ("\t", @locations);
print OUT "\n";
print OUT "Gene_count\t";
print OUT join("\t", @count);
print OUT "\n";
#print OUT join("\t", @loc_found);
#print OUT "\n";
print OUT "Codon_position\t";
print OUT join("\t", @codon);
print OUT "\n";
print OUT "Codon_variant_1\t";
print OUT join("\t", @codon_variant_1);
print OUT "\n";
print OUT "Codon_variant_2\t";
print OUT join("\t", @codon_variant_2);
print OUT "\n";
print OUT "Protein_variant_1\t";
print OUT join("\t", @protein_variant_1);
print OUT "\n";
print OUT "Protein_variant_2\t";
print OUT join("\t", @protein_variant_2);
print OUT "\n";
print OUT "SNP_type\t";
print OUT join("\t", @SNP_type);
print OUT "\n";
#print OUT join("\t", @DNA_sequences);
#print OUT "\n";
#print OUT join("\t", @lines);
#print OUT "\n";



close OUT;
