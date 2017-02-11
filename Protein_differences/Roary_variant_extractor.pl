###
#	Roary_variant_extractor.pl
#	Script to extract genetic variant sequences from a Roary output
#	Author: Samuel Bloomfield
###


use warnings;
use strict;

#Location of working directory
my $wkdir="/working/directory/location/";

#Name of roary output file
my $input_fasta="Roary_output_example.txt";

#Name of gene sequences as they appear in the roary output file
my @genes = ('Gene_1', 'Gene_2', 'Gene_3');

#Length of flanks
my $flank_length = 500;


my $wkdir_output_1 = ($wkdir . "Gene_sequences/");
my $wkdir_output_2 = ($wkdir . "Gene_sequences_reduced/");


system "echo mkdir $wkdir_output_1";
system "mkdir $wkdir_output_1";

system "echo mkdir $wkdir_output_2";
system "mkdir $wkdir_output_2";


my $gene_numbers = scalar @genes;

my $line;
my $line2;
my $line3;
my $isolate;
my @locations = ();
my @locations_genes = ();
my $first_check;
my $last_check;
my $gene_output;
my $gene_output2;
my @genes_found =();
my @genes_location_found=();
my $size;
my $string;
my $substring;
my $first;
my $last;
my $difference;

my $R_flank;
my $L_flank;
my $gene_size;

my $file1=($wkdir . $input_fasta);
my $file2;
my $outfile;
my $outfile_summary=($wkdir . 'Gene_location_output.txt');
my $outfile_summary_3=($wkdir . 'Gene_location_reduced_output.txt');
my $direction;
my $R_substring;

my $ind;
my $cnt;
my $max;
my $val;

my $sequence;
my @gene_names;
my @multiple_genes;
my @multiple_genes_length;
my $GATC;
my $gene_sequence;
my $GATCN_total;
my $N_start;
my $N_end;
my @start_boundary =();
my @end_boundary = ();
my @specific_names = ();
my @specific_gene_isolate_names = ();

my $loop_variable;
my @gene_length =();



open FILE1, $file1 or die "FILE $file1 NOT FOUND - $!\n";
foreach $_ (<FILE1>) {
	$line = $_;
	@locations = split /\t/, $line;
	@locations_genes = ();
	
	push @locations_genes, @locations[1..$#locations-1];
	
	print $locations[0];
	$file2=($wkdir . $locations[0] . '.gff');
	$outfile=($wkdir . $locations[0] . '.fa');

	#Forms a fasta file from the .gff file

	open FILE2, $file2 or die "FILE $file2 NOT FOUND - $!\n";
	open(OUT, ">$outfile");
	foreach $_ (<FILE2>) {
		$line2 = $_;
		if ($line2=~m/sequence-region\s(.*?)\s\d+\s(\d+)\n/) { 
			print "$1\t$2";
			print OUT ">$1\n";
			$size = $2;
			print $size;
		} elsif ($line2=~/ |\t|#|>/) {
		} elsif ($line2=~/([GATCN]+)/) {
			print OUT "$1";
		}
	}
	close OUT;
	close FILE2;

	
	my $gene_count = 0;
	foreach (@locations_genes) {
		my $loop_variable =$_;
		if ($loop_variable eq "XXXXX") {
			} else {
			open FILE2, $file2 or die "FILE $file2 NOT FOUND - $!\n";
			foreach $_ (<FILE2>) {
				$line2 = $_;
				if ($line2 =~ /($loop_variable)/) {

					#Searches .gff file for gene sequence and records sequence length and flank size

					if ($line2 =~ /CDS\t(\d+)\t(\d+)\t\.\t([\+-])\t/) {
						$first = $1;
						$last = $2;
						$L_flank = $flank_length;
						$R_flank = $flank_length;
						$direction = $3;
						$difference = $last - $first;
						$first_check = $first-$flank_length;
						$last_check = $last+$flank_length+1;
						$gene_size = $last - $first + 1;

						#adjusts flank size if sequence is at the start or end of the genome sequence

						if ($first_check < 1) {
							$first_check = 1;
							$L_flank = $first;
						}
						if ($last_check > $size) {
							$last_check = $size;
							$R_flank = $size - $last;
						}
						open FILE3, $outfile or die "FILE3 $outfile NOT FOUND - $!\n";
						foreach $_ (<FILE3>) {
							$line3 = $_;

							#adjust flank size if the sequence at the start or end of a contig

							if ($line3=~/[GATCN]{50,}/) {
								if ($line3=~/([GATCN]+)/) {
									$string = $1;
									$substring = substr($string, ($first_check - 1), ($last_check - $first_check));
									if ($direction eq "-") {
										$R_substring = reverse($substring);
										$R_substring =~ tr/ACGTN/TGCAN/;
										$substring = $R_substring;
										($L_flank, $R_flank) = ($R_flank, $L_flank);
									}
									$GATCN_total = $substring =~ tr/[GATCN]//;
									$N_start = substr($substring, 0, $L_flank) =~ tr/[N]//;
									$N_end = substr($substring, ($GATCN_total - $R_flank), $R_flank) =~ tr/[N]//;
									push @start_boundary, ($L_flank - $N_start);
									push @end_boundary, ($R_flank - $N_end);
									push @specific_names, $genes[$gene_count];
									push @specific_gene_isolate_names, $loop_variable;
									
									$substring =~ s/N//g;

									push @gene_length, (length $substring);
									
									#Outputs gene sequence in a fasta file

									$gene_output = ($wkdir . 'Gene_sequences/' . 'fasta_' . $genes[$gene_count] . "__" .  $loop_variable . '.fa');
									

									open(OUT2, ">$gene_output") or die "Couldn't open OUT2 $gene_output $!\n";
									print OUT2 ">";
									print OUT2 $genes[$gene_count];
									print OUT2 "__";
									print OUT2 $loop_variable;
									print OUT2 "\n";
									print OUT2 $substring;
								}
							}
						}
					}
				}	
			}
		}
		$gene_count++;
	}
	close FILE2;
}

open(OUT3, ">$outfile_summary") or die "Couldn't open OUT2 $outfile_summary $!\n";
print OUT3 join("\t", @start_boundary);
print OUT3 "\n";
print OUT3 join("\t", @end_boundary);
print OUT3 "\n";
print OUT3 join("\t", @specific_names);
print OUT3 "\n";
print OUT3 join("\t", @specific_gene_isolate_names);
close OUT3;
print OUT3 join("\t", @gene_length);
close OUT3;


my @gene_searcher;
my $gene_files;
my @files;
my @index;
my @different_sequences;
my $location_variable;
my $gene_variable;
my $gene_value;
my $location_name;
my @location_array;
my $gene_print;
my $location_print;
my $gene_index;
my $number;
my $total;

my @final_start_boundary =();
push @final_start_boundary, "Left_flank";
my @final_end_boundary = ();
push @final_end_boundary, "Right_flank";
my @final_specific_names = ();
push @final_specific_names, "Gene_names";
my @final_specific_gene_isolate_names = ();
push @final_specific_gene_isolate_names, "Isolates_gene_names";
my @final_length =();
push @final_length, "Gene_length";

foreach (@genes) {
	$loop_variable = $_;
	$gene_files =($wkdir . "Gene_sequences/fasta_" . $loop_variable . "__");
	@files = glob "$gene_files*.fa";
	@gene_searcher = ();
	@multiple_genes = ();
	@multiple_genes_length = ();
		
	#Places all the sequences of each gene into a single array

	for (0..$#files){
		if ($files[$_]=~/$loop_variable\_\_(.*?)\.fa/) {
			push @gene_searcher, $1;
			open FILE1, $files[$_] or die "FILE $files[$_] NOT FOUND - $!\n";
			foreach $_ (<FILE1>) {
				$line = $_;
				if ($line=~/([GATCN]{20,})/) {
					$sequence = $1;
					push @multiple_genes, $sequence;
					$GATC = $sequence =~ tr/[GATC]//;
					push @multiple_genes_length, $GATC;
				}
				close FILE1;
			}
		}
	}

	#Sorts sequence array by size

	my @index = sort { $multiple_genes_length[$b] <=> $multiple_genes_length[$a] } 0 .. $#multiple_genes_length;

	@gene_searcher = @gene_searcher[@index];
	@multiple_genes = @multiple_genes[@index];
	@multiple_genes_length = @multiple_genes_length[@index];
						
	@different_sequences=();
	@location_array=();

	#Removes sequences in array that are duplicates or contained within a larger sequence
						
	for (0..$#multiple_genes) {
		$location_variable = $multiple_genes[$_];
		$location_name = $gene_searcher[$_];
		if ((scalar @different_sequences) == 0) {
			push @different_sequences, $location_variable;
			push @location_array, $location_name;
		} else {
			$gene_value = 0;
			for (@different_sequences) {
				$gene_variable = $_;
				if ($gene_variable=~/$location_variable/) {
				$gene_value++;
				}
			} 
			if ($gene_value == 0) {
				push @different_sequences, $location_variable;
				push @location_array, $location_name;
			}
		}
	}

	#Places gene variant sequences into a separate directory and records their length and flank lengths

	for (0..$#different_sequences) {
		$gene_print = $different_sequences[$_];
		$location_print = $location_array[$_];
		$number = $_;
		my( $idx )= grep { $specific_gene_isolate_names[$_] eq $location_print } 0..$#specific_gene_isolate_names;
		push @final_start_boundary, $start_boundary[$idx];
		push @final_end_boundary, $end_boundary[$idx];
		push @final_specific_names, $specific_names[$idx];
		push @final_specific_gene_isolate_names, $specific_gene_isolate_names[$idx];
		push @final_length, $gene_length[$idx];
		
		$gene_output2 = ($wkdir . "Gene_sequences_reduced/fasta_" . $loop_variable . "__" .  $location_print . ".fa");
		open(OUT4, ">$gene_output2") or die "Couldn't open OUT4 $gene_output2 $!\n";
		print OUT4 ">";
		print OUT4 $loop_variable;
		print OUT4 "__";
		print OUT4 $location_print;
		print OUT4 "\n";
		print OUT4 $gene_print;
		close OUT4;
	}
}

open(OUT5, ">$outfile_summary_3") or die "Couldn't open OUT5 $outfile_summary_3 $!\n";
print OUT5 join("\t", @final_start_boundary);
print OUT5 "\n";
print OUT5 join("\t", @final_end_boundary);
print OUT5 "\n";
print OUT5 join("\t", @final_specific_names);
print OUT5 "\n";
print OUT5 join("\t", @final_specific_gene_isolate_names);
print OUT5 "\n";
print OUT5 join("\t", @final_length);
	
close OUT5;
