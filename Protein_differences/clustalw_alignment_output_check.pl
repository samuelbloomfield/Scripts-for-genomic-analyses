###
#	clustalw_alignemnt_output_check.pl
#	Script to check clustalw alignments
#	Author: Samuel Bloomfield
###


use warnings;
use strict;

#Array of specific gene sequence alignments
my @genes = ('pan_genome_results_Gene_1',
'pan_genome_results_Gene_2',
'pan_genome_results_Gene_3');

#Pathway to directory containing alignment output files
my $alignment_folder = "/direction/to/directory containing/alignment/files/";

#Location of working directory
my $wkdir="/working/directory/location/";

#Number of isolates analysed
my $isolate_number = 3;

my $output = ($wkdir . "alignment_summary.txt");
my $file1;
my $input;
my $line;
my $check;
my @align = ();
my @no_align = ();
my $isolate_check;
my $score_check;
my $number;
my $score_number;

foreach (@genes) {
	$input = $_;
	$file1 = ($alignment_folder . $input . ".output");
	$isolate_check = "Yes";
	$score_check = "Yes";
	
	open FILE1, $file1 or die "FILE $file1 NOT FOUND - $!\n";
	
	#Checks whether there are any sequences that do not perfectly align

	foreach $_ (<FILE1>) {
		$line = $_;
		if ($line =~/Score:\s+(\d+)/) {
			$score_number = $1;
			if ($score_number < 100) {
				$score_check = "No";
			}
		}

		#Checks whether there are any sequences missing

		if ($line =~/There\s+are\s+(\d+)\s+groups/) {
			$number = $1;
			if ($number < ($isolate_number - 1)) {
				$isolate_check = "No";
			}
		}
	}
	
	if ($isolate_check eq "No") {
		push @no_align, $input;
	} elsif ($score_check eq "No") {
		push @no_align, $input;
	} else {
		push @align, $input;
	}
	close FILE1;
}

open(OUT, ">$output") or die "Couldn't open OUT $output $!\n";
print OUT "Protein sequences that align perfectly:\n";
print OUT join("\t", @align);
print OUT "\nProtein sequences that do not align perfectly:\n";
print OUT join("\t", @no_align);
close OUT;



