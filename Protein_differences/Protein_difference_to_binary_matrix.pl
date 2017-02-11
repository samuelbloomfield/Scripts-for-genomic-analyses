###
#	Protein_difference_to_binary_matrix.pl
#	Script to convert Protein differences into a binary matrix
#	Author: Samuel Bloomfield
###


use warnings;
use strict;


#Location of working directory
my $wkdir="/working/directory/location/";

#Location of protein differences file
my $input="Protein_differences.txt";

#Array of isolates with comma after isolate name
my @isolates = ('Isolate_1,', 'Isolate_2,', 'Isoalte_3,');

my $file1=($wkdir . $input);
my $outfile=($wkdir . 'binary.txt');

my @presence;
my $line;
my $isolate;
my $check;


open(OUT2, ">$outfile") or die "Couldn't open OUT2 $outfile $!\n";
print OUT2 join("\t", @isolates);
print OUT2 "\n";


#Checks whether each isolate contains the protein difference
open FILE1, $file1 or die "FILE $file1 NOT FOUND - $!\n";
foreach $_ (<FILE1>) {
	$line = $_;
	@presence = ();
	
	foreach (@isolates) {
		$isolate = $_;
		if ($line=~/$isolate/){
			$check = 1;
		} else {
			$check = 0;
		}
		push @presence, $check;
	}
	print OUT2 join("\t", @presence);
	print OUT2 "\n";
}

close OUT2;	
		
	
