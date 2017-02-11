###
#	Flank_remover.pl
#	Script to convert .cns files into amino acid sequences and remove flanks
#	Author: Samuel Bloomfield
###


use warnings;
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

#Location of working directory
my $wkdir="/working/directory/location/";

#Array of isolates
my @isolates = ('Isolate_1', 'Isolate_2', 'Isoalte_3');


#Array of specific gene sequences
my @isolates = ('Gene_1__Isolate_1','Gene_1__Isolate_2', 'Gene_2__Isolate_2', 'Gene_2__Isolate_3', 'Gene_3__Isolate_3');

#Array of gene sizes
my @gene_size = (1005, 1006, 1932, 2018, 576);

#Array of left flank lengths
my @left_bound = (500, 500, 452, 500, 450);

#Array of right flank lengths
my @right_bound = (500, 490, 500, 500, 36);


my $gene;
my $line;
my $location;
my $consensus;
my $alteration;
my @sequence;
my $right_border;
my $gene_sequence;
my $alteration_sequence;
my $deletion_number;
my $check;
my $gene_size_altered;
my $left_bound_altered;
my $right_bound_altered;
my $gene_output;
my $gene_name;
my $gene_string;
my $gene_substring;
my $gene_substring_faa;
my @sequences;
my $gene_number;
my $file1;
my $newdir;
my $isolate_faa;
my $isolate;
my $reference;
my $seq_obj;
my $check_number;

foreach (@isolates) {
    $isolate = $_;

    #Makes a new directory for each isolate    

    $newdir = ($wkdir . $isolate . "/");
    
    system "echo mkdir $newdir";
    system "mkdir $newdir";

    $file1 = ($wkdir . $isolate . ".cns");
    $gene_number= 0;

    foreach (@genes) {
	    $gene = $_;
	    @sequences = ();
	    $deletion_number = 0;
	    $gene_size_altered = $gene_size[$gene_number];
	    $left_bound_altered = $left_bound[$gene_number];
	    $right_bound_altered = $right_bound[$gene_number];
	
	    open FILE1, $file1 or die "FILE $file1 NOT FOUND - $!\n";
	    foreach $_ (<FILE1>) {
		    $line = $_;

		    #Reads .pileup file and finds results for gene

		    if ($line =~ /$gene\s+(\d+)\s+(\w+)\s+(.*?)\s+(.*?):/) {
			    $location = $1;
			    $reference = $2;
			    $alteration = $3;
			    $check = $4;
			    if ($deletion_number > 0) {
	    			$consensus = "";
	    			$deletion_number--;
	    			$check_number=0;
    			} else {
	    			$consensus = $check;
	    			$check_number=1;
	    		}
				
			#Alters sequence for indels and removed flanks

			    if ($alteration =~ /\./) {
				    push @sequences, $consensus;
				} elsif ($alteration =~ /\+(\w+)/) {
				    $alteration_sequence = $1;
				    push @sequences, $reference;
				    push @sequences, $alteration_sequence;
				    if ($location < ($left_bound_altered + 1)) {
					    $left_bound_altered += (length $alteration_sequence);
				    } elsif ($location < ($gene_size_altered + 1 - $right_bound_altered)) {
					    $gene_size_altered += (length $alteration_sequence);
				    } elsif ($location < ($gene_size_altered + 1)) {
					    $right_bound_altered += (length $alteration_sequence);
				    }
				
			    } elsif ($alteration =~ /\-(\w+)/) {
				    $alteration_sequence = $1;
				    if($check_number < 1) {
				        push @sequences, $consensus;
				    } else {
				        push @sequences, $reference;
				        $deletion_number += (length $alteration_sequence);
				        if ($location < ($left_bound_altered + 1)) {
					        $left_bound_altered -= (length $alteration_sequence);
				        } elsif ($location < ($gene_size_altered + 1 - $right_bound_altered)) {
					        $gene_size_altered -= (length $alteration_sequence);
				        } elsif ($location < ($gene_size_altered + 1)) {
					        $right_bound_altered -= (length $alteration_sequence);
				        }
				    }
			    } elsif ($alteration =~ /\w/) {
				    push @sequences, $consensus;
			    }		
		    }		
	    }
	    close FILE1;
	
	    $gene_string = join("", @sequences);
	    $gene_substring = substr($gene_string, $left_bound_altered, ($gene_size_altered - $left_bound_altered - $right_bound_altered));
	
	    $gene_output = ($newdir . $gene . "__" . $isolate . ".fa");
	    $gene_name = ($gene);

	    #Prints out nucleotide sequence as a fasta file
	
	    open(OUT2, ">$gene_output") or die "Couldn't open OUT2 $gene_output $!\n";
	    print OUT2 ">";
	    print OUT2 $gene_name;
	    print OUT2 "__";
	    print OUT2 $isolate;
	    print OUT2 "\n";
	    print OUT2 $gene_substring;
	    close OUT2;
	    
	  
            #Converts nucleotide sequence to an amino acid sequence and prints out amino acid sequence as a fasta file

            $isolate_faa = ($newdir . $gene . "__" . $isolate . ".faa");

	    $seq_obj = Bio::Seq->new(-seq => $gene_substring, -alphabet =>'dna');
        
            $gene_substring_faa = $seq_obj->translate;
        
        
       
   	    open(OUT3, ">$isolate_faa") or die "Couldn't open OUT2 $isolate_faa $!\n";
	    print OUT3 ">";
	    print OUT3 $gene_name;
	    print OUT3 "__";
	    print OUT3 $isolate;
	    print OUT3 "\n";
	    print OUT3 $gene_substring_faa->seq;
	    close OUT3;

	    $gene_number++;
    }
}

