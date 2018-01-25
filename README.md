# Scripts for genomic analysis

The scripts contained in this repository were used for analysing, rearranging and compiling Genomic data. 


### Prerequisites

The scripts contained in this repository consist of Perl and Python scripts. One script also also utilises the Velvet assembler:

```
Perl v5.18 or higher
BioPerl v1.6 or higher
Python v2.7 or higher
Velvet v1.1 or higher
```

### Installing

The scripts in this depository can be downloaded directly from GitHub:


```
git clone https://github.com/samuelbloomfield/Scripts-for-genomic-analyses
```


## Running the scripts

Most of the scripts in this depository will have paths and arrays that will need to be re-written before they can be run:


### De novo assembly

The Assembly_scipt.pl script was written by Patrick Biggs and used to assemble genomes de novo using Miseq sequenced reads. 
The Assembly_scipt.pl script will need to be re-written to include the location of the root directory, the directory of the reads, the isolates to be assembled, the contig sizes to assemble, and the subsets of reads to use, and the assemblies to extract.
The script will return the QC data on all of the isolates reads, the QC results for all the assemblies, and the extracted assemblies.
The Assembly_script.pl script also relied on two other scripts: 'fastq2fasta.pl' by Brian J. Knaus and 'getSubset.py' by Mauro Truglio.
From the terminal the file can be run as:

```
perl Assembly_scipt.pl -NZGL Project_name -kmerStart Starting_kmer_value -kmerEnd Ending_kmer_value -mode assemble|analyse|retrieve
```


### kSNP read count

The Read_searcher_python.py script was used to count the number of reads that contain SNPs identified in kSNP.
The Read_searcher_python.py script requires a file containing all the kmer variants that kSNP identified. I have included 'kSNP_SNPs_example.txt' as an example of what the kmer variant file should look like.
The Read_searcher_python.py script also requires a tab-delaminated file with the location of the kmer variant file in the first column, the location of the first processed read files in the second column, the location of the second processed read files in the third column and the names of the output files for the isolate in the fourth column. I have included 'SNP_and_read_example.txt' as an examble of what this file should look like.
The script will return a file for each isolate that contains the number of reads that contain each SNP variant.
From the terminal the file can be run as:

```
python Read_searcher_python.py SNP_and_read_file.txt
```

### SNP codon and position

The SNP_location.pl script was used to determine the type of SNPs (e.g. synonymous/non-synonymous, intergenic/found within one or more genes, and the codon position within a gene that it is associated with) based on their location within a reference genome.
The SNP_location.pl requires two arguments: the location of the reference genome in .gff format, and the location of a text-delaminated file containing the SNP information.
The text-delaminated should contain: the position of the SNP on the reference file, whether the SNPs are forward (F) or reverse compliment (R), and the variations of the SNPs. An example of the text-delaminated file can be found in the example directory. Most SNP-identification programs that use a reference (e.g. Snippy) will give all the SNPs as forward, but for programmes that identify SNPs via kmers (e.g. kSNP), the SNPs may be reverse complement.
The script will return: the number of genes the SNP is located within, the codon position, and whether the SNP is synonymous or non-synonymous in the "*_SNP_types.txt" file; a summary of the number of synonymous and non-synonymous SNPs, and the number of SNPs that are intergenic, or associated with one or multiple genes in the "*_SNP_summary.txt" file; and for each gene that contains a SNP, the number of SNPs associated with that gene in the "*_gene_SNP_counts.txt" file.
From the terminal the file can be run as:

```
perl SNP_location.pl -GFF Reference.gff -SNP Text_delaminated_SNP_file.txt
```


### Clustalw check

The clustalw_alignment_output_check.pl script was to check if there were any sequences missing or if there were any nucleotides that did not perfectly align in the clustalw alignments.
The clustalw_alignment_output_check.pl script will need to be re-written to include the location of the working directory, location of the directory containing the clustalw alignments, number of isolates aligned, and the names of genes aligned.
The script will form an 'alignment_summary.txt' file with information on the sequences that align perfectly and those that do not. 
From the terminal the file can be run as:

```
perl clustalw_alignment_output_check.pl
```

### Roary variant extractor

The Roary_variant_extractor.pl script was used to extract the genetic variant sequences from the roary output. 
The Roary_variant_extractor.pl script will need to be re-written to include the location of the working directory, the name of the Roary output file, the flank length required and an array of the gene names as they are located in the Roary output file.
The working directory should contain a copy of all the .gff files used in roary along with the Roary output file.
The Roary output file should be a tab-demilinated file with the names of the isolates in the first columns and then the names of the genes in the following columns. If roary did not find the gene then the blank space should be replaced with 'XXXXX'. This information can be extracted from the 'gene_presence_absence.csv' output file, but the rows and columns will need to be transposed. I have included 'Roary_output_example.txt' as an example of what the Roary output file should look like.
The script will return a copy of all gene variants in the "Gene_sequences_reduced" directory and the flank lengths, names and lengths of the gene variants in the "Gene_location_reduced_output.txt" file.
From the terminal the file can be run as:

```
perl Roary_variant_extractor.pl
```

### Flank remover

The Flank_remover.pl script was used convert the .cns file produced from SRST2 and Samtools into a nucleotide sequence, remove flanks and convert the sequence into an amino acid sequence.
The Flank_remover.pl script will need to be re-written to include the location of the working directory, gene names, the lengths of the genes, and the lengths of the left and right flanks.
The working directory should contain a copy of all the .cns files to be analysed.
The script will form a directory for each isolate and fill the directory with a nucleotide and amino acid copy of each gene sequence with the flanks removed. 
From the terminal the file can be run as:

```
perl Flank_remover.pl
```

### Protein differences to binary

The Protein_differences_to_binary.pl script was used to convert protein differences data into a binary table.
The Protein_differences_to_binary.pl script will need to be re-written to include the location of the working directory, the Protein differences file, and the names of the isolates.
The working directory should contain a copy of Protein differences file.
Each line in the Protein difference file should represent a separate difference and should contain a list of the isolates that contain the protein difference followed by a comma. I have included 'Protein_difference_example.txt' as an example of what the Protein differences file should look like.
The script will return a binary table regarding the presence of the protein differences in the "binary.txt" file.
From the terminal the file can be run as:

```
perl Protein_differences_to_binary.pl
```
