#!/usr/bin/perl
#
#	Copy (C) 2016 - 2017 Massey University.
#	Written by Patrick Biggs PhD, and modified by Samuel Bloomfield
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY. See the GNU General Public License for
#	more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  This software is provided 'as-is', without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.
#
#
#############################
#
#	This is a Perl script called 'Assembly_scipt.pl', and is provided as
#	part of the manuscript for Bloomfield et al. via the Github repository.
#	It uses the assembler Velvet to assemble random subsets of reads to
#	maximise the kmer and kmer coverage space to try to find the optimum
#	assembly. In this versin, three isolates are analysed from two
#	different directory of Illumina fastq reads.
#
#	It should be noted that this is a highly individualised script for this
#	work only, and whilst the concepts and workflow described here can be
#	applied again, it will require a significant amount of customisation by
#	the user for their local setup.
#
#	Required software:
#		Scripts:

#			getSubset.py - by Mauro Truglio.  Also included in the
#				repository with this script.

#			fastq2fasta.pl - by Brian J Knaus.  Also included in the
#				repository with this script.
#		Software:
#			Velvet - the most recent version can be found at:
#				https://www.ebi.ac.uk/~zerbino/velvet/
#			SolexaQA++ - the most recent version can be found at:
#				http://solexaqa.sourceforge.net
#
#	Description of input options:
#		mode:
#			Chooses which one of the three main modes is to be run.
#			have to be run in the order: 'assemble' then 'analyse'
#			then 'retrieve'
#		kmerStart:
#			The lower end of a kmer range over which assembly metrics can
#			be outputted for inspection
#		kmerEnd:
#			The upper end of a kmer range over which assembly metrics can
#			be outputted for inspection
#		NZGL:
#			A local folder used as is due to the local setup
#
#	Description of working modes:
#		assemble:
#			as defined in the 'tableCreate' and 'random' subroutines,
#			including the further nested 'velLoad' subroutine
#		analyse:
#			as defined in the 'summaryTable' subroutine
#		retrieve:
#			as defined in the 'retrieveContigs' subroutine
#
#####################################################################

use warnings;
use strict;
use Getopt::Long;
use DBI;
use Statistics::Descriptive;
use Bio::Seq;
use Bio::SeqIO;


##################
#	Declaration of input options.
##################

my ($mode, $kmerStart, $kmerEnd, $NZGL);

GetOptions ('mode:s'		=> \$mode,
			'kmerStart:s'	=> \$kmerStart,
			'kmerEnd:s'		=> \$kmerEnd,
			'NZGL:s'		=> \$NZGL);


##################
#	Location of base directory.
##################

my $base			= ("/pathway/to/base/directory/");
my $root			= ($base . $NZGL . "/");
my $scriptFolder	= ($base . "scripts/");
my $otherFolder		= ($root . "otherOutput/");


##################
#	Array of kmer lengths to form assemblies from
#	when using Velvet, and other required arrays.
##################

my @kmer		= (245, 235, 225, 215, 205, 195, 185, 175, 165, 155, 145, 135, 125, 115, 105, 95, 85, 75, 65, 55);
my @trimming	= ('processed');
my @runs		= ('A', 'B');
my @genomes;

##################
#	Common variables.
##################

my ($dbh, $sth, $datasource, $count, $rowcount, $statement, $joiner, $querystring);
my ($logTable, $logTableID, $contigTable, $contigTableID, $realFolder, $curRun, $curTrimmed, $curGenome, $f1, $f2, $i, $n, $pair1, $pair2, $newFolder, $single, $shuffle, $prob, $sumDir, $curKmer, $curContigs, $outFile, $empty, $summary, $newGenome, $other);


##################
#	The location of other scripts required.
##################

my $pySubset	= ($scriptFolder. "getSubset.py");
my $fastq2fasta	= ($scriptFolder . "fastq2fasta.pl");
my $log			= ($root . $mode . "_LogBrief.txt");


##################
#	Record output to a process log.
##################

open (LOG, ">$log") or die ("couldn't open $log: $!\n");

print ("Process started at " . scalar(localtime) . ".\n");
print LOG ("Process started at " . scalar(localtime) . ".\n");


##################
#	Create the analysis tables.
##################

my $version		= ($NZGL . "_");

$logTable		= ($version . "Log");
$logTableID		= ($logTable . "ID");

$contigTable 	= ($version . "Contigs");
$contigTableID 	= ($contigTable . "ID");

$summary		= ($version . "Summary");

if (-e $otherFolder) {	print ("$otherFolder already exists.\n");
} else {				system "mkdir $otherFolder";
}


##################
#	Connect to the database.
##################

&dbConnect();


##################
#	Perform the work depending on the mode.
##################

if ($mode eq 'assemble') {

	##################
	#	Create core tables.
	##################

	&tableCreate();

	##################
	#	Start the work.
	##################

	foreach $curRun (@runs) {

		#	Name of run
		if ($curRun eq 'A') {

			#	Base directory containing reads
			$realFolder = ("Read_directory_1/");

			#	Array of isolates in the directory
			@genomes	= ('Isolate_1', 'Isolate_2', 'Isolate_3');

		} elsif ($curRun eq 'B') {
			$realFolder = ("Read_directory_2");
			@genomes	= ('Isolate_4', 'Isolate_5', 'Isolate_6');
		}

		print ("You are working with run $curRun...\n");
		print LOG ("You are working with run $curRun...\n");

		foreach $curGenome (@genomes) {
			foreach $curTrimmed (@trimming) {
				$newFolder	= ($root . $realFolder . $curTrimmed . "/");

				if ($curTrimmed eq 'processed') {

					#	Pathway to read fastq files in read directory
					$f1		= ($newFolder . "processed_" . $curGenome . "_L001_R1_001.fastq");
					$f2		= ($newFolder . "processed_" . $curGenome . "_L001_R2_001.fastq");

#					#	Probability cut-off for accepting nucleotides in reads
					# 	as required by SolexaQA
					$prob	= '0.01';
				}

				##################
				#	Define a set of subset sizes - input-dependent on the
				#	number of reads, and needs to be emperically set. We
				#	then perform assemblies with these subsets for each kmer.
				##################

				for ($i = 5; $i <= 8; $i++) {
					$n 	= (150000 * $i);


					##################
					#	Run the 'random' subroutine with these parameters.
					##################

					&random($n, $newFolder, $curGenome, $prob, $curTrimmed);
				}
			}
		}
	}


} elsif ($mode eq 'analyse') {

	##################
	#	Generate a summary table using this subroutine.
	##################

	&summaryTable();

} elsif ($mode eq 'retrieve') {

	##################
	#	Extract the chosen contigs using this subroutine.
	##################

	&retrieveContigs();

} else {
	print ("\nThere is something is wrong with the mode name, please try again.\n\n");

}

print ("Process finished at " . scalar(localtime) . ".\n");
print LOG ("Process finished at " . scalar(localtime) . ".\n");

close LOG;


##########################
##########################
##						##
##		subroutines		##
##						##
##########################
##########################

sub random {

	#############################################################
	#															#
	#	The main body of the code accepting variables to		#
	#	generate assembly contigs under a variety of			#
	#	parameters, in this case kmer and random reas subset.	#
	#															#
	#############################################################

	($n, $newFolder, $curGenome, $prob, $curTrimmed) = @_;

	my $subFolder	= ($newFolder . "subsets/");
	my $otherFolder	= ($newFolder . "assemblies/");
	my $nowFolder	= ($subFolder . $curGenome . "/");
	my $curSeq1		= ($nowFolder . $curGenome . "_1_" . $n . ".fastq");
	my $curSeq2		= ($nowFolder . $curGenome . "_2_" . $n . ".fastq");

	if (-e $subFolder) {	print ("$subFolder already exists.\n");
	} else {				system "mkdir $subFolder $otherFolder";
	}

	if (-e $nowFolder) {	print ("$nowFolder already exists.\n");
	} else {				system "mkdir $nowFolder";
	}

	system "python getSubset.py $f1 $f2 $n";
	system "mv $f1.subset $curSeq1";
	system "mv $f2.subset $curSeq2";

	print ("--and have made a random subset of $n at " . scalar(localtime) . ".\n");
	print LOG ("--and have made a random subset of $n at " . scalar(localtime) . ".\n");


	##################
	#	Perform DynamicTrim at the required probability
	#	threshold on processed reads only.
	##################

	my $DT1		= ($curSeq1 . ".trimmed");
	my $DT2		= ($curSeq2 . ".trimmed");
	my $moved	= ("./*fastq.trimmed");
	$pair1		= ($DT1 . ".paired1");
	$pair2		= ($DT1 . ".paired2");
	$single		= ($DT1 . ".single");
	my $shuffle	= ($subFolder . $curGenome . "/" . $n . "_shuffle.fastq");

	system "SolexaQA -p $prob $curSeq1 $curSeq2 -d $nowFolder";
	system "DynamicTrim -p $prob $curSeq1 $curSeq2";
	system "mv $moved $nowFolder";
	system "LengthSort -l 50 $DT1 $DT2";
	system "shuffleSequences_fastq.pl $pair1 $pair2 $shuffle";

	print ("---and have run DynamicTrim and LengthSort at " . scalar(localtime) . ".\n");
	print LOG ("---and have run DynamicTrim and LengthSort at " . scalar(localtime) . ".\n");


	##################
	# 	Perform the assemblies.
	##################

	print ("-On to the assembly at " . scalar(localtime) . ".\n");
	print LOG ("-On to the assembly at " . scalar(localtime) . ".\n");

	$sumDir	= ($newFolder . "assemblies/" . $curGenome . "/");

	if (-e $sumDir) {	print ("$sumDir already exists.\n");
	} else {			system "mkdir $sumDir";
	}

	foreach $curKmer (@kmer) {
		my $resultsDir		= ($sumDir . "kmer_" . $curKmer . "_" . $n);
		my $contigNewName 	= ($n . "_contigs". $curKmer . ".fa");
		my $statsNewName 	= ($n . "_stats". $curKmer . ".txt");
		my $logNewName		= ($n . "_log". $curKmer . ".txt");

		system "mkdir $resultsDir";
		system "velveth $resultsDir $curKmer -fastq -shortPaired $shuffle -short $single";
		system "velvetg $resultsDir -exp_cov auto -cov_cutoff auto -ins_length 500 -min_contig_lgth 200 -scaffolding no";

		system "cp $resultsDir/contigs.fa $sumDir/$contigNewName";
		system "cp $resultsDir/stats.txt $sumDir/$statsNewName";
		system "cp $resultsDir/Log $sumDir/$logNewName";
		system "rm -r $resultsDir";

		##################
		# 	Process the Velvet log and contigs with the subroutine.
		##################

		&velLoad($curKmer, $sumDir, $curRun, $n, $curGenome, $curTrimmed);

		print LOG ("\tAssembly on a kmer of $curKmer complete at " . scalar(localtime) . ".\n");
		print ("\tAssembly on kmer of $curKmer complete at " . scalar(localtime) . ".\n");
	}

	return ($n, $newFolder, $curGenome, $prob, $curTrimmed);

}


sub velLoad {

	#####################################################
	#													#
	#	This subroutine converts the Velvet assembly	#
	#	contigs into a MySQL parsable form and loads	#
	#	them into MySQL tables							#
	#													#
	#####################################################

	($curKmer, $sumDir, $curRun, $n, $curGenome, $curTrimmed) = @_;

	##################
	#	Process the Velvet log files for MySQL loading
	##################

	my $log 	= ($sumDir . '/' . $n . '_log' . $curKmer . '.txt');
	my $logMod	= ($sumDir . '/' . $n . '_log' . $curKmer . 'Mod.txt');

	open (IN1, "<$log") 		or die ("couldn't open $log: $!\n");
	open (OUT2, ">$logMod") 	or die ("couldn't open $logMod: $!\n");

	while (<IN1>) {
		chomp;
		if ($_ =~ m/Final graph has (\d+) nodes and n50 of (\d+), max (\d+), total (\d+), using (\d+)\/(\d+) reads/) {
			my $nodes 			= $1;
			my $N50 			= $2;
			my $max 			= $3;
			my $totalLength 	= $4;
			my $usedReads 		= $5;
			my $totalReads 		= $6;
			my @logArray = ($curTrimmed, $curGenome, $n, $curKmer, $nodes, $N50, $max, $totalLength, $usedReads, $totalReads);
		print OUT2 (join ("\t", @logArray), "\n");
		}
	}

	close IN1;
	close OUT2;

	$sth = $dbh->prepare (qq{load data local infile '$logMod' into table $logTable});	$sth->execute();

	##################
	#	Process the Velvet contigs for MySQL loading.
	##################

	$curContigs = ($sumDir . '/' . $n . '_contigs' . $curKmer . '.fa');
	$outFile 	= ($sumDir . '/' . $n . '_contigs' . $curKmer . '.txt');
	$empty	= ("");

	open (OUT, ">$outFile") or die ("couldn't open $outFile: $!\n");

	my $in = Bio::SeqIO->new(	-file 	=> "$curContigs",	-format => 'Fasta');

	while (my $seq = $in->next_seq()) {
		my $thisSeq 	= $seq->seq();
		my $thisHeader 	= $seq->primary_id();
		if ($thisHeader =~ m/NODE_(\d+)_length_(\d+)_cov_(\d+.\d+)/) {
			my $nodeID			= $1;
			my $nodeLength 		= $2;
			my $nodeCoverage 	= $3;
			my $realLength 		= $nodeLength + $curKmer - 1;
			my @results = ($curTrimmed, $curGenome, $n, $curKmer, $nodeID, $nodeLength, $nodeCoverage, $thisHeader, $realLength, $thisSeq, $empty);
		print OUT (join ("\t", @results), "\n");
		}
	}

	close OUT;

	print ("Contigs made from a kmer of $curKmer processed at " . scalar(localtime) . ".\n");

	system "dos2unix $outFile";

	$sth = $dbh->prepare (qq{load data local infile '$outFile' into table $contigTable});	$sth->execute();

	system "rm $logMod $outFile";

	return ($curKmer, $sumDir, $curRun, $n, $curGenome, $curTrimmed);
}


sub tableCreate {

	#########################################
	#										#
	#	A small subroutine to set up some	#
	#	commonly needed MySQL tables		#
	#										#
	#########################################

	##################
	#	Create a contigs table.
	##################

	$sth = $dbh->prepare (qq{drop table if exists $contigTable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $contigTable (processing varchar(20), genome varchar(20), subset mediumint, kmer smallint, nodeID mediumint, nodeLength mediumint, nodeCoverage decimal(12,6), subjectName varchar(40), realLength mediumint, sequence mediumtext, empty char(1))});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $contigTable add index index1(processing, genome, subset, kmer)});	$sth->execute();

	##################
	#	Create a log table.
	##################

	$sth = $dbh->prepare (qq{drop table if exists $logTable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $logTable (processing varchar(20), genome varchar(20), subset mediumint, kmer smallint, nodes mediumint, N50 int, max int, totalLength int, usedReads int, totalReads int)});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $logTable add index index1(processing, genome, subset, kmer)});	$sth->execute();
}


sub summaryTable {

	#####################################################
	#													#
	#	A major subroutine to geenrate summary tables	#
	#	both for the overall dataset and individual		#
	#	isolates.  It also generates a ranking summary	#
	#	allowing the 'best' acontigs from an assembly	#
	#	subset to be chosen for subsequent analysis.	#
	#													#
	#####################################################

	##################
	#	Update the contig table.
	##################

	$sth = $dbh->prepare (qq{alter table $logTable add column $logTableID mediumint primary key auto_increment first});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $contigTable add column $contigTableID mediumint primary key auto_increment first});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $contigTable add column $logTableID mediumint after $contigTableID});	$sth->execute();
	$sth = $dbh->prepare (qq{update $contigTable c, $logTable l set c.$logTableID = l.$logTableID where l.processing = c.processing and l.genome = c.genome and l.subset = c.subset and l.kmer = c.kmer});	$sth->execute();


	##################
	#	Create and work on a summary table that can summarise
	#	contigs across the whole assembly, as well as those
	#	greater than 500bp and 1000bp.
	##################

	$sth = $dbh->prepare (qq{drop table if exists $summary});	$sth->execute();

	$sth = $dbh->prepare (qq{create table $summary select v.*, a.contigsAll, a.lengthAll, b.contigsGE500, b.lengthGE500, c.contigsGE1000, c.lengthGE1000 from
	(SELECT *, 100*(usedReads/totalReads) as pCreads FROM $logTable) v left join
	(SELECT genome, processing, subset, kmer, count(*) as contigsAll, sum(realLength) as lengthAll FROM $contigTable v where kmer between '$kmerStart' and '$kmerEnd' group by subset, processing, kmer, genome) a
	on v.kmer = a.kmer and v.subset = a.subset and v.processing = a.processing and v.genome = a.genome left join
	(SELECT genome, processing, subset, kmer, count(*) as contigsGE500, sum(realLength) as lengthGE500 FROM $contigTable v where kmer between '$kmerStart' and '$kmerEnd' and realLength >= 500 group by subset, processing, kmer, genome) b
	on v.kmer = b.kmer and v.subset = b.subset and v.processing = b.processing and v.genome = b.genome left join
	(SELECT genome, processing, subset, kmer, count(*) as contigsGE1000, sum(realLength) as lengthGE1000 FROM $contigTable v where kmer between '$kmerStart' and '$kmerEnd' and realLength >= 1000 group by subset, processing, kmer, genome) c
	on v.kmer = c.kmer and v.subset = c.subset and v.processing = c.processing and v.genome = c.genome});	$sth->execute();

	my $outSum	= ($otherFolder . "All_VelContig_kS"  . $kmerStart . "_kE" . $kmerEnd . ".txt");

	open (OUT, ">$outSum") or die ("couldn't open $outSum: $!\n");

	my @header = ('processing', 'genome', 'subset', 'kmer', 'nodes', 'N50', 'max contig', 'total length', 'reads used', 'number reads', '% reads', 'all contigs', 'all length', 'GE500 contigs', 'GE500 length', 'GE1000 contigs', 'GE1000 length');

	print OUT (join("\t", @header), "\n");

	$statement = ("select * from $summary");

	&statementPull ($statement, "\t");

	close OUT;

	print ("Summary table made at " . scalar(localtime) . ".\n");
	print LOG ("Summary table made at " . scalar(localtime) . ".\n");

	##################
	#	Create a summary table for each isolate.
	##################

	my $outGenomes	= ($otherFolder . "genomes4ThisRun.txt");

	open (OUT, ">$outGenomes") or die ("couldn't open $outGenomes: $!\n");

	$statement = ("select genome, processing from $summary group by genome");

	&statementPull ($statement, "\t");

	close OUT;

	print ("Genome list made at " . scalar(localtime) . ".\n");
	print LOG ("Genome list made at " . scalar(localtime) . ".\n");


	##################
	#	Analyse the isolate by using an equal weighted ranking process
	#	for the assemblies for a given isolate usng a combination of:
	#	fewest contigs; longest overall length; longest contig length;
	#	and largest N50 value.
	##################

	open (IN, "<$outGenomes") or die ("couldn't open $outGenomes: $!\n");

	while (<IN>) {
		chomp;
		($newGenome, $other)	= split;

		my $rankingT	= ($version . $newGenome . "sumRankingNew");
		my $rankOut		= ($otherFolder . $rankingT . ".txt");

		$sth = $dbh->prepare (qq{drop table if exists $rankingT});	$sth->execute();
		$sth = $dbh->prepare (qq{create table $rankingT select $logTableID, processing, genome, subset, kmer, N50, max, contigsAll, lengthAll from $summary where genome = '$newGenome'});	$sth->execute();

		my @stats	= ('N50', 'max', 'contigsAll', 'lengthAll');

		foreach my $curStats (@stats) {
			my $uniqT	= ("uniqData4_" . $curStats);
			my $uniqTID	= ($uniqT . "ID");
			my $newCol	= ("counts_" . $curStats);

			$sth = $dbh->prepare (qq{drop table if exists $uniqT});	$sth->execute();

			if ($curStats eq 'contigsAll') {
				$sth = $dbh->prepare (qq{create table $uniqT select $curStats, count(*) as $newCol from $summary where genome = '$newGenome' group by $curStats});	$sth->execute();
			} else {
				$sth = $dbh->prepare (qq{create table $uniqT select a.*, count(*) as $newCol from (select $curStats from $summary where genome = '$newGenome' order by $curStats desc) a group by a.$curStats order by $curStats desc});	$sth->execute();
			}

			$sth = $dbh->prepare (qq{alter table $uniqT add column $uniqTID mediumint primary key auto_increment});	$sth->execute();
			$sth = $dbh->prepare (qq{alter table $rankingT add column $uniqTID mediumint});	$sth->execute();
			$sth = $dbh->prepare (qq{update $rankingT r, $uniqT t set r.$uniqTID = t.$uniqTID where r.$curStats = t.$curStats});	$sth->execute();

			$sth = $dbh->prepare (qq{drop table if exists $uniqT});	$sth->execute();
		}

		$sth = $dbh->prepare (qq{alter table $rankingT add column overallRank mediumint});	$sth->execute();
		$sth = $dbh->prepare (qq{update $rankingT set overallRank = (uniqData4_N50ID + uniqData4_maxID + uniqData4_contigsAllID + uniqData4_lengthAllID)});	$sth->execute();

		open (OUT, ">$rankOut") or die ("couldn't open $rankOut: $!\n");

		my @header1 = ($logTableID, 'processing', 'genome', 'subset', 'kmer', 'N50', 'max', 'contigsAll', 'lengthAll',  'N50Rank', 'maxRank', 'contigsAllRank', 'lengthAllRank', 'overallRank');

		print OUT (join("\t", @header1), "\n");

		$statement = ("select * from $rankingT order by overallRank");

		&statementPull ($statement, "\t");

		close OUT;

		$sth = $dbh->prepare (qq{drop table if exists $rankingT});	$sth->execute();

		print ("--Analysis of $newGenome complete at " . scalar(localtime) . ".\n");
		print LOG ("--Analysis of $newGenome complete at " . scalar(localtime) . ".\n");

	}

	close IN;
}


sub retrieveContigs {

	#############################################################
	#															#
	#	Mode to extract identified contigs from tables after	#
	#	after manual inspection and list generation.  Contigs	#
	#	are currently hardcoded after inspection.				#
	#															#
	#############################################################

	my @complete	 = ('Isolate_1', 'Isolate_2', 'Isolate_3', 'Isolate_4', 'Isolate_5', 'Isolate_6');
	my @contigs;

	foreach my $curComplete (@complete) {

		print ("You are working with sample $curComplete at " . scalar(localtime) . ".\n");
		print LOG ("You are working with sample $curComplete at " . scalar(localtime) . ".\n");

		##################
		#	Assemblies to extract for each isolate with example contig numbers,
		#	here with two contigs per initial isolate
		##################

		if ($curComplete eq 'Isolate_1') 		{	@contigs	= (480, 440);
		} elsif ($curComplete eq 'Isolate_2') {		@contigs	= (520, 540);
		} elsif ($curComplete eq 'Isolate_3') {		@contigs	= (559, 572);
		} elsif ($curComplete eq 'Isolate_4') {		@contigs	= (575, 592);
		} elsif ($curComplete eq 'Isolate_5') {		@contigs	= (601, 615);
		} elsif ($curComplete eq 'Isolate_6') {		@contigs	= (632, 640);
		}

		foreach my $curContig (@contigs) {
			my $seqFile	 = ($otherFolder . $curComplete . "_" . $curContig . ".fa");

			open (OUT, ">$seqFile") or die ("couldn't open $seqFile: $!\n");

			$statement = ("select concat('>', '$curComplete', '_', $contigTableID), sequence from $contigTable where $logTableID = '$curContig' order by realLength desc");

			&statementPull ($statement, "\n");

			close OUT;

			print ("\t and have pulled the contigs from contigID $curContigs at " . scalar(localtime) . ".\n");
			print LOG ("\t and have pulled the contigs from contigID $curContigs at " . scalar(localtime) . ".\n");
		}
	}
}


sub statementPull {

	#############################################################
	#															#
	#	DBI-based subroutine to extract data from MySQL tables.	#
	#															#
	#############################################################

	($statement, $joiner) = @_;

	$sth = $dbh->prepare (qq{$statement});	$sth->execute();

	$count++;

	while (my @row_items = $sth->fetchrow_array ()) {
		$rowcount++;
		print OUT (join ("$joiner", @row_items), "\n");
		} unless ($rowcount) {
		print OUT ("No data to display\n");
	}

	return ($statement, $joiner);
}


sub dbConnect {

	#################################################################
	#																#
	#	Subroutine for database connection via DBI.					#
	#	Adjust with your username and password for your connection	#
	#																#
	#################################################################

	$count = 0;
	$rowcount = 0;

	# 	Name of database
	$datasource = "Database_name";

	# 	Username and password for database
	$dbh = DBI->connect($datasource, 'username', 'password');
	$querystring = '';
}
