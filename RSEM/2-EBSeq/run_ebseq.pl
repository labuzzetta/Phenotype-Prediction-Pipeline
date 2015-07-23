#!/usr/bin/perl

use strict;
use warnings;

my @DFSamples = (76, 77, 78, 79, 81, 93, 97, 109, 110, 111, 112);

my @RSamples = (80, 83, 84, 86, 87, 88, 89, 95, 98, 99);

for(my $i = 0; $i < scalar @DFSamples; $i++){

	for(my $j = 0; $j < scalar @RSamples; $j++){

		my $directory = "EBSeq_L20_Drop.".$DFSamples[$i]."DF.".$RSamples[$j]."R";
	
		chdir($directory);

		my $run_ebseq = "~/bin/rsem-1.2.21/rsem-run-ebseq "."EBSeq-TRAIN-Matrix_L2O-".$DFSamples[$i]."DF-".$RSamples[$j]."R.matrix.txt 10,9 EBSeq-TRAIN-results.genes.txt";

		system($run_ebseq);

		#system('cd ..');

		chdir("..");		

	}

}


