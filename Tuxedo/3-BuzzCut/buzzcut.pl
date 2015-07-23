#!usr/bin/perl

use warnings;

use Getopt::Long;

my $cuff_format = 0;

my $p_value = 0.05;
my $fold = 0.6;
my $filterP = 0;
my $filterFold = 0;
my $help = 0;

GetOptions ( 	
				"p=f" => \$p_value,
				"f=f" => \$fold,
				"filterP" => \$filterP,
				"filterFold" => \$filterFold,
				"help" => \$help,
				"d|diff_file=s" => \$diff_file,
				"c|count_file=s" => \$count_file,
				"o|output=s" => \$output
			);


if($help){

	&printHelp;
	exit;

} 

if(!$diff_file or !$count_file or !$output){

	&printHelp;
	die "Please run buzzcut with the correct arguments\n";

}

#my($diff_file, $count_file, $output) = @ARGV;
open $diff_file, '<', $diff_file or die "Can't open $diff_file: $!\n";
open $count_file, '<', $count_file or die "Can't open $count_file: $!\n";
open $output, '>', $output or die "Can't open $output: $!\n";


my $num_sig = 0;
my @sig_genes;

while(<$diff_file>){
	
	chomp;
	my @line = split(/\t/);

	#Detect File Format...if test_id is first header then CuffDiff
	if($line[0] =~ /test_id/){ 
		$cuff_format = 1; 
		next;
	} elsif( $line[0] =~ /"PPEE"/){
		next;
	}

	#Filter on P-Value
	if($filterP and !$filterFold){

		#If File is in CuffDiff Format
		if($cuff_format){

			#Store Gene_ID in array if P-value is significant
			if($line[11] <= $p_value) { $sig_genes[scalar @sig_genes] = "$line[0]"; }

		} else {

			#Store Gene_ID in array if PPEE value is significant in EBSeq
			if($line[1] <= $p_value) { $sig_genes[scalar @sig_genes] = "$line[0]"; }

		}

	} elsif(!$filterP and $filterFold){

		#If File is in CuffDiff Format
		if($cuff_format){

			#Store Gene_ID in array if Fold Change is significant
			if(abs($line[9]) >= $fold) { $sig_genes[scalar @sig_genes] = "$line[0]"; }

		} else {

			#Store Gene ID in array if Fold Change is significant in EBSeq
			if(abs(log($line[4])/log(2)) >= $fold) { $sig_genes[scalar @sig_genes] = "$line[0]"; }

		}
	
	} elsif($filterP and $filterFold){

		#If File is in CuffDiff Format
		if($cuff_format){

			#Store Gene_ID in array if both P-value and Fold Change are significant
			if($line[11] <= $p_value and abs($line[9]) >= $fold) { $sig_genes[scalar @sig_genes] = "$line[0]"; }

		} else {

			#Store Gene_ID in array if both P-value and Fold Change are significant in EBSeq
			if($line[1] <= $p_value and abs(log($line[4])/log(2)) >= $fold){ $sig_genes[scalar @sig_genes] = "$line[0]"; }

		}

	} else {

		die "Error: At least one criteria must be selected to filter\n";

	}

}

$num_sig = scalar @sig_genes;
print "$num_sig Genes/Isoforms Selected\n";
print "Begin Parsing Data\n";

my $j = 0;
my $foundFirstMatch = 0;
my $finishFirstMatch = 0;
my @firstMatchCounts;
my @firstMatchSamples;
my $num_samples = 0;
my $sample_index = 0;

while(<$count_file>){

	chomp;
	@line = split(/\t/);

	#If analyzing CuffDiff Format
	if($cuff_format){

		#If on the first significant gene
		if($j == 0){

			#Skip the first line in count file
			if($line[0] =~ /tracking_id/){
			
				next;
			
			} else {
				
				#Iterate through file until finding match of first significant gene
				if(!$foundFirstMatch and $line[0] ne $sig_genes[$j]){

					next;
				
				} elsif($foundFirstMatch and $line[0] ne $sig_genes[$j]){
					
					$finishFirstMatch = 1;

#					#Print line to output file
#					print $output "tracking_ID\t";
#
#					for(my $i = 0; $i < scalar @firstMatchSamples; $i++){
#						print $output "$firstMatchSamples[$i]";
#						if($i < scalar @firstMatchSamples - 1){
#							print $output "\t";
#						} else {
#							print $output "\n";
#						}
#					}

					print $output "$sig_genes[$j]\t";

					for(my $i = 0; $i < scalar @firstMatchCounts; $i++){
						print $output "$firstMatchCounts[$i]";
						if($i < scalar @firstMatchCounts - 1){
							print $output "\t";
						} else {
							print $output "\n";
						}
					}

					$j++;

					#Move on to regular parsing
				
				} else {
					
					$foundFirstMatch = 1;

					$firstMatchCounts[$num_samples] = $line[6];
					$firstMatchSamples[$num_samples] = "$line[1]$line[2]";

					$num_samples++;

				}

			}

		} else {

			if($j == scalar @sig_genes){

				close($output);
				print "All significant genes/isoforms parsed to $output!\n";
				exit;

			} elsif($line[0] eq $sig_genes[$j] and $j < scalar @sig_genes){

				if($sample_index == $num_samples - 1){

					print $output "$line[6]\n";

					$sample_index = 0;

					$j++;

				} else {

					if($sample_index == 0){

						print $output "$line[0]\t";
					}

					print $output "$line[6]\t";

					$sample_index++;

				}

			}

		}

	#If analyzing EBSeq Format
	} else {

		#Form Sample List from Header Line
		if($line[0] =~ /\A\Z/){

			#Print line to output file
#			print $output "Gene_ID\t";
#
#			for(my $i = 0; $i < scalar @line; $i++){
#				print $output "$line[$i]";
#				if($i < scalar @line - 1){
#					print $output "\t";
#				} else {
#					print $output "\n";
#				}
#			}

		} else {

			#Match Significant Gene Names in the files
			if($j == scalar @sig_genes){

				close($output);
				print "All significant genes/isoforms parsed to $output!\n";
				exit;

			} elsif($line[0] eq $sig_genes[$j] and $j < scalar @sig_genes){

				for(my $z = 0; $z < scalar @line; $z++){
					print $output "$line[$z]";
					if($z < scalar @line - 1){
						print $output "\t";
					} else {
						print $output "\n";
					}
				}

				$j++;

			} elsif($line[0] gt $sig_genes[$j]){

				$j++;

			}

		}

	}

}

sub printHelp {

	print "\n";
	print "Usage: buzzcut.pl <Options> -d <Diff File> -c <Counts File> -o <Output File>\n";
	print "\n";
	print "Required Inputs:\n";
	print "-d\tDifferential Expression Data (CuffDiff: .diff or EBSeq results)\n";
	print "-c\tTranscript Count Data (CuffDiff: .read_group_tracking or EBSeq matrix)\n";
	print "-o\tOutput Path\n";
	print "\n";
	print "Options:\n";
	print "-p\tSpecify <= p-value (Default: 0.05)\n";
	print "-f\tSpecify >= fold change (Default: 0.6)\n";
	print "-filterP\tFilter significant transcripts by p-value\n";
	print "-filterFold\tFilter significant transcripts by fold change\n";
	print "\n";

}

