# Phenotype-Prediction-Pipeline
Alternative Splicing and Differential Expression Variable Selection for Predictive Phenotype Analysis Package Scripts

All scripts involved in the creation of the pipeline are included in this repository. There are also descriptions of the commands used to run various softwares which were steps along the way, such as Trimmomatic, TopHat, RSEM, iLOO, CuffLink, and EBSEQ. The organization is divided into three main folders with the various directions of research described: Tuxedo, RSEM, SeqGSEA. These folders are organized as follows:

* Tuxedo
	+ 1-TopHat
	+ 2-Cufflinks
	+ 3-BuzzCut
	+ 4-ElasticNet
	+ 5-SPLS

* RSEM
	+ 1-RSEM
	+ 2-EBSeq
	+ 3-BuzzCut
	+ 4-ElasticNet
	+ 5-SPLS

* SeqGSEA
	+ 1-GenerateMatrices
	+ 2-SeqGSEA
	+ 3-BuzzCut
	+ 4-ElasticNet
	+ 5-SPLS

There are also two other folders which contain scripts for other analyses that are used in the pipeline.

* Prepare-Fastq
	+ 1-RunTrimmomatic

* iLOO
	+ 1-iLOO
	+ 2-FilterGeneList
