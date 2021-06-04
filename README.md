# AMRPlusPlus and bait-capture target detection
These are the scripts used to run AMRPlusPlus and bait-capture target detection outside of Galaxy. The method for detection follows the same basic pipeline as [AMRPlusPlus](https://github.com/meglab-metagenomics/amrplusplus_v2). Reads are trimmed with Trimmomatic, then there is an optional host-filtering step with BWA MEM, and then reads are aligned to MegaRes genes and bait-capture targets with BWA MEM, and custom scripts are used to cluster and summarize results.

## Getting set up
OLC's bait-capture target sequences and clusters are included in this repository. To run AMRPlusPlus, you need to download MEGARes v1.0.1 [here](https://megares.meglab.org/download/index.php).

You can build a conda environment for these tools using conda-spec-file.txt:

	conda create --name myenv --file conda-spec-file.txt

To run AMRPlusPlus, you also need to install the resistome tool. You can find it [here](https://github.com/cdeanj/resistomeanalyzer).

### amrplusplus.sh
amrplusplus.sh runs the AMRPlusPlus pipeline.

	Usage: amrplusplus.sh -1 in1.fq -2 in2.fq -o outfolder
	use -t if reads are already trimmed
	use -c to provide a contaminant genome (either fasta or a bwa index), or 
	-b to provide a bam file if you've already aligned reads against the contaminant genomes
	use -n to specify number of CPUs. default is 8
	use -a to specify the path to adapters if you are trimming reads
	use -d to specify the path to the megares database
	use -x if you want the pipeline to stop before aligning reads to MEGARes

### aligntargets.sh
Once amrplusplus.sh is finished running (including host contamination filtering), you can use aligntargets.sh to align reads to bait-capture targets.

	Usage: aligntargets.sh infolder cpus

### readcounts.py
After running aligntargets.sh, you can use this script to generate read counts from the bam output file.

	usage: readcounts.py [-h] [-t TARGETS] [-b BAM]
	
	optional arguments:
	  -h, --help  show this help message and exit
	  -t TARGETS  fasta file to which reads were aligned
	  -b BAM      bam alignment file

### clusterandmerge_multi.py
Given a folder containing multiple subfolders where AMRPlusPlus and bait-capture target detection has been done, clusterandmerge_multi.py will combine results from multiple samples into a single table. This script will also cluster targets according to one or more cluster specification files (two specification files, allthethings_amr.tsv and allthethings_plasmid.tsv, are included in this repository, but any tab-delimited files with "ID" and "Cluster" columns could be used).

	usage: clusterandmerge_multi.py [-h] -i INFOLDER [INFOLDER ...]
                                [-s SUFFIX [SUFFIX ...]] [-n FNAME]
                                [-flagstat FLAGSTAT]
                                [-c CLUSTERS [CLUSTERS ...]] [-p MINPROP]
                                [-o OUTPREFIX] [-drop DROP]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -i INFOLDER [INFOLDER ...]
	                        Input Folder(s)
	  -s SUFFIX [SUFFIX ...]
	                        Suffix for samples coming from each folder
	  -n FNAME              Name of target alignment summary file within each
	                        subfolder
	  -flagstat FLAGSTAT    Optional flagstat file name within each subfolder. If
	                        you set this, the script will output a column for
	                        proportion of reads within each sample.
	  -c CLUSTERS [CLUSTERS ...]
	                        Tab-delimited files with genes and clusters
	  -p MINPROP            minimum proportion of gene covered
	  -o OUTPREFIX          output file prefix
	  -drop DROP            drop Gene Fraction column from output
