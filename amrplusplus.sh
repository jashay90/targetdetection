#!/bin/bash
# Julie Shay
# February 13, 2018
# Changing my original amrplusplus script to match Rahat's pipeline instead of the original
# amrplusplus pipeline

DDIR="/mnt/nas/users/julie/megares_v1.01"
adapters="/mnt/nas/users/julie/envs/amrplusplus/share/trimmomatic/adapters/TruSeq3-PE.fa"
n=8

while getopts '1:2:o:c:tb:n:ha:d:' flag; do
	case $flag in
		1)
			IN1=$OPTARG
			;;
		2)
			IN2=$OPTARG
			;;
		o)
			OUT=$OPTARG
			;;
		c)
			CONTAM=$OPTARG
			;;
		t)
			TRIMMED=1
			;;
		b)
			CONTAMBAM=$OPTARG
			;;
		n)
			n=$OPTARG
			;;
		h)
			h=1
			;;
		a)
			adapters=$OPTARG
			;;
		d)
			DDIR=$OPTARG
	esac
done

if [ -n "$h" ]; then
	echo "Usage: $0 -1 in1.fq -2 in2.fq -o outfolder"
	echo "use -t if reads are already trimmed"
	echo "use -c to provide a contaminant genome (either fasta or a bwa index), or -b to provide a bam file if you've already aligned reads against the contaminant genomes"
	echo "use -n to specify number of CPUs. default is 8"
	echo "use -a to specify the path to adapters if you are trimming reads"
	echo "use -d to specify the path to the megares database"
	exit
fi
if [ -z $(which bwa) ]; then
	echo "You need to have bwa installed to use this script."
	exit
fi
if [ -z $(which samtools) ]; then
	echo "You need to have samtools installed to use this script."
	exit
fi
if [ -z $(which resistome) ]; then
	echo "You need to have resistome installed to use this script."
	echo "Look for it at https://github.com/cdeanj/resistomeanalyzer"
	exit
fi
if [ -z "$OUT" ]; then
	echo "You need to provide an output folder name. Use -h to get help message."
	exit
fi

DPREFIX="${DDIR}/megares"
DB="${DPREFIX}_database_v1.01.fasta"
AN="${DPREFIX}_annotations_v1.01.csv"
if [ ! -f $DB ]; then
	echo "Are you sure you input the correct path to the megares database?"
	exit
fi
if [ ! -f "${DPREFIX}.amb" ]; then
	mkdir -p $OUT
	DPREFIX="$OUT/megares"
	bwa index -p $DPREFIX $DB
fi



difchar=$(cmp -l <(echo -n $IN1) <(echo -n $IN2) | wc -l)
if [ ${#difchar} -ne 1 ]; then
	echo "did you really input paired fastq files?"
	exit
fi

mkdir -p $OUT

FQ1="$OUT/contam_1.fq"
FQ2="$OUT/contam_2.fq"
if [ -z "$CONTAMBAM" ]; then
	difchar=$(cmp -l <(echo -n $IN1) <(echo -n $IN2) | wc -l)
	if [ ${#difchar} -ne 1 ]; then
		echo "did you really input paired fastq files?"
		exit
	fi
	if [ -z "$TRIMMED" ]; then
		if [ ! -f $adapters ]; then
			echo "Please specify an adapter file with the -a option to trim reads with this script."
			exit
		fi
		if [ -z $(which trimmomatic) ]; then
			echo "You need to have trimmomatic installed to trim reads with this script."
			echo "If your reads are already trimmed, use the -t option."
			exit
		fi
		mkdir $OUT/trimmomatic
		BASENAME=$OUT/trimmomatic/$(basename $IN1 | sed -e 's/\.gz$//' -e 's/\.fq$//' -e 's/\.fastq$//' -e 's/1$//' -e 's/R$//' -e 's/R1_00$//' -e 's/_$//')
		trimmomatic PE -threads $n -trimlog $OUT/trimmomatic/log.txt $IN1 $IN2 -baseout $BASENAME.fq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		gzip ${BASENAME}_*U.fq &
		IN1=${BASENAME}_1P.fq
		IN2=${BASENAME}_2P.fq
		
	fi
	if [ -n "$CONTAM" ]; then
		if [ -f ${CONTAM}.amb ]; then
			CINDEX=$CONTAM
		else
			CINDEX=$OUT/$(basename $CONTAM)
			bwa index -p $CINDEX $CONTAM
		fi
		CONTAMBAM=$OUT/bwa_contam.bam
		bwa mem -t $n -a $CINDEX $IN1 $IN2 | samtools view -b -f 0x0004 -f 0x0008 -f 0x0001 | samtools sort > $CONTAMBAM
	fi

fi

if [ -z "$CONTAMBAM" ]; then
	FQ1=$IN1
	FQ2=$IN2
else
# I added in the conversion to fastq because bwa sometimes uses inconsistent flags for aligned mates,
# which led to unmatched reads being included, and all of the reads appearing afterwards
# in the file being misaligned
	samtools fastq -1 $FQ1 -2 $FQ2 -s $OUT/shouldnotexist.fq $CONTAMBAM
	rm $OUT/shouldnotexist.fq &
	# print % reads passing contamination filter to an outfile
	pass=$(wc -l $FQ1 | awk '{print $1}')
	pass=$(expr $pass / 2)
	if [[ $IN1 == *.gz ]]; then
		total=$(gunzip -c $IN1 | wc -l | awk '{print $1}')
	else
		total=$(wc -l $IN1 | awk '{print $1}')
	fi
	total=$(expr $total / 2)
	prop=$(echo "$pass/$total" | bc -l)
	echo "Input Reads,Reads Passing Contamination Filter,Fraction Reads Passing Contamination Filter" > $OUT/contaminfo.txt
	echo ${total},${pass},$prop >> $OUT/contaminfo.txt
w
fi

bwa mem $DPREFIX -t $n -a $FQ1 $FQ2 | samtools view -bT $DB | samtools sort | samtools view > $OUT/bwa_amr.sam &
wait

if [ -n "$CONTAMBAM" ]; then
	gzip $FQ1 &
	gzip $FQ2 &
fi

if [ -z "$TRIMMED" ]; then
	gzip $IN1 &
	gzip $IN2 &
fi

resistome -ref_fp $DB -annot_fp $AN -sam_fp $OUT/bwa_amr.sam -gene_fp $OUT/resistome_genes -group_fp $OUT/resistome_groups -mech_fp $OUT/resistome_mechs -class_fp $OUT/resistome_classes -t 1 &
wait

samtools view -bT $DB $OUT/bwa_amr.sam > $OUT/bwa_amr.bam
samtools flagstat $OUT/bwa_amr.bam > $OUT/amr_flagstat.txt
rm $OUT/bwa_amr.sam
