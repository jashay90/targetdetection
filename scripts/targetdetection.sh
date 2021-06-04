#!/bin/bash
# Julie Shay
# June 4, 2021
# I'm trying to make it a bit easier to run target detection without running all of AMRPlusPlus first

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
DB=$(echo $SCRIPT_DIR | sed 's/scripts/targets\/targets.fa/')
AMRPLUSPLUS="$SCRIPT_DIR/amrplusplus.sh"
RCOUNTS="$SCRIPT_DIR/readcounts.py"
n=8

amrplusplusoptions="-x"
while getopts '1:2:o:c:tb:n:ha:d:' flag; do
	case $flag in
		1)
			IN1=$OPTARG
			amrplusplusoptions="$amrplusplusoptions -1 $IN1"
			;;
		2)
			IN2=$OPTARG
			amrplusplusoptions="$amrplusplusoptions -2 $IN2"
			;;
		o)
			OUT=$OPTARG
			amrplusplusoptions="$amrplusplusoptions -o $OUT"
			;;
		c)
			amrplusplusoptions="$amrplusplusoptions -c $OPTARG"
			;;
		t)
			amrplusplusoptions="$amrplusplusoptions -t"
			;;
		b)
			amrplusplusoptions="$amrplusplusoptions -b $OPTARG"
			;;
		n)
			n=$OPTARG
			amrplusplusoptions="$amrplusplusoptions -n $n"
			;;
		h)
			h=1
			;;
		a)
			amrplusplusoptions="$amrplusplusoptions -a $OPTARG"
			;;
		d)
			DB=$OPTARG
			;;
	esac
done

if [ -n "$h" ]; then
	echo "Usage: $0 -1 in1.fq -2 in2.fq -o outfolder"
	echo "use -t if reads are already trimmed"
	echo "use -c to provide a contaminant genome (either fasta or a bwa index), or -b to provide a bam file if you've already aligned reads against the contaminant genomes"
	echo "use -n to specify number of CPUs. default is 8"
	echo "use -a to specify the path to adapters if you are trimming reads"
	echo "use -d to specify the path to the target fasta file"
	exit
fi

if [ ! -f $DB ]; then
	echo "Are you sure you input the correct path to the target sequences?"
	exit
fi

# run the trimming/decontam steps from amrplusplus.sh
$AMRPLUSPLUS $amrplusplusoptions

if [ ! -d $OUT ]; then
	echo "Something went wrong!"
	exit
fi

# use fastq files from amrplusplus if applicable
if [ -f "$OUT/contam_1.fq.gz" ]; then
	IN1="$OUT/contam_1.fq.gz"
	IN2="$OUT/contam_2.fq.gz"
elif [ -f "$OUT/trimmomatic/$(basename $IN1 | sed -e 's/\.gz$//' -e 's/\.fq$//' -e 's/\.fastq$//' -e 's/1$//' -e 's/R$//' -e 's/R1_00$//' -e 's/_$//')_1P.fq.gz" ]; then
	BASENAME="$OUT/trimmomatic/$(basename $IN1 | sed -e 's/\.gz$//' -e 's/\.fq$//' -e 's/\.fastq$//' -e 's/1$//' -e 's/R$//' -e 's/R1_00$//' -e 's/_$//')"
	IN1="${BASENAME}_1P.fq.gz"
	IN2="${BASENAME}_2P.fq.gz"
fi


if [ ! -f "${DB}.amb" ]; then
	DPREFIX="$OUT/targets.fa"
	cp $DB $DPREFIX
	bwa index -p $DPREFIX $DB
else
	$DPREFIX=$DB
fi

bwa mem $DPREFIX -t $n -a $IN1 $IN2 | samtools view -bT $DPREFIX | samtools sort > $OUT/targets.bam
samtools index $OUT/targets.bam
samtools flagstat $OUT/targets.bam > $OUT/targets_flagstat.txt
python $RCOUNTS -b $OUT/targets.bam -t $DPREFIX > $OUT/target_counts.txt
