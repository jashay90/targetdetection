#!/bin/bash
# Julie Shay
# February 13, 2018
# short script to align targets in fewer lines...

INDIR=$1
n=$2
MEDIR="/mnt/nas/users/julie"
REF="$MEDIR/baitinginfo/targets.fa"
bwa mem $REF -t $n -a $INDIR/contam_1.fq.gz $INDIR/contam_2.fq.gz | samtools view -bT $REF | samtools sort > $INDIR/targets_new.bam
samtools index $INDIR/targets_new.bam
samtools flagstat $INDIR/targets_new.bam > $INDIR/targets_new_flagstat.txt
python $MEDIR/scripts/readcounts.py -b $INDIR/targets_new.bam -t $REF > $INDIR/target_new_counts.txt

