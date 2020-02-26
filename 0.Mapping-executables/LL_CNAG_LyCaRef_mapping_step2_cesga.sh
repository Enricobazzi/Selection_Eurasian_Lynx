#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

# The following softwares are used:
# BWA
# SAMtools
# Picard-tools
# GATK

module load gcc/6.4.0 bwa/0.7.17
module load gcc/6.4.0 samtools/1.8
module load gatk/3.7-0-gcfedb67
module load picard/2.18.14

THREADS=24
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa
SAMPLE=($(echo "${1}"))
OUT=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_bams



echo " - Merging all ${SAMPLE} BAMs and Re-Sorting -"

ls $OUT/${SAMPLE}_*_sorted_rg.bam  > $OUT/${SAMPLE}.bam.list

# check if sample has more than one bam, merge and resort if it does, change name if it doesn't :

if [[ $(wc -l < "${OUT}/${SAMPLE}.bam.list") -ge 2 ]]
 then

  samtools merge  -@ $THREADS -r $OUT/${SAMPLE}_merged.bam -b $OUT/"${SAMPLE}".bam.list

  BAMARRAY=($(cat ${OUT}/${SAMPLE}.bam.list))
  for k in ${BAMARRAY[@]}
   do
    echo " - Removing ${k} -"
    rm ${k}
  done

  echo " - Re-Sorting ${SAMPLE} -"
  samtools sort  -@ $THREADS $OUT/${SAMPLE}_merged.bam -o $OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam && rm $OUT/${SAMPLE}_merged.bam;

 else
  mv $(cat ${SAMPLE}.bam.list) ${SAMPLE}_LyCa_ref_sorted_rg.bam

fi

# Marking Duplicates
echo " - Marking Duplicates of $fastqid of ${SAMPLE} -"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
METRICS_FILE=$OUT/${SAMPLE}_rmdup.txt \
I=$OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam \
O=$OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg.bam

# Re-Sorting
echo " - Re-Sorting ${SAMPLE} -"
samtools sort $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam \
-@ $THREADS -o $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam
rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup.bam

# Indexing for GATK
echo " - Indexing ${SAMPLE} -"
samtools index $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam

# Realigning:
# RealignerTargetCreator
echo " - Realigner Target Creator on ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-nt $THREADS -R $REF -I $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam \
-o $OUT/${SAMPLE}_realignertargetcreator.intervals
# IndelRealigner
echo " - Realigning INDELS of ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
-R $REF -targetIntervals $OUT/${SAMPLE}_realignertargetcreator.intervals \
-I $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam \
-o $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
rm $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted.bam

# Get stats
echo " - Getting stats of ${SAMPLE} -"
samtools flagstat $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
> $OUT/${SAMPLE}_LyCa_ref_sorted_rg_rmdup_sorted_indelrealigner.stats
