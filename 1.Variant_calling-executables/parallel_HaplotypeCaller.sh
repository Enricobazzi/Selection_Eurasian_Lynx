#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH -n 1

module load gatk/4.1.4.1

###################################
## VARIABLE and PATHS definition ##
###################################

scaffold=($(echo ${1}))

i=($(echo ${2}))

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa

# Input bams - c_ll_ca_0249 and c_ll_ca_0253 didn't pass QC and are excluded:
INbamARRAY=($(ls /mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_bams/*_indelrealigner.bam | grep -vE "c_ll_ca_0249|c_ll_ca_0253"))


###############################
## GATK 4.1.1.0 CombineGVCFs ##
###############################

for k in {1..8}
 do
  echo "HaplotypeCaller of ${scaffold}_w${i}_w${k}.bed"
  OUTvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_vcfs/${scaffold}_w${i}_w${k}_LyCa_ref.vcf
  BED=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Scaffold_Partitions/${scaffold}_${i}/${scaffold}_w${i}_w${k}.bed
  gatk HaplotypeCaller \
     -R $REF \
     $(for bam in ${INbamARRAY[@]}; do echo "-I ${bam}";done) \
     -L $BED \
     -O $OUTvcf \
     --native-pair-hmm-threads 3 &
done

wait
