#!/bin/bash
#SBATCH -t 0-12:30:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

##############################################
## BWA mem and BAM post-processing launcher ##
##############################################

# This program runs BWA mem and then standard filtering steps for per-sample read alignment
# to the Canada lynx reference genome for all of the Lynx lynx samples from project LYNX_21.

# It will be run on the finis terrae II server of CESGA.

# Some of these samples have been sequenced in more than one lane, so the script
# will be split in two parts.
# This first part will generate a bam file for each fastqID (BARCODEID), sort it
# and add its read groups.
# The BARCODEID and the folder where the FASTQ is have to be given as a commands
# in the input while launching the script:

# ./LYNX21_LyCaRef_mapping_step1_cesga.sh <BARCODEID> <PATH/TO/FOLDER>


# The following softwares are used:
# BWA
# SAMtools
# Picard-tools
# GATK

module load gcc/6.4.0 bwa/0.7.17
module load gcc/6.4.0 samtools/1.8
module load gatk/3.7-0-gcfedb67
module load picard/2.18.14

###################################
## VARIABLE and PATHS definition ##
###################################

# BARCODE LYNX_21
BARCODEID+=(["HY5FMDSXX_1_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5FMDSXX_1_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5FMDSXX_1_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5FMDSXX_1_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5FMDSXX_1_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5FMDSXX_1_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5FMDSXX_1_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5FMDSXX_2_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5FMDSXX_2_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5FMDSXX_2_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5FMDSXX_2_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5FMDSXX_2_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5FMDSXX_2_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5FMDSXX_2_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5FMDSXX_3_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5FMDSXX_3_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5FMDSXX_3_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5FMDSXX_3_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5FMDSXX_3_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5FMDSXX_3_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5FMDSXX_3_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5FMDSXX_4_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5FMDSXX_4_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5FMDSXX_4_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5FMDSXX_4_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5FMDSXX_4_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5FMDSXX_4_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5FMDSXX_4_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5JLDSXX_1_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5JLDSXX_1_250UDI-idt-UMI"]="c_ll_ca_0249" ["HY5JLDSXX_1_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5JLDSXX_1_196UDI-idt-UMI"]="c_ll_ba_0226" ["HY5JLDSXX_1_208UDI-idt-UMI"]="c_ll_ba_0229" ["HY5JLDSXX_2_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5JLDSXX_2_250UDI-idt-UMI"]="c_ll_ca_0249" ["HY5JLDSXX_2_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5JLDSXX_2_196UDI-idt-UMI"]="c_ll_ba_0226" ["HY5JLDSXX_2_208UDI-idt-UMI"]="c_ll_ba_0229" ["HY5JLDSXX_3_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5JLDSXX_3_250UDI-idt-UMI"]="c_ll_ca_0249" ["HY5JLDSXX_3_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5JLDSXX_3_196UDI-idt-UMI"]="c_ll_ba_0226" ["HY5JLDSXX_3_208UDI-idt-UMI"]="c_ll_ba_0229" ["HY5JLDSXX_4_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5JLDSXX_4_250UDI-idt-UMI"]="c_ll_ca_0249" ["HY5JLDSXX_4_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5JLDSXX_4_196UDI-idt-UMI"]="c_ll_ba_0226" ["HY5JLDSXX_4_208UDI-idt-UMI"]="c_ll_ba_0229" ["HY5WLDSXX_1_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5WLDSXX_1_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5WLDSXX_1_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5WLDSXX_1_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5WLDSXX_1_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5WLDSXX_1_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5WLDSXX_1_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5WLDSXX_1_202UDI-idt-UMI"]="c_ll_ca_0244" ["HY5WLDSXX_1_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5WLDSXX_1_226UDI-idt-UMI"]="c_ll_ca_0247" ["HY5WLDSXX_1_238UDI-idt-UMI"]="c_ll_ca_0248" ["HY5WLDSXX_1_262UDI-idt-UMI"]="c_ll_ca_0252" ["HY5WLDSXX_1_274UDI-idt-UMI"]="c_ll_ca_0253" ["HY5WLDSXX_1_286UDI-idt-UMI"]="c_ll_ca_0254" ["HY5WLDSXX_1_203UDI-idt-UMI"]="c_ll_ca_0259" ["HY5WLDSXX_1_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5WLDSXX_2_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5WLDSXX_2_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5WLDSXX_2_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5WLDSXX_2_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5WLDSXX_2_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5WLDSXX_2_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5WLDSXX_2_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5WLDSXX_2_202UDI-idt-UMI"]="c_ll_ca_0244" ["HY5WLDSXX_2_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5WLDSXX_2_226UDI-idt-UMI"]="c_ll_ca_0247" ["HY5WLDSXX_2_238UDI-idt-UMI"]="c_ll_ca_0248" ["HY5WLDSXX_2_262UDI-idt-UMI"]="c_ll_ca_0252" ["HY5WLDSXX_2_274UDI-idt-UMI"]="c_ll_ca_0253" ["HY5WLDSXX_2_286UDI-idt-UMI"]="c_ll_ca_0254" ["HY5WLDSXX_2_203UDI-idt-UMI"]="c_ll_ca_0259" ["HY5WLDSXX_2_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5WLDSXX_3_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5WLDSXX_3_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5WLDSXX_3_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5WLDSXX_3_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5WLDSXX_3_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5WLDSXX_3_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5WLDSXX_3_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5WLDSXX_3_202UDI-idt-UMI"]="c_ll_ca_0244" ["HY5WLDSXX_3_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5WLDSXX_3_226UDI-idt-UMI"]="c_ll_ca_0247" ["HY5WLDSXX_3_238UDI-idt-UMI"]="c_ll_ca_0248" ["HY5WLDSXX_3_262UDI-idt-UMI"]="c_ll_ca_0252" ["HY5WLDSXX_3_274UDI-idt-UMI"]="c_ll_ca_0253" ["HY5WLDSXX_3_286UDI-idt-UMI"]="c_ll_ca_0254" ["HY5WLDSXX_3_203UDI-idt-UMI"]="c_ll_ca_0259" ["HY5WLDSXX_3_215UDI-idt-UMI"]="c_ll_ca_0260" ["HY5WLDSXX_4_33UDI-idt-UMI"]="c_ll_ba_0227" ["HY5WLDSXX_4_45UDI-idt-UMI"]="c_ll_ba_0228" ["HY5WLDSXX_4_57UDI-idt-UMI"]="c_ll_ba_0230" ["HY5WLDSXX_4_69UDI-idt-UMI"]="c_ll_ba_0233" ["HY5WLDSXX_4_81UDI-idt-UMI"]="c_ll_ca_0241" ["HY5WLDSXX_4_93UDI-idt-UMI"]="c_ll_ca_0242" ["HY5WLDSXX_4_10UDI-idt-UMI"]="c_ll_ca_0243" ["HY5WLDSXX_4_202UDI-idt-UMI"]="c_ll_ca_0244" ["HY5WLDSXX_4_214UDI-idt-UMI"]="c_ll_ca_0245" ["HY5WLDSXX_4_226UDI-idt-UMI"]="c_ll_ca_0247" ["HY5WLDSXX_4_238UDI-idt-UMI"]="c_ll_ca_0248" ["HY5WLDSXX_4_262UDI-idt-UMI"]="c_ll_ca_0252" ["HY5WLDSXX_4_274UDI-idt-UMI"]="c_ll_ca_0253" ["HY5WLDSXX_4_286UDI-idt-UMI"]="c_ll_ca_0254" ["HY5WLDSXX_4_203UDI-idt-UMI"]="c_ll_ca_0259" ["HY5WLDSXX_4_215UDI-idt-UMI"]="c_ll_ca_0260")

THREADS=24
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/Ref_Genome_LyCa/lc4.fa
SAMPLE=($(echo "${BARCODEID[${1}]}"))
fastqid=($(echo "${1}"))
OUT=/mnt/lustre/scratch/home/csic/ebd/jg2/LL_selection/LyCaRef_bams
INfastq_PATH=($(echo "${2}"))

# Mapping
echo " - Mapping $fastqid of ${SAMPLE} -"
bwa mem $REF $INfastq_PATH/${fastqid}_trimmomatic_s1_pe $INfastq_PATH/${fastqid}_trimmomatic_s2_pe \
-t $THREADS | samtools view -hbS -@ $THREADS - -o $OUT/${SAMPLE}_${fastqid}.LyCa_ref.bam

# Sorting
echo " - Sorting ${SAMPLE} -"
samtools sort -@ $THREADS $OUT/${SAMPLE}_${fastqid}.LyCa_ref.bam -o $OUT/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam \
&& rm $OUT/${SAMPLE}_${fastqid}.LyCa_ref.bam

# Adding READ Groups
echo " - Adding READ Groups of $fastqid of ${SAMPLE} -"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=$OUT/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam \
O=$OUT/${SAMPLE}_${fastqid}_LyCa_ref_sorted_rg.bam \
RGID=${fastqid} RGLB=${SAMPLE}_lib \
RGPL=Illumina RGPU=${fastqid} RGSM=${SAMPLE} \
VALIDATION_STRINGENCY=SILENT && rm $OUT/${SAMPLE}_${fastqid}.LyCa_ref.sorted.bam
