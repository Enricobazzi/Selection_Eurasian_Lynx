#!/bin/bash
#SBATCH -t 0-02:30:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

##############################################
## BWA mem and BAM post-processing launcher ##
##############################################

# This program runs BWA mem and then standard filtering steps for per-sample read alignment
# to the Canada lynx reference genome for all of the Lynx lynx samples.

# It will be run on the finis terrae II server of CESGA.

# Some of these samples have been sequenced in more than one lane, so the script will be split in two parts.
# This first part will generate a bam file for each fastqID (BARCODEID), sort it and add its read groups.
# The BARCODEID and the folder where the FASTQ is have to be given as a commands in the input while launching the script:

# ./LL_CNAG_LyCaRef_mapping_step1_cesga.sh <BARCODEID> <PATH/TO/FOLDER>


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

# BARCODES:
declare -A BARCODEID
# BARCODE LYNX_06_08:
BARCODEID+=(["C5RR4ACXX_7_14nf"]="c_ll_no_0075" ["C5T3DACXX_1_14nf"]="c_ll_no_0075" ["C3553ACXX_4_15nf"]="c_ll_no_0076" ["C5TMDACXX_5_15nf"]="c_ll_no_0076" ["C4592ACXX_3_16nf"]="c_ll_no_0077" ["C5TMDACXX_3_16nf"]="c_ll_no_0077" ["C3553ACXX_4_18nf"]="c_ll_no_0078" ["C5TMDACXX_1_18nf"]="c_ll_no_0078" ["C3553ACXX_4_19nf"]="c_ll_no_0079" ["C5TMDACXX_1_19nf"]="c_ll_no_0079" ["C3553ACXX_5_20nf"]="c_ll_no_0080" ["C5TMDACXX_2_20nf"]="c_ll_no_0080" ["C4592ACXX_3_21nf"]="c_ll_no_0081" ["C5TMDACXX_2_21nf"]="c_ll_no_0081" ["C3553ACXX_5_22nf"]="c_ll_no_0082" ["C5TMDACXX_3_22nf"]="c_ll_no_0082" ["C3553ACXX_5_23nf"]="c_ll_ki_0091" ["C5TMDACXX_5_23nf"]="c_ll_ki_0091" ["C3553ACXX_6_25nf"]="c_ll_ki_0092" ["C5TMDACXX_4_25nf"]="c_ll_ki_0092" ["C3553ACXX_6_27nf"]="c_ll_ki_0093" ["C5TMDACXX_4_27nf"]="c_ll_ki_0093" ["C3553ACXX_6_1nf"]="c_ll_ki_0094" ["C5T3DACXX_4_1nf"]="c_ll_ki_0094" ["C5RRAACXX_6_8nf"]="c_ll_ki_0095" ["C5RRAACXX_6_9nf"]="c_ll_ki_0096" ["C5TMDACXX_8_9nf"]="c_ll_ki_0096" ["C5TMDACXX_4_10nf"]="c_ll_ki_0097" ["C5TMUACXX_4_10nf"]="c_ll_ki_0097" ["C5TMUACXX_4_11nf"]="c_ll_ki_0098" ["C5TN1ACXX_7_11nf"]="c_ll_ki_0098" ["C5TMDACXX_8_12nf"]="c_ll_ki_0099" ["C5TMUACXX_5_12nf"]="c_ll_ki_0099" ["C5TMDACXX_8_13nf"]="c_ll_ki_0100" ["C5TMUACXX_5_13nf"]="c_ll_ki_0100" ["C5TMDACXX_7_14nf"]="c_ll_ki_0101" ["C5TMTACXX_8_14nf"]="c_ll_ki_0101" ["C5TMDACXX_7_15nf"]="c_ll_ki_0102" ["C5TMTACXX_8_15nf"]="c_ll_ki_0102" ["C5TMUACXX_2_1nf"]="c_lp_sm_0134" ["C5TN1ACXX_7_1nf"]="c_lp_sm_0134" ["C5TMUACXX_2_2nf"]="c_lp_do_0141" ["C5RRAACXX_5_3nf"]="c_lp_do_0144" ["C5T3DACXX_2_3nf"]="c_lp_do_0144" ["C5RRAACXX_5_4nf"]="c_lp_sm_0155" ["C5TMUACXX_3_5nf"]="c_lp_sm_0156" ["C5TN1ACXX_7_5nf"]="c_lp_sm_0156" ["C5TMUACXX_3_6nf"]="c_lp_sm_0161" ["C5TMTACXX_6_7nf"]="c_lp_do_0162" ["C5TN1ACXX_7_7nf"]="c_lp_do_0162" ["C5RR4ACXX_1_1nf"]="c_lp_do_0163" ["C5TN1ACXX_8_1nf"]="c_lp_do_0163" ["C5RR4ACXX_1_2nf"]="c_lp_sm_0206" ["C5TN1ACXX_8_2nf"]="c_lp_sm_0206" ["C5RR4ACXX_2_3nf"]="c_lp_sm_0208" ["C5T3DACXX_1_3nf"]="c_lp_sm_0208" ["C5RR4ACXX_2_4nf"]="c_lp_sm_0213" ["C5T3DACXX_1_4nf"]="c_lp_sm_0213" ["C5RR4ACXX_3_5nf"]="c_lp_sm_0226" ["C5T3DACXX_1_5nf"]="c_lp_sm_0226" ["C5RR4ACXX_3_6nf"]="c_lp_sm_0276" ["C5T3DACXX_1_6nf"]="c_lp_sm_0276" ["C5RR4ACXX_4_7nf"]="c_lp_do_0300" ["C5T3DACXX_2_7nf"]="c_lp_do_0300" ["C5RR4ACXX_4_8nf"]="c_lp_sm_0320" ["C5T3DACXX_1_8nf"]="c_lp_sm_0320" ["C5RR4ACXX_5_9nf"]="c_lp_sm_0325" ["C5T3DACXX_2_9nf"]="c_lp_sm_0325" ["C5RR4ACXX_5_10nf"]="c_lp_do_0333" ["C5TMDACXX_5_10nf"]="c_lp_do_0333" ["C5RR4ACXX_6_11nf"]="c_lp_do_0335" ["C5T3DACXX_2_11nf"]="c_lp_do_0335" ["C5RR4ACXX_6_12nf"]="c_lp_do_0444" ["C5T3DACXX_2_12nf"]="c_lp_do_0444" ["C5RR4ACXX_7_13nf"]="c_lp_sm_0450" ["C5T3DACXX_2_13nf"]="c_lp_sm_0450")
# BARCODE LYNX_09:
BARCODEID+=(["C6DUUANXX_2_12nf"]="c_ll_po_0001" ["C6DV6ANXX_1_12nf"]="c_ll_po_0001" ["C6DUUANXX_3_13nf"]="c_ll_po_0002" ["C6DV6ANXX_1_13nf"]="c_ll_po_0002" ["C6DUUANXX_3_14nf"]="c_ll_po_0003" ["C6DV6ANXX_1_14nf"]="c_ll_po_0003" ["C6DUUANXX_3_15nf"]="c_ll_po_0011" ["C6DV6ANXX_1_15nf"]="c_ll_po_0011" ["C6DUUANXX_2_16nf"]="c_ll_po_0014" ["C6DV6ANXX_1_16nf"]="c_ll_po_0014" ["C6DUUANXX_3_18nf"]="c_ll_po_0019" ["C6DV6ANXX_1_18nf"]="c_ll_po_0019" ["C6DUUANXX_4_19nf"]="c_ll_po_0105" ["C6DV6ANXX_1_19nf"]="c_ll_po_0105" ["C6DUUANXX_4_20nf"]="c_ll_po_0106" ["C6DV6ANXX_1_20nf"]="c_ll_po_0106" ["C6DUUANXX_2_21nf"]="c_ll_vl_0107" ["C6DV6ANXX_1_21nf"]="c_ll_vl_0107" ["C6DV6ANXX_7_5nf"]="c_ll_vl_0108" ["C6DUUANXX_2_22nf"]="c_ll_vl_0109" ["C6DUUANXX_1_23nf"]="c_ll_vl_0110" ["C6DV6ANXX_1_23nf"]="c_ll_vl_0110")
# BARCODES LYNX_PROYECT:
BARCODEID+=(["B09HCABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_1_0"]="c_lp_do_0153" ["B0B5KABXX_2_0"]="c_lp_do_0153" ["B09HCABXX_5_0"]="c_lp_do_0173" ["B09HCABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_6_0"]="c_lp_do_0173" ["B0B5KABXX_5_0"]="c_lp_do_0173" ["D0D6JABXX_4_0"]="c_lp_do_0443" ["D0D6JABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_3_0"]="c_lp_do_0443" ["B0999ABXX_4_0"]="c_lp_do_0443" ["B09HCABXX_3_0"]="c_lp_sm_0138" ["B09HCABXX_4_0"]="c_lp_sm_0138" ["B0B5KABXX_3_0"]="c_lp_sm_0138" ["B0B5KABXX_4_0"]="c_lp_sm_0138" ["C02CHABXX_1_0"]="c_lp_sm_0140" ["C02CHABXX_2_0"]="c_lp_sm_0140" ["C02CHABXX_4_0"]="c_lp_sm_0140" ["C02CHABXX_3_0"]="c_lp_sm_0140" ["D0D6JABXX_2_0"]="c_lp_sm_0185" ["D0D6JABXX_1_0"]="c_lp_sm_0185" ["B0999ABXX_2_0"]="c_lp_sm_0185" ["B0999ABXX_1_0"]="c_lp_sm_0185" ["D0D6JABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_6_0"]="c_lp_sm_0186" ["B0999ABXX_5_0"]="c_lp_sm_0186" ["D0D6JABXX_8_0"]="c_lp_sm_0298" ["D0D6JABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_7_0"]="c_lp_sm_0298" ["B0999ABXX_8_0"]="c_lp_sm_0298" ["C02CHABXX_6_0"]="c_lp_sm_0359" ["C02CHABXX_7_0"]="c_lp_sm_0359" ["C02CHABXX_5_0"]="c_lp_sm_0359" ["C02CHABXX_8_0"]="c_lp_sm_0359" ["B09HCABXX_8_0"]="h_lp_do_0007" ["B09HCABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_7_0"]="h_lp_do_0007" ["B0B5KABXX_8_0"]="h_lp_do_0007")
# BARCODE LYNX_13:
BARCODEID+=(["C9KJ0ANXX_1_22nf"]="c_ll_vl_0113" ["CA5U3ANXX_2_22nf"]="c_ll_vl_0113" ["C9KH1ANXX_7_25nf"]="c_ll_vl_0128" ["C9KH3ANXX_8_25nf"]="c_ll_vl_0128" ["C9KJ0ANXX_2_23nf"]="c_ll_vl_0132" ["CA5U3ANXX_2_23nf"]="c_ll_vl_0132" ["C9KH1ANXX_5_21nf"]="c_ll_ya_0138" ["C9KH3ANXX_7_21nf"]="c_ll_ya_0138" ["C9KH1ANXX_5_10nf"]="c_ll_ya_0139" ["C9KH3ANXX_7_10nf"]="c_ll_ya_0139" ["C9KH1ANXX_5_11nf"]="c_ll_ya_0140" ["C9KH3ANXX_7_11nf"]="c_ll_ya_0140" ["C9KH1ANXX_5_12nf"]="c_ll_ya_0142" ["C9KH3ANXX_7_12nf"]="c_ll_ya_0142" ["C9KH1ANXX_5_15nf"]="c_ll_ya_0143" ["C9KH3ANXX_7_15nf"]="c_ll_ya_0143" ["C9KH1ANXX_7_14nf"]="c_ll_ya_0145" ["C9KH3ANXX_7_14nf"]="c_ll_ya_0145" ["C9KH1ANXX_5_13nf"]="c_ll_ya_0147" ["C9KH3ANXX_7_13nf"]="c_ll_ya_0147" ["C9KJ0ANXX_4_19nf"]="c_ll_cr_0205" ["CA5U3ANXX_2_19nf"]="c_ll_cr_0205" ["CA2W6ANXX_4_16nf"]="c_ll_cr_0206" ["CA5U3ANXX_2_16nf"]="c_ll_cr_0206" ["C9KJ0ANXX_3_20nf"]="c_ll_cr_0207" ["CA5U3ANXX_2_20nf"]="c_ll_cr_0207" ["CA2W6ANXX_4_18nf"]="c_ll_cr_0208" ["CA5U3ANXX_2_18nf"]="c_ll_cr_0208" ["C9KH1ANXX_7_27nf"]="c_ll_cr_0209" ["C9KH3ANXX_8_27nf"]="c_ll_cr_0209")
# BARCODE LYNX_14:
BARCODEID+=(["C9KH6ANXX_5_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KN6ANXX_7_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KNWANXX_1_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KNWANXX_2_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KNWANXX_3_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KNWANXX_4_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KNWANXX_5_LYNX7-702ii5-4"]="h_lp_mt_0884" ["C9KH6ANXX_5_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KN6ANXX_7_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KNWANXX_1_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KNWANXX_2_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KNWANXX_3_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KNWANXX_4_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KNWANXX_5_LYNX7-710ii5-4"]="h_lp_mt_0885" ["C9KH6ANXX_5_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KN6ANXX_7_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KNWANXX_1_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KNWANXX_2_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KNWANXX_3_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KNWANXX_4_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KNWANXX_5_LYNX7-625ii5-2"]="h_lp_mt_0976" ["C9KH6ANXX_5_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KN6ANXX_7_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KNWANXX_1_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KNWANXX_2_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KNWANXX_3_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KNWANXX_4_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KNWANXX_5_LYNX7-709ii5-4"]="h_lp_mt_0978" ["C9KH6ANXX_5_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KN6ANXX_7_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KNWANXX_1_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KNWANXX_2_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KNWANXX_3_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KNWANXX_4_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KNWANXX_5_LYNX7-604ii5-2"]="h_lp_mt_0979" ["C9KH6ANXX_5_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KN6ANXX_7_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KNWANXX_1_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KNWANXX_2_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KNWANXX_3_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KNWANXX_4_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KNWANXX_5_LYNX7-706ii5-4"]="h_lp_mt_1025" ["CA3D2ANXX_3_LYNX7-706ii5-4"]="h_lp_mt_1025" ["C9KH6ANXX_5_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KN6ANXX_7_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KNWANXX_1_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KNWANXX_2_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KNWANXX_3_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KNWANXX_4_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KNWANXX_5_LYNX7-687ii5-4"]="h_lp_mt_1087" ["CA3D2ANXX_3_LYNX7-687ii5-4"]="h_lp_mt_1087" ["C9KH6ANXX_5_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KN6ANXX_7_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KNWANXX_1_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KNWANXX_2_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KNWANXX_3_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KNWANXX_4_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KNWANXX_5_LYNX7-700ii5-4"]="h_lp_mt_1117" ["CA3D2ANXX_3_LYNX7-700ii5-4"]="h_lp_mt_1117" ["C9KH6ANXX_5_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KN6ANXX_7_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KNWANXX_1_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KNWANXX_2_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KNWANXX_3_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KNWANXX_4_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KNWANXX_5_LYNX7-693ii5-4"]="h_lp_mt_1141" ["CA3D2ANXX_3_LYNX7-693ii5-4"]="h_lp_mt_1141" ["C9KH6ANXX_5_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KN6ANXX_7_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KNWANXX_1_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KNWANXX_2_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KNWANXX_3_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KNWANXX_4_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KNWANXX_5_LYNX7-686ii5-4"]="h_lp_mt_1272" ["C9KH6ANXX_5_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KN6ANXX_7_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KNWANXX_1_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KNWANXX_2_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KNWANXX_3_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KNWANXX_4_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KNWANXX_5_LYNX7-695ii5-4"]="h_lp_mt_1295" ["C9KH6ANXX_5_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KN6ANXX_7_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KNWANXX_1_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KNWANXX_2_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KNWANXX_3_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KNWANXX_4_LYNX7-708ii5-4"]="h_lp_mt_1305" ["C9KNWANXX_5_LYNX7-708ii5-4"]="h_lp_mt_1305" ["CA3D2ANXX_3_LYNX7-708ii5-4"]="h_lp_mt_1305")
# BARCODE CANDILES (6 samples that I added lately + (12/04/2017) I also added finally all the CANDILES here, because I need to have them in a separate array for the stats checkings):
BARCODEID+=(["6220RAAXX_lane3_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane4_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane6_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane7_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane8_sequence_0"]="c_lp_sm_0221" ["62AHEAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["621CYAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane1_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane5_sequence_0"]="c_lp_sm_0221" ["6220RAAXX_lane2_sequence_0"]="c_lp_sm_0221")
# BARCODE LYNX_15:
BARCODEID+=(["C9KH1ANXX_7_1nf"]="c_ll_to_0191" ["C9KH1ANXX_7_7nf"]="c_ll_la_0053" ["C9KH1ANXX_7_10nf"]="c_ll_la_0054" ["C9KH1ANXX_8_11nf"]="c_ll_la_0044" ["C9KH1ANXX_8_12nf"]="c_ll_tu_0157" ["C9KH1ANXX_8_13nf"]="c_ll_tu_0158" ["C9KH1ANXX_8_14nf"]="c_ll_tu_0159" ["C9KH1ANXX_8_15nf"]="c_ll_tu_0165" ["C9KH1ANXX_8_18nf"]="c_ll_tu_0166" ["C9KH1ANXX_8_19nf"]="c_ll_tu_0153" ["C9KH3ANXX_7_1nf"]="c_ll_to_0191" ["C9KH3ANXX_8_7nf"]="c_ll_la_0053" ["C9KH3ANXX_8_10nf"]="c_ll_la_0054" ["C9KH3ANXX_8_11nf"]="c_ll_la_0044" ["C9KH3ANXX_8_12nf"]="c_ll_tu_0157" ["C9KH3ANXX_8_13nf"]="c_ll_tu_0158" ["CA2W6ANXX_1_1nf"]="c_ll_og_0181" ["CA2W6ANXX_1_2nf"]="c_ll_og_0187" ["CA2W6ANXX_1_3nf"]="c_ll_ka_0184" ["CA2W6ANXX_1_4nf"]="c_ll_ka_0186" ["CA2W6ANXX_2_5nf"]="c_ll_ka_0189" ["CA2W6ANXX_2_7nf"]="c_ll_to_0190" ["CA2W6ANXX_2_8nf"]="c_ll_la_0047" ["CA2W6ANXX_4_6nf"]="c_ll_ka_0188" ["CA3D2ANXX_5_1nf"]="c_ll_og_0181" ["CA3D2ANXX_5_2nf"]="c_ll_og_0187" ["CA3D2ANXX_5_3nf"]="c_ll_ka_0184" ["CA3D2ANXX_5_4nf"]="c_ll_ka_0186" ["CA3D2ANXX_5_5nf"]="c_ll_ka_0189" ["CA3D2ANXX_5_6nf"]="c_ll_ka_0188" ["CA3D2ANXX_6_7nf"]="c_ll_to_0190" ["CA3D2ANXX_6_8nf"]="c_ll_la_0047" ["CA3D2ANXX_6_9nf"]="c_ll_la_0048" ["CA3D2ANXX_6_21nf"]="c_ll_la_0052" ["CAABGANXX_5_14nf"]="c_ll_tu_0159" ["CAABGANXX_5_15nf"]="c_ll_tu_0165" ["CAABGANXX_5_18nf"]="c_ll_tu_0166" ["CAABGANXX_5_19nf"]="c_ll_tu_0153" ["CAABGANXX_5_9nf"]="c_ll_la_0048" ["CAABGANXX_5_21nf"]="c_ll_la_0052")
# BARCODE LYNX_16:
BARCODEID+=(["CA2W6ANXX_4_23nf"]="h_ll_ba_0214" ["CA2W6ANXX_4_25nf"]="c_ll_ba_0216" ["CA2W6ANXX_4_27nf"]="h_ll_ba_0215" ["CA3D2ANXX_6_23nf"]="h_ll_ba_0214" ["CA3D2ANXX_6_25nf"]="c_ll_ba_0216" ["CA3D2ANXX_6_27nf"]="h_ll_ba_0215")
# BARCODE MACROGEN:
BARCODEID+=(["LC1"]="c_lc_zz_0001" ["LL112"]="c_ll_vl_0112" ["LL146"]="c_ll_ya_0146" ["LL212"]="c_ll_cr_0212" ["LL90"]="c_ll_ki_0090" ["LR1"]="c_lr_zz_0001")
# BARCODE LYNX_18:
BARCODEID+=(["CB84PANXX_7_152nfDI"]="c_ll_ur_0194" ["CB84PANXX_7_164nfDI"]="c_ll_ur_0195" ["CBCA6ANXX_7_176nfDI"]="c_ll_ur_0200" ["CBCA6ANXX_7_188nfDI"]="c_ll_ur_0196" ["CBCA6ANXX_8_105nfDI"]="c_ll_ur_0203" ["CBCA6ANXX_8_117nfDI"]="c_ll_ur_0199")
# BARCODE LYNX_20
BARCODEID+=(["HY5FMDSXX_1_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5FMDSXX_1_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5FMDSXX_2_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5FMDSXX_2_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5FMDSXX_3_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5FMDSXX_3_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5FMDSXX_4_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5FMDSXX_4_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5WLDSXX_1_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5WLDSXX_1_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5WLDSXX_2_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5WLDSXX_2_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5WLDSXX_3_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5WLDSXX_3_9UDI-idt-UMI"]="c_ll_ba_0224" ["HY5WLDSXX_4_21UDI-idt-UMI"]="c_ll_ca_0240" ["HY5WLDSXX_4_9UDI-idt-UMI"]="c_ll_ba_0224")
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
bwa mem $REF $INfastq_PATH/${fastqid}_1.fastq.gz $INfastq_PATH/${fastqid}_2.fastq.gz \
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
