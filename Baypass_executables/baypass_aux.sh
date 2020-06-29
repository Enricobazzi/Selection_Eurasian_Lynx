#!/bin/bash

WD=/home/ebazzicalupo/BayPass

g_baypass -npop 8 -gfile $WD/AlleleCounts/all.allelecounts.${1} -efile $WD/Covariate_Data/${2}_data.txt \
-auxmodel -omegafile $WD/OutPut/CORE_1_mat_omega.out -outprefix $WD/OutPut/AUX_${2}_${1}
