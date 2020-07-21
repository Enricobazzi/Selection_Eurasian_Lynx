#!/bin/bash

WD=/home/ebazzicalupo/BayPass

for n in {1..50}
 do
  echo "dataset ${n}"
  g_baypass -npop 8 -gfile $WD/AlleleCounts/all.allelecounts.${n} -efile $WD/Covariate_Data/${1}_data.txt \
  -auxmodel -scalecov -omegafile $WD/OutPut/CORE_1_mat_omega.out -outprefix $WD/OutPut/AUX_${1}_${n}
done
