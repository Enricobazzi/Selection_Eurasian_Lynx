#!/bin/bash

WD=/home/ebazzicalupo/BayPass

g_baypass -npop 8 -gfile $WD/AlleleCounts/all.allelecounts.${1} -outprefix $WD/OutPut/CORE_${1}
