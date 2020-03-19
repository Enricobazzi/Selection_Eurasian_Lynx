#!/bin/bash
scaffold=($(echo ${1}))

for i in {1..10}
 do
  echo "sbatching parallel_HaplotypeCaller for $scaffold w $i"
  sbatch parallel_HaplotypeCaller.sh $scaffold $i
done
