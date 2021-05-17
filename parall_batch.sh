#!/bin/bash

# $1 = size of the bootstrap chunks
# $2 = nb of bootstrap chunks
# $3 = nb of threads
# overall, $1 * $2 = nb of bootstraps

# trick, pass on 4 25 1 on a first shot to build the rarefied communities
# then on a 2 50 2 to pass the analyses to not overload the RAM

chk_sz=$1

for i in `seq 1 $2`
do
  for j in `seq $((i*chk_sz-(chk_sz-1))) $((i*chk_sz))`
  do
    Rscript 01_analyses_on_boot.r $j 10000 $3 F &
  done
  wait
done
