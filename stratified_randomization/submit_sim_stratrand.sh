#!/bin/bash

##test for theta0, theta1, and use simple vs stratified randomization
##so in total 4 combinations
njobs=`expr $6 / $7 \* 4`

qsub -cwd -e iotrash/ -o iotrash/ -t 1-$njobs -tc 20 ./call_sim_stratrand.sh $1 $2 $3 $4 $5 $6 $7 