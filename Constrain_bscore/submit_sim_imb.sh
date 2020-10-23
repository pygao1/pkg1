#!/bin/bash

## num_p is the number of ps for the simulation
num_theta=2
num_p=10
njobs=`expr $6 / $7 \* $num_p \* $num_theta`

qsub -cwd -e iotrash/ -o iotrash/ -t 1-$njobs -tc 20 ./call_sim_imb.sh $1 $2 $3 $4 $5 $6 $7