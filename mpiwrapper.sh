#!/bin/bash

HOSTS=${1}
PROCS=${2}

export OMP_NUM_THREADS=4
mpirun -display-map -n $PROCS -H $HOSTS --mca btl_tcp_if_exclude docker0 ./PiCoPiC
