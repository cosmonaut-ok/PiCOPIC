#!/bin/bash

# HOSTS=${1}
PROCS=${1}

if test -z ${PROCS}; then
    PROCS=2
fi

TMPFILE=$(mktemp /tmp/.picopicXXXX)

echo "localhost slots=${PROCS}" > ${TMPFILE}

rm -f ./data.h5

export OMP_NUM_THREADS=8
mpirun -display-map -x OMP_NUM_THREADS=8 --hostfile ${TMPFILE} -n $PROCS -H localhost --mca btl_tcp_if_exclude docker0 ./PiCoPiC ## xterm -e gdb ./PiCoPiC

rm -f ${TMPFILE}
