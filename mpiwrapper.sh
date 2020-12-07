#!/bin/bash

set -e

# HOSTS=${1}
PROCS=${1}

if test -z ${PROCS}; then
    PROCS=2
fi

TMPFILE=$(mktemp /tmp/.picopicXXXX)

PROCESSORS=$(cat /proc/cpuinfo  | grep processor | tail -n1 | awk '{print $3}')

echo "localhost slots=${PROCS}" > ${TMPFILE}

rm -f ./data.h5
make
export OMP_NUM_THREADS=$(( ($PROCESSORS + 1) / $PROCS))
mpirun -display-map -x OMP_NUM_THREADS=8 --hostfile ${TMPFILE} -n $PROCS -H localhost --mca btl_tcp_if_exclude docker0 ./PiCoPiC ## xterm -e gdb ./PiCoPiC
# mpirun -display-map -x OMP_NUM_THREADS=8 --hostfile ${TMPFILE} -n $PROCS -H localhost --mca btl_tcp_if_exclude docker0 valgrind --tool=massif --depth=10 ./PiCoPiC
rm -f ${TMPFILE}
