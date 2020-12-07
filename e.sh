#!/bin/sh

mpic++ mpi_send_recv.cpp -I./lib/mpp/include/

echo "RUN"
mpirun -n 4 ./a.out
echo "END"
