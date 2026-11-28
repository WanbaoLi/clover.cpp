#!/bin/bash
#PJM -L rscunit=fx
#PJM -L rscgrp=fx-debug
#PJM -L node=1
#PJM -L elapse=1:00:00
#PJM --mpi proc=1
#PJM -j

export OMP_NUM_THREADS=48
export FLIB_BARRIER=HARD

mpiexec ./bridge.elf

