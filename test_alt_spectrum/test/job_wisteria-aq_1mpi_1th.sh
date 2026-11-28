#!/bin/bash
#PJM -L rscgrp=share-debug
#PJM -L gpu=1
#PJM -L elapse=0:30:00
#PJM --omp thread=1
####PJM -g wo22i060
#PJM -g ga42
#PJM -j

module load nvidia
module load ompi

mpiexec ${NQSII_MPIOPTS} -np 1 -npernode 1 ./bridge.elf

