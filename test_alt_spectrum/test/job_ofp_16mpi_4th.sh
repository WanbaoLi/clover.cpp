#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L node=1
#PJM -L elapse=04:00:00
#PJM --mpi proc=16
#PJM --omp thread=4
#PJM -g xg18i014
#PJM -j

### This is a sample script file to submit on Oakforest-PACS.
### Change the bugjet group before submission.
### Also confirm that the parameters in main.yaml are
### corresponding to the run parameters aobve.
###                                    [27 Jul 2020 H.Matsufuru]

export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
#export KMP_AFFINITY=verbose,scatter

mpiexec.hydra -n ${PJM_MPI_PROC} ./bridge.elf
