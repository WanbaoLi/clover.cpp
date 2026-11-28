#!/bin/bash
#PBS -A BRIDGE
#PBS -q gpu
#PBS -b 1
#PBS -l elapstim_req=0:30:00
#PBS -T openmpi
###PBS -v NQSV_MPI_VER=4.1.7/nvhpc24.5-cuda12.4
#PBS -v NQSV_MPI_VER=nvhpc-hpcx-cuda12/24.5

module load openmpi/${NQSV_MPI_VER}

export MV2_ENABLE_AFFINITY=0

export OMP_NUM_THREADS=2


cd ${PBS_O_WORKDIR}/

mpiexec ${NQSII_MPIOPTS} -np 1 -npernode 1 -bind-to none ${PBS_O_WORKDIR}/bridge.elf
