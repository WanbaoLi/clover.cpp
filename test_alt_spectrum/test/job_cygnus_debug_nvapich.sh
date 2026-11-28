#!/bin/bash
#PBS -A BRIDGE
#PBS -q gpu
#PBS -b 1
#PBS -l elapstim_req=1:00:00
#PBS -T mvapich
#PBS -v NQSV_MPI_VER=gdr/2.3.5/nvhpc20.9-cuda11.0

module load mvapich/${NQSV_MPI_VER}

export MV2_ENABLE_AFFINITY=0

cd ${PBS_O_WORKDIR}/

mpiexec ${NQSII_MPIOPTS} -np 4 -ppn 4 ${PBS_O_WORKDIR}/bridge.elf
