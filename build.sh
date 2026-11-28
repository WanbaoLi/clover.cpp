#!/bin/bash
set -e
cd "$(dirname "$0")"

export PATH="/public2/software/mpi/openmpi/4.1.4/gnu/bin:/usr/bin:$HOME/miniconda3/bin:$PATH"

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build build -j "$(nproc)"
