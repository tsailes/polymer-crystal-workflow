#!/usr/bin/env bash
# Source this file for the verified CPU environment; it starts no calculation.
export PATH="/opt/pe_dft/envs/pe-gpaw-20261020/bin:${PATH}"
export GPAW_SETUP_PATH="/opt/pe_dft/envs/pe-gpaw-20261020/lib/python3.12/site-packages/gpaw_data/setups"
export GPAW_MPI_BACKEND=cgpaw
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# WSL distribution is run as root.
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
