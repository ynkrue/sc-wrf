#!/bin/bash

source /home/yrfenach/Fromager/setup-env.sh
fgr mpi load openmpi-5.0.10-gcc-15.2.0
fgr compiler load gcc-15.2.0

NC=/scratch/yrfenach/.fromager/cellars/netcdf

export NETCDF="$NC"
export PNETCDFPATH="$NC"
export HDF5="$NC"
export WRF_CUDA_PATH=/scratch/yrfenach/.fromager/cellars/system

WRF_DIR=/scratch/yrfenach/applications/sc26/wrf
cd $WRF_DIR/sc-wrf

./clean -aa

./configure

./compile em_real

SHARE_DIR=$WRF_DIR/share/gcc
mkdir -p "$SHARE_DIR" && cp -rL sc-wrf/run/. "$SHARE_DIR/" && for e in wrf real ndown tc; do [ -s "$SHARE_DIR/$e.exe" ] || { echo "FATAL: $e.exe missing/empty" >&2; exit 1; }; done
