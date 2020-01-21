#!/bin/bash

set -ie

if [[ -z $1 ]]; then echo -e "USAGE\n\t$0 BUILD_NAME" && exit 256; fi
BUILD_NAME=$1

module purge
module load cmake
module load cuda/10.0
module load gcc/7.1.0


rm -rf $BUILD_NAME
mkdir -p $BUILD_NAME && cd $BUILD_NAME

cmake ../src/ \
-DCMAKE_C_COMPILER=$(which gcc)  \
-DCMAKE_CXX_COMPILER=$(which g++)  \
-DCUDA_INCLUDE_DIRS=${CUDA_HOME}/include \
-DBUILD_TESTING=OFF \


make -j 12