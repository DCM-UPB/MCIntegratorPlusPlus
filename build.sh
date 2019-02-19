#!/bin/sh
. ./config.sh
mkdir -p build
cd build
cmake -DCMAKE_CXX_COMPILER="${CXX_COMPILER}" -DUSER_CXX_FLAGS="${CXX_FLAGS}" ..
make
