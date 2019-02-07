#!/bin/sh

source ./config.sh
mkdir -p build
cd build
cmake -DUSER_CXX_FLAGS="${CXX_FLAGS}" ..
make
cd ..
