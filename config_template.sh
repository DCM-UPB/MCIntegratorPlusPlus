#!/bin/sh

#C++ compiler
CXX_COMPILER="g++"

# C++ flags
CXX_FLAGS="-O3 -flto -march=native -Wall -Wno-unused-function"

# compile with MPI
USE_MPI=0

# add coverage flags
USE_COVERAGE=0
