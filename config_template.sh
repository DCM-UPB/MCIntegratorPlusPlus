#!/bin/bash

# Library name
LIBNAME="mci"

# C++ compiler
CC="g++"

# MPI compiler wrapper
MPICC="mpic++"

# C++ flags (std=c++11 is necessary)
FLAGS="-std=c++11 -Wall -Werror"

# Optimization flags
OPTFLAGS="-O3 -flto"

# Debuggin flags
DEBUGFLAGS="-g -O0"
