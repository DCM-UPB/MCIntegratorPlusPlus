#!/bin/bash

# Library name
LIBNAME="mci"

# C++ compiler
CC="g++"

# C++ flags (std=c++11 is necessary)
FLAGS="-std=c++11 -Wall -Werror"

# Optimization flags
OPTFLAGS="-O3 -flto"

# Debuggin flags
DEBUGFLAGS="-g -O0"

#FFNN Library (used in ex3)
FFNN_FOLDER="/...../FeedForwardNeuralNetwork/lib/.libs"
IFFNN="-I${FFNN_FOLDER}/../../include/"
LFFNN="-L${FFNN_FOLDER}"
LIBNAMEFFNN="ffnn"
LIBFFNN="-lffnn"
