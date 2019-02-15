#!/bin/bash
ORIGDIR=$(pwd)
cd ../../build/examples
mpirun -np $1 ./ex2.exe
cd "${ORIGDIR}"
