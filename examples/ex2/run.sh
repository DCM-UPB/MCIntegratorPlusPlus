#!/bin/sh
cd ../../build/examples
mpirun -np $1 ./ex2.exe
