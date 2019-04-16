[![Build Status](https://travis-ci.com/DCM-UPB/MCIntegratorPlusPlus.svg?branch=master)](https://travis-ci.com/DCM-UPB/MCIntegratorPlusPlus)
[![Coverage Status](https://coveralls.io/repos/github/DCM-UPB/MCIntegratorPlusPlus/badge.svg?branch=master)](https://coveralls.io/github/DCM-UPB/MCIntegratorPlusPlus?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4fb9f98862c7474b86c8ef88b501b454)](https://www.codacy.com/app/NNVMC/MCIntegratorPlusPlus?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DCM-UPB/MCIntegratorPlusPlus&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/dcm-upb/mcintegratorplusplus/badge)](https://www.codefactor.io/repository/github/dcm-upb/mcintegratorplusplus)

# MCIntegrator++

C++ Library for computing numerical integrals with the Monte Carlo method. Includes some convenient
(optional) methods for automatic step calibration, decorrelation and error estimation. Provides a simple
interface to execute the MC integration in parallel, via Message Passing Interface (MPI).

In `doc/` there is a user manual in pdf (not accurate for current master!) and a config for doxygen.

In `examples/`, `test/` and `benchmark/` there are examples, tests and benchmarks for the library.

In `script/` we collect useful scripts and in `res/` we provide resources, like a true random seed file.


Some subdirectories come with an own `README.md` file which provides further information.


# Supported Systems

Currently, we automatically test the library on Arch Linux (GCC 8) and MacOS (with clang as well as brewed GCC 8).
However, in principle any system with C++11 supporting compiler should work.


# Requirements

- CMake, to use our build process
- (optional) a MPI implementation, to use parallelized integration
- (optional) valgrind, to run `./run.sh` in `test/`
- (optional) pdflatex, to compile the tex file in `doc/`
- (optional) doxygen, to generate doxygen documentation in `doc/doxygen`


# Build the library

Copy the file `config_template.sh` to `config.sh`, edit it to your liking and then simply execute the command

   `./build.sh`

Note that we build out-of-tree, so the compiled library and executable files can be found in the directories under `./build/`.


# First steps

You may want to read `doc/user_manual.pdf` to get a quick overview of the libraries functionality. However, it is not guaranteed to be perfectly up-to-date and accurate. Therefore, the best way to get your own code started is by studying the examples in `examples/`. See `examples/README.md` for further guidance.


# Multi-threading: MPI

This library supports multi-threaded MC integration with a distributed-memory paradigm, thanks to Message Passing interface (MPI).

To be able to use this feature, just compile the library with a MPI implementation present on your system. The header `MPIMCI.hpp` provides convenient functions
for using MCI++ with MPI. For example usage, look into example ex2.
