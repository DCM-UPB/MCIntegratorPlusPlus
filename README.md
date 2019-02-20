[![Build Status](https://travis-ci.com/DCM-UPB/MCIntegratorPlusPlus.svg?branch=master)](https://travis-ci.com/DCM-UPB/MCIntegratorPlusPlus)
[![Coverage Status](https://coveralls.io/repos/github/DCM-UPB/MCIntegratorPlusPlus/badge.svg?branch=travis)](https://coveralls.io/github/DCM-UPB/MCIntegratorPlusPlus?branch=travis)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4fb9f98862c7474b86c8ef88b501b454)](https://www.codacy.com/app/NNVMC/MCIntegratorPlusPlus?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DCM-UPB/MCIntegratorPlusPlus&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/dcm-upb/mcintegratorplusplus/badge)](https://www.codefactor.io/repository/github/dcm-upb/mcintegratorplusplus)

# MCIntegrator++

C++ Library for computing numerical integrals with the Monte Carlo method.
Some basic tools for finding the average and standard deviation of an array of data are included.

In `doc/` there is a user manual in pdf.

In `examples/` there are examples.



# Build the library

We use the CMake build system, so you need to have it on your system to build the library out of the box.
Then copy the file `config_template.sh` to `config.sh`, edit it to your liking and then simply execute the command

   `./build.sh`

Note that we build out-of-tree, so the compiled library and executable files can be found in the directories under `./build/`.
