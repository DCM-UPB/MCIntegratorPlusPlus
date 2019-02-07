
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
