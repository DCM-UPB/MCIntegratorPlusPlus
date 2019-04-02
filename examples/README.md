# LEGEND OF THE EXAMPLES

Make sure the examples are compiled, by running `./build.sh` in the project root folder.
Execute an example by switching into one of the example folders and running `./run.sh`.
Note that the actual example executables reside inside the `build/examples/` folder under the project's root.

## Common

In the folder `common/` you will find a file named `ExampleFunctions.hpp` that contains some simple sampling
function and observable implementations that are shared between the examples.


## Example 1

`ex1/`: integration of a 1-dimensional quadratic function, without and with a sampling function.


## Example 2

`ex2/`: like ex2, but using the MPI wrapper for parallel integration. Pass the wanted number of threads to run.sh.
