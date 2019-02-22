# LEGEND OF THE UNIT TESTS

Use `./run.sh` inside the test directory to run the check program and unit tests
with valgrind or use `make test` inside the build directory, to run unit tests without valgrind.


## Unit Test 1

`ut1/`: check that calling the integrator again without doing the findMRT2step and initialDecorrelation gives the same results.


## Unit Test 2

`ut2/`: Like ut1, but with one one-dimensional and one three-dimensional observable.


## Unit Test 3

`ut3/`: Like ut1, but checking with fixed number of findMRT2step/decorrelation steps and fixed blocking technique
