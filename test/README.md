# LEGEND OF THE UNIT TESTS

Use `./run.sh` inside the test directory to run the check program and unit tests
with valgrind or use `make test` inside the build directory, to run unit tests without valgrind.

## Unit Test 1

`ut1/`: Check that the accumulators and estimators are working correctly


## Unit Test 2

`ut2/`: check that calling the integrator again without doing the findMRT2step and initialDecorrelation gives the same results.


## Unit Test 3

`ut3/`: Like ut2, but with one one-dimensional and one three-dimensional observable.


## Unit Test 4

`ut4/`: Like ut2, but checking with fixed number of findMRT2step/decorrelation steps and fixed blocking.


## Unit Test 5

`ut5/`: Like ut3, but testing with all the available trial moves (including elementary updates in sampling fun).
