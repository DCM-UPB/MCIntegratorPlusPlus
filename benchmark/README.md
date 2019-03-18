# Benchmarks

This directory contains benchmarks to test the performance of certain parts of the library.
The `common` subfolder contains common source code and script files that are used by the individual benchmarks in `bench_*` folders.

A conda environment that can be used to run the python script is available as `conda-env.yml`

Currently there are the following benchmarks:

   `bench_estimators`: Benchmark of MCI's estimator implementations, for different settings.
   `bench_integrate_mixed`: Benchmark of MC integration (uni-all-moves) in 3D, for a fast PDF and a small mix of observables.
   `bench_throughput_nmc`: Benchmark of maximal MC sampling throughput in 1D, depending on NMC, using near-zero cost PDF&observable.
   `bench_throughput_ndim`: Like the previous, but with fixed NMC and varying number of dimensions.
   `bench_throughput_3G`: Like the first throughput bench, but a single run of 3 Giga-Samples (also a test regarding integer overflow).

# Using the benchmarks

Just provide the script `run.sh` the desired benchmark's name, e.g.:
   `./run.sh bench_estimators`

Alternatively, you can run all benchmarks sequentially by calling:
   `./run_all.sh`

The benchmark results will be written to a file named `benchmark_new.out` under the respective benchmark folder.
You may visualize the result by entering that directory and using:
   `python plot.py benchmark_new.out`

To let the plot compare the new result versus an older one, you have to provide the old output file like:
   `python plot.py benchmark_old.out benchmark_new.out`.

You may also change new/old to more meaningful labels, anything like benchmark_*.out is allowed (except extra _ or . characters). The
provided labels will be used automatically to create the plot legends.


# Profiling

If you want to performance profile the library under execution of a benchmark,
you just need to provide gperftools's libprofiler.so library to `run_prof.sh` as second argument, e.g.:
   `./run_prof.sh bench_actfs_ffprop /usr/lib/libprofiler.so`

Note that this script does not save any benchmark results.
