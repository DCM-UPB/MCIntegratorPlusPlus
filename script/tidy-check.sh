#!/bin/sh
# Scans the library using clang-tidy, with a large set of checks enabled
# First argument can be used to add extra arguments, like -fix
for file in src/*.*pp test/main.cpp test/*/*.*pp examples/*/*.*pp; do
    clang-tidy -p build $file $1 -checks=bugprone-*,cppcoreguidelines-*,clang-analyzer-*,google-*,llvm-*,misc-*,modernize-*,mpi-*,performance-*,portability-*,readability-*,-cppcoreguidelines-owning-memory,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-google-runtime-int,-llvm-header-guard -header-filter=.* -- -Iinclude/ -Itest/common/ -DUSE_MPI=1
done
