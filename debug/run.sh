#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

OS_NAME=$(uname)
echo "The Operating System is: "${OS_NAME}  # here we consider only Linux and Darwin (Mac Os X)
echo

source ../config.sh
DEBUGFLAGS="-g -O0"

#runtime dynamic library path
RPATH="$(pwd)/.."

\rm -f exe
\rm -f *.o

# Build the library using the debugging flags
echo "Build the library using the debugging flags . . ."
cd ..
   LIBFOLDER=$(pwd)
   ./build_debug_library.sh
cd debug
echo "Rebuilt the library with the debugging flags"


## Build the debugging main executable
echo ""
echo "Build the executable . . ."

echo "$CC $FLAGS $DEBUGFLAGS -Wall -I$(pwd)/../src/ -c *.cpp"
$CC $FLAGS $DEBUGFLAGS -Wall -I$(pwd)/../src/ -c *.cpp

case ${OS_NAME} in
      "Darwin")
         echo "$CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -o exe *.o -l${LIBNAME}"
         $CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -o exe *.o -l${LIBNAME}
         ;;
      "Linux")
         echo "$CC $FLAGS $DEBUGFLAGS -L$(pwd)/../ -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}"
         $CC $FLAGS $DEBUGFLAGS -L$(pwd)/../ -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}
         ;;
      *)
         echo "The detected operating system is not between the known ones (Linux and Darwin)"
         ;;
esac
echo "Rebuilt the debugging executable"
echo ""
echo ""

# Run the debugging executable
valgrind --track-origins=yes ./exe

