#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

OS_NAME=$(uname)
source ../../../config.sh

\rm -f exe
\rm -f *.o

ROOT_FOLDER=$(dirname $(dirname $(dirname $(pwd))))


## Build the debugging main executable
$CC $FLAGS $DEBUGFLAGS -Wall -I${ROOT_FOLDER}/src/ -c *.cpp

case ${OS_NAME} in
    "Linux")
        $CC $FLAGS $DEBUGFLAGS -L${ROOT_FOLDER} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}
        ;;
    "Darwin")
        $CC $FLAGS $DEBUGFLAGS -L${ROOT_FOLDER} -o exe *.o -l${LIBNAME}
        ;;
    *)
        echo "The detected operating system is not between the known ones (Linux and Darwin)"
        ;;
esac

# Run the debugging executable
./exe
