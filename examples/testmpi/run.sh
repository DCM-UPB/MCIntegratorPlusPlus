#!/bin/bash
source ../../config.sh
OS_NAME=$(uname)

# Delete old compiled files
\rm -f exe
\rm -f *.o

#runtime dynamic library path
RPATH="$(pwd)/../.."

FLAGS_TO_USE=$OPTFLAGS

# Build the main executable
echo "$MPICC $FLAGS $FLAGS_TO_USE -I$(pwd)/../../src/ -c *.cpp"
$MPICC $FLAGS $FLAGS_TO_USE -Wall -I$(pwd)/../../src/ -c *.cpp

case ${OS_NAME} in
    "Darwin")
        echo "$MPICC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. -o exe *.o -l${LIBNAME}"
        $MPICC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. -o exe *.o -l${LIBNAME}
        ;;
    "Linux")
        echo "$MPICC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}"
        $MPICC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME}
        ;;
esac

echo "Rebuilt the executable file"
echo ""
echo ""

# Run the debugging executable
echo "Ready to run!"
echo ""
echo "--------------------------------------------------------------------------"
echo ""
echo ""
echo ""
mpirun -np 4 exe
