#!/bin/bash

# After using this script it is necessary to run again the build.sh script
# for generating again the library with the optimization flags

OS_NAME=$(uname)
echo "The Operating System is: "${OS_NAME}  # here we consider only Linux and Darwin (Mac Os X)
echo

source ../config.sh
DEBUGFLAGS="-g -O0"

\rm -f exe
\rm -f *.o
\rm -f ../src/*.o
\rm -f ../*.so

# Build the library using the debugging flags
echo "Build the library using the debugging flags . . ."
cd ..
LIBFOLDER=$(pwd)
cd src

case ${OS_NAME} in
      "Linux")
         echo "$CC $DEBUGFLAGS -std=c++11 -shared -o lib${LIBNAME}.so *.o"
         $CC $DEBUGFLAGS $FLAGS -shared -o lib${LIBNAME}.so *.o
         echo "$CC $DEBUGFLAGS -std=c++11 -shared -o lib${LIBNAME}.so *.o"
         $CC $DEBUGFLAGS $FLAGS -shared -o lib${LIBNAME}.so *.o
         mv lib*.so ../
         ;;
      "Darwin")
         echo "$CC $DEBUGFLAGS $FLAGS -dynamiclib -c *.cpp"
         $CC $DEBUGFLAGS $FLAGS -dynamiclib -c *.cpp
         echo "$CC $DEBUGFLAGS $FLAGS -dynamiclib -install_name ${LIBFOLDER}/lib${LIBNAME}.dylib -o lib${LIBNAME}.dylib *.o"
         $CC $DEBUGFLAGS $FLAGS -dynamiclib -install_name ${LIBFOLDER}/lib${LIBNAME}.dylib -o lib${LIBNAME}.dylib *.o
         mv lib*.dylib ../
         ;;
      *)
         echo "The detected operating system is not between the known ones (Linux and Darwin)"
         ;;
esac

cd ../debug
echo "Rebuilt the library with the debugging flags"


## Build the debugging main executable
echo ""
echo "Build the executable . . ."

case ${OS_NAME} in
      "Linux")
         echo "$CC $FLAGS $DEBUGFLAGS -Wall -I$(pwd)/../src/ -c *.cpp"
         $CC $FLAGS $DEBUGFLAGS -Wall -I$(pwd)/../src/ -c *.cpp
         echo "$CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -Wl,-rpath=$(pwd)/../ -o exe *.o -l${LIBNAME}"
         $CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -Wl,-rpath=$(pwd)/../ -o exe *.o -l${LIBNAME}
         ;;
      "Darwin")
         echo "$CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -c *.cpp"
         $CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -c *.cpp
         echo "$CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -o exe *.o -l${LIBNAME}"
         $CC $FLAGS $DEBUGFLAGS -I$(pwd)/../src/ -L$(pwd)/../ -o exe *.o -l${LIBNAME}
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

