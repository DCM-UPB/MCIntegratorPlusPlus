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
echo "$CC $FLAGS $FLAGS_TO_USE -I$(pwd)/../../src/ ${IFFNN} -c *.cpp"
$CC $FLAGS $FLAGS_TO_USE -Wall -I$(pwd)/../../src/ ${IFFNN} -c *.cpp

case ${OS_NAME} in
   "Darwin")
      echo "$CC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. ${LFFNN} -o exe *.o -l${LIBNAME} ${LIBFFNN}"
      $CC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. ${LFFNN} -o exe *.o -l${LIBNAME} ${LIBFFNN}
      ;;
   "Linux")
      echo "$CC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. ${LFFNN} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} ${LIBFFNN}"
      $CC $FLAGS $FLAGS_TO_USE -L$(pwd)/../.. ${LFFNN} -Wl,-rpath=${RPATH} -o exe *.o -l${LIBNAME} ${LIBFFNN}
      ;;
esac


./exe
