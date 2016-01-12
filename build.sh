#!/bin/bash

OS_NAME=$(uname)
echo "The Operating System is: "${OS_NAME}  # here we consider only Linux and Darwin (Mac Os X)

source config.sh
ACTUAL_FOLDER=$(pwd)

case ${OS_NAME} in
   "Linux")
      \rm -f *.so
      cd src/
         \rm -f *.o *.so
         echo ""
         echo "$CC $FLAGS $OPTFLAGS -fpic -c *.cpp"
         $CC $FLAGS $OPTFLAGS -fpic -c *.cpp
         echo "$CC $FLAGS $OPTFLAGS -shared -o lib${LIBNAME}.so *.o"
         $CC $FLAGS $OPTFLAGS -shared -o lib${LIBNAME}.so *.o
         echo ""
         mv lib${LIBNAME}.so ../
      cd ..
      ;;
   "Darwin")
      \rm -f *.so
      cd src/
         \rm -f *.o *.so
         echo ""
         echo "$CC $FLAGS $OPTFLAGS -dynamiclib -c *.cpp"
         $CC $FLAGS $OPTFLAGS -dynamiclib -c *.cpp
         echo "$CC $FLAGS $OPTFLAGS -shared -o lib${LIBNAME}.so *.o"
         $CC $FLAGS $OPTFLAGS -dynamiclib -install_name ${ACTUAL_FOLDER}/lib${LIBNAME}.dylib -o lib${LIBNAME}.dylib *.o
         echo ""
         mv lib${LIBNAME}.dylib ../
      cd ..
      ;;
   *)
      echo "The detected operating system is not between the known ones (Linux and Darwin)"
      ;;
esac

echo
echo "Library ready!"
echo
echo "Help, how can I use it?"
case ${OS_NAME} in
   "Linux")
      echo "1)   $CC $FLAGS -I$(pwd)/src/ -c example.cpp"
      echo "     $CC $FLAGS -L$(pwd) example.o -l${LIBNAME}" 
      echo "2)   $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
      ;;
   "Darwin")
      echo "1)   $CC $FLAGS -I$(pwd)/src/ -c example.cpp"
      echo "     $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.o -l${LIBNAME}"
      echo "2)   $CC $FLAGS -I$(pwd)/src/ -L$(pwd) example.cpp -l${LIBNAME}"
      ;;
esac
