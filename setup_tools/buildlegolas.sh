#!/bin/bash

if [ "$1" == "clean" ]; then
  # to clean, $LEGOLASDIR must be set in order to find compiled files
  if [ -z "${LEGOLASDIR}" ]; then
    echo "Environment variable \$LEGOLASDIR is not defined, can't clean."
    exit
  fi
  # if local executable is found, remove it
  if [ -f legolas ]; then
    rm legolas
    echo "Local legolas executable removed."
  else
    echo "No local executable found, skipping."
  fi
  if [ -d build ] && [ ! "${LEGOLASDIR}" == "$(pwd)" ]; then
    rm -r build
    echo "Local build directory removed."
  else
    echo "No local build directory found, skipping."
  fi
  # clean build directory
  rm -rf "$LEGOLASDIR/build"
  echo "Build directory ${LEGOLASDIR}/build removed."
  # clean source directory executable
  rm -f "$LEGOLASDIR/legolas"
  echo "Executable in main legolas directory removed."
  exit
fi
if [ ! -f CMakeLists.txt ]; then
    echo "CMakeLists.txt not found! Can not start build process."
    exit
fi
mkdir build
cd build || exit
cmake ..
make -j 2
cd ..
