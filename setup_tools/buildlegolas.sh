if [ ! -f CMakeLists.txt ]; then
    echo "CMakeLists.txt not found! Can not start build process."
    exit
fi
if [ "$1" == "clean" ]
then
  rm -f legolas
  echo "legolas executable removed."
  if [[ -z "${LEGOLASDIR}" ]]; then
    echo "Environment variable \$LEGOLASDIR is not defined, can't clean."
  else
    rm -rf "$LEGOLASDIR/build"
    echo "build directory removed."
  fi
  exit
fi
mkdir build
cd build || exit
cmake ..
make -j 2
cd ..
rm -rf build
