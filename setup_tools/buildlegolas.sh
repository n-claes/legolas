rm -rf build
rm -f legolas
if [ "$1" == "clean" ]
then
  echo "build directory removed."
  echo "legolas executable removed."
  exit
fi
mkdir build
cd build || exit
cmake ..
make -j 2
cd ..
