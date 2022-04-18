#!/bin/bash

cd $LEGOLASDIR

# Build Legolas
if [[ ! -d build ]]; then
    mkdir build
fi
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DDebug=ON
make -j 2
cd ..

# Run regression tests
cd tests/regression_tests
pytest -v regression.py test*

