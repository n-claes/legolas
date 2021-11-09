#!/bin/bash

cd $LEGOLASDIR

# Build Legolas
if [[ ! -d build ]]; then
    mkdir build
fi
cd build
cmake ..
make -j 2
cd ..

# Run regression tests
cd tests/regression_tests
pytest -v regression.py test*

