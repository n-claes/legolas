
![legolas-logo](docs/assets/images/logo_legolas_1280_trans.png)
**L**arge **E**igensystem **G**enerator for **O**ne-dimensional p**LAS**mas

### Table of contents
1. [CI status](#ci-status)
2. [Using Legolas](#using-legolas)
3. [Obtaining Legolas](#obtaining-legolas)
4. [Requirements](#requirements)
5. [Compilation](#compilation)
6. [Contents](#directory-contents)

## CI status
- Branch `master` <br>
  [![regression](https://github.com/n-claes/legolas/workflows/regression/badge.svg?branch=master)](https://github.com/n-claes/legolas/actions?query=workflow%3Aregression+branch%3Amaster)
  [![legolas-core](https://github.com/n-claes/legolas/workflows/legolas-core/badge.svg?branch=master)](https://github.com/n-claes/legolas/actions?query=workflow%3Alegolas-core+branch%3Amaster)
  [![pylbo-core](https://github.com/n-claes/legolas/workflows/pylbo-core/badge.svg?branch=master)](https://github.com/n-claes/legolas/actions?query=workflow%3Apylbo-core+branch%3Amaster)
- Branch `develop` <br>
  [![regression](https://github.com/n-claes/legolas/workflows/regression/badge.svg?branch=develop)](https://github.com/n-claes/legolas/actions?query=workflow%3Aregression+branch%3Adevelop)
  [![legolas-core](https://github.com/n-claes/legolas/workflows/legolas-core/badge.svg?branch=develop)](https://github.com/n-claes/legolas/actions?query=workflow%3Alegolas-core+branch%3Adevelop)
  [![pylbo-core](https://github.com/n-claes/legolas/workflows/pylbo-core/badge.svg?branch=develop)](https://github.com/n-claes/legolas/actions?query=workflow%3Apylbo-core+branch%3Adevelop)

## Using Legolas
:warning: **Please note** :warning:

Legolas is the result of months and months of developing, testing, fixing issues, testing again,
thinking bugs are fixed, further development, discovering that bugs weren't fixed, headscratching, testing again, etc.
In short, a typical development process of a brand new code. Since this took (and still takes) a lot of effort and time,
we therefore kindly ask that _**the first published peer-reviewed paper from applying Legolas is done in co-authorship with at least one of the original authors**_.
Since the code is brand new we would like to know how it is used and provide guidance if possible.

## Obtaining Legolas
You can obtain the latest stable version of the code through
```bash
git clone https://github.com/n-claes/legolas.git
cd legolas
git checkout master
```
We ensure that the `master` branch always contains a stable release. If you would
like to check out the latest development goodies of the code, you can checkout the `develop` branch instead
and build from that one. Our typical workflow is branching off of `develop`, extend or create a feature and merge it back
into `develop` when the automated tests are passing. Every now and then we merge `develop` into `master` and bump the
release version. We hence (try to) ensure that `develop` is stable, but beware that this may not always be the case.

Also note that due to this development workflow there will be multiple branches present in the repository at any
given time. Feel free to check these out, but note that we will frequently rebase those so know that merge conflicts may
be common when pulling updates from those branches.

## Requirements
We make extensive use of submodules (Fortran 2008) and object-oriented features (Fortran 2003),
meaning that Legolas requires relatively recent compilers. We build using [CMake](https://cmake.org).
To make use of the post-processing framework [Pylbo](post_processing) (which is highly recommended since we
use a dedicated file format to store information) you need Python 3.6 or higher.
The list below gives an indication of what is needed, note that lower versions _may_ work (but we haven't tested that).
- gfortran v8.x+
- CMake v3.12+
- Python v3.6+
- make
- BLAS and LAPACK v3.5+

For more information see the website.

## Compilation
For compilation instructions we refer to the website.

## Directory contents
Below is an overview of what's in the directory:
- `.github/workflows`: all `.yml` files used for automatic testing.
- `docs`: everything related to documentation, contains the Jekyll/Sphinx/Ford configuration files. The website is built from this folder.
- `post-processing`: contains the Python package Pylbo.
- `setup_tools`: scripts related to Legolas setup, also contains the `findBLAS` and `findLAPACK` CMake files.
- `src`: all legolas source files and subdirectories.
- `tests`: everything related to testing
    - `core_tests`: tests for the Legolas core subroutines, uses [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).
    - `pylbo_tests`: tests for the Pylbo framework, does checks on internal Pylbo routines + Legolas interface for multiruns.
    - `regression_tests`: tests for Legolas output, compares results from newly compiled code to expected (stored) output.
- `CMakeLists.txt`: used for build configuration
- `legolas_config.par`: example parameter file for legolas configuration.
- `pylbo_wrapper.py`: wrapper script to interface with Pylbo.
