---
title: Testing the core routines
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-10-27
---

Legolas has multiple _test suites_ that are used to test the code and make sure that new features or changes that are
introduced do not break previously working pieces of code.
To that end we have automated core tests that use [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to
perform unit tests of the core subroutines of Legolas itself.

To run these locally you need to have pFUnit installed, which has some
[prerequisites](https://github.com/Goddard-Fortran-Ecosystem/pFUnit#prerequisites). Make sure the legolas code is
compiled, and then navigate to the `tests/core_tests` folder and run
```bash
mkdir build
cd build
cmake ..
make
```
to compile them. To speed things up you can always compile using `make -j 8` (or any amount of cpus you can spare).
The executable to run the tests will be placed in the `core_tests` folder, run it through
```bash
./test_legolas
```

These tests are run automatically for each commit and pull request to the `master` and `develop` branches.
