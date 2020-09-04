---
title: Testing the core routines
layout: single
sidebar:
  nav: "leftcontents"
---

Legolas has multiple _test suites_ that are used to test the code and make sure that new features or changes that are
introduced do not break previously working pieces of code. 
To that end we have automated core tests that use [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) to
perform unit tests of the core subroutines of Legolas itself.

To run these locally you need to have pFUnit installed, which has some 
[prerequisites](https://github.com/Goddard-Fortran-Ecosystem/pFUnit#prerequisites). Make sure the legolas code is
compiled, and then navigate to the `tests/core_tests` folder and run
```bash 
make
./test_legolas
```

These tests are run automatically for each commit and pull request to the `master` and `develop` branches.