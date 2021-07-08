---
title: Legolas history & changelog
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-11-30
---

### [release 1.0.0](https://github.com/n-claes/legolas/releases/tag/v1.0.0)
**What's new:**
- CMake can now also link to ARPACK, SCALAPACK and MUMPS
- Complete overhaul of the solver modules:
  - Dedicated submodules for each solver
  - Implemented a QZ-direct solver, which solves Ax = wBx in its general form
- ARPACK-based solvers are implemented:
  - Solves the standard or general eigenvalue problem for selected eigenvalues
  - Includes a shift-invert method to probe specific spectral regions

### [version 0.9.0](https://github.com/n-claes/legolas/releases/tag/v0.9.0)
- First public beta release
