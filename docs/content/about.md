---
title: About
layout: single
last_modified_at: 2020-08-27
toc: true
last_modified_at: 2020-10-27
---

Legolas, which stands for **L**arge **E**igensystem **G**enerator for **O**ne-dimensional p**LAS**mas,
is a novel finite element code specifically developed to do MHD spectroscopy of 1D Cartesian or cylindrical equilibria.
These setups balance pressure gradients, Lorentz forces, centrifugal effects and gravity, and are enriched with
a multitude of physical effects ranging from optically thin radiative losses with modern cooling curves
to anisotropic thermal conduction and resistivity.

## Aims
The main goal of Legolas is to provide a new and modern user-friendly code, making use of the most recent developments
in linear algebra and theoretical (linear) MHD. This allows for new and detailed investigations of the full (M)HD
linear spectrum through fully realistic setups and parametric studies.

## Features
Legolas is written in modern Fortran and makes extensive use of recent features.
To that extent the code is highly modularised, with all pre-implemented equilibria contained
in dedicated submodules. This allows for easy maintenance, and makes the code ready to be extended with additional
physics and modern algorithmic requirements (mesh refinement, for example).

## Development
The Legolas code is being developed at the Centre for mathematical Plasma-Astrophysics
[(CmPA)](https://wis.kuleuven.be/CmPA), KU Leuven, Belgium. The current development team consists of _Niels Claes_ and
_Jordi De Jonghe_, two PhD students, under the supervision of prof. Rony Keppens.
