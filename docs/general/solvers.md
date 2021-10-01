---
title: Solvers
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "chevron-circle-down"
last_modified_at: 2021-07-26
---

Legolas has interfaces implemented to various BLAS, LAPACK and ARPACK routines.
Below is an overview of which routines you can call, which problems are supported and how
you can configure the parfile to select the solver you want.
Note that in (most) cases we have a general eigenvalue problem of the form

$$ A\mathbf{x} = \omega B\mathbf{x} $$

where $A$ is a non-symmetric and non-Hermitian complex matrix. The $B$-matrix is
symmetric, real and positive definite. Both matrices are block-tridiagonal, meaning
they are very sparse.


## QR-invert
{% capture pros %}
**Pros:**
- Fast
- Calculates complete spectrum and eigenfunctions
{% endcapture %}
<div class="notice--success">
  {{ pros | markdownify }}
</div>

{% capture cons %}
**Cons:**
- Higher resolution needed to resolve some modes
- Datfiles become very large at high resolution if eigenfunctions are included
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>

This is the default solver that Legolas uses, which relies on an inversion of the B-matrix to write
the eigenvalue problem in the form

$$ B^{-1}A\mathbf{x} = \omega\mathbf{x}. $$

The LAPACK routine `dgetrf` is used to calculate the LU factorisation of B, followed by a call
to `dgetri` which uses that factorisation to invert the B-matrix.
Finally a call to LAPACK's `zgeev` is made which returns all eigenvalues and optionally
the right eigenvectors.

Note that the B-matrix is always nicely conditioned, such that inverting it does not yield problems.
This solver can be explicitly specified in the `solvelist` through
```fortran
&solvelist
  solver = "QR-invert"
/
```
and is called by default if no `solvelist` is supplied.

## QZ-direct
{% capture pros %}
**Pros:**
- No inversion of the B-matrix needed
- Calculates complete spectrum
{% endcapture %}
<div class="notice--success">
  {{ pros | markdownify }}
</div>

{% capture cons %}
**Cons:**
- Noticably slower than QR-invert, especially for large matrices
- Higher resolution needed to resolve some modes
- Returns _generalised_ eigenvectors instead of standard ones
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>

This is a variant of the QR-invert solver, with as main difference that the B-matrix is not inverted
such that the eigenvalue problem is kept in its general form.
The LAPACK routine `zggev` is used to solve the general eigenvalue problem, returning all
eigenvalues and (optionally) the _generalised_ right eigenvectors.

This solver can be specified in the `solvelist` through
```fortran
&solvelist
  solver = "QZ-direct"
/
```

## ARPACK Reverse Communication Interface
Legolas has various solvers implemented which interface with the ARPACK package to
solve the eigenvalue problem. ARPACK is a reverse communication interface specifically designed to
solve large-scale, sparse matrix eigenvalue problems, and is hence perfectly suited for Legolas.
ARPACK can run in various modes, most notably a shift-invert method to probe
various parts of the spectrum, only returning eigenvalues of regions you are interested in.

The main difference with the LAPACK solvers is that one can query for only a number of eigenvalues
instead of the full spectrum. This is essentially the Fortran analog of SciPy's
[`eigs`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html)
method in Python, which is a wrapper to ARPACK in itself.

When using Arnoldi-based solvers the solvelist can be set as follows:
```fortran
&solvelist
  solver = "arnoldi"
  arpack_mode = "standard" | "general" | "shift-invert"
  number_of_eigenvalues = 100
  which_eigenvalues = "LM" | "SM" | "LR" | "SR" | "LI" | "SI"
  maxiter = 2500
/
```

- `arpack_mode`: which mode to use, see the subsections below.
- `number_of_eigenvalues`: this is an integer denoting the $k$ eigenvalues to be computed.
   Note that this number (obviously) has to be positive and should be smaller than the dimension of
   the eigenvalue problem (that is, `matrix_gridpts`).
- `which_eigenvalues`: denotes the type of eigenvalues that should be calculated
   - `"LM"`: eigenvalues with largest magnitude
   - `"SM"`: eigenvalues with smallest magnitude
   - `"LR"`: eigenvalues with largest real part
   - `"SR"`: eigenvalues with smallest real part
   - `"LI"`: eigenvalues with largest imaginary part
   - `"SI"`: eigenvalues with smallest imaginary part
- `maxiter`: integer which limits the maximum iterations the Arnoldi solver may take when
   doing the reverse communication. This defaults to 10 times the size of the eigenvalue problem,
   so for 100 gridpoints `maxiter` will be 10 x 100 x 16 = 16000, which is usually more than sufficient.
   However, sometimes (especially for small eigenvalues) this may not be enough,
   in which case you can increase this number.

If the iterative solver reaches `maxiter`, only a number $j < k$ eigenvalues will be converged.
Legolas will notify you how many are converged, and you can still plot these $j$ eigenvalues and their eigenfunctions.

Note that ARPACK is better at finding large eigenvalues. We recommend using the shift-invert mode
if you want better performance for smaller eigenvalues. Ideally a combination of both is used, where
one first solves for all eigenvalues using QR-invert or the standard/general Arnoldi solver, locate
spectral regions of interest, and then follow-up with shift-invert at those locations.

### Standard / general mode
{% capture pros %}
**Pros:**
- Calculates only specific eigenvalues and eigenvectors
- Reduced datfile size, only eigenvectors for requested eigenvalues are calculated
{% endcapture %}
<div class="notice--success">
  {{ pros | markdownify }}
</div>

{% capture cons %}
**Cons:**
- Fast for the largest eigenvalues (`"LM", "LR", "LI"`), significantly slower for the smaller ones
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>

_Standard mode_: set `arpack_mode = "standard"`.
This is analogeous to the QR algorithm, inverts the B-matrix and passes the eigenvalue
problem in standard form to the iterative solver.

_General mode_: set `arpack_mode = "general"`.
Solves the eigenvalue problem in its general form, however an inversion of the
B-matrix is still needed. The main difference with the standard mode is that $B\textbf{x}$ is also
used during the iteration.

### Shift-invert
{% capture pros %}
**Pros:**
- Ability to probe specific parts of the spectrum by shifting $\sigma$
- Better performance for small eigenvalues
{% endcapture %}
<div class="notice--success">
  {{ pros | markdownify }}
</div>

{% capture cons %}
**Cons:**
- May take a very long time if the shift is badly chosen
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>
Running ARPACK in shift-invert mode allows one to set a certain shift $\sigma$ and calculate
the shifted eigenvalues. Note that for this mode, the setting `which_eigenvalues` in the parfile
refers to the shifted eigenvalues

$$ \dfrac{1}{\omega_i - \sigma} $$

The value of $\sigma$ can be specified by adding it to the solvelist, like so
```fortran
&solvelist
  sigma = (1.0d0, 0.05d0)
/
```
and should be a complex tuple (standard Fortran notation for complex numbers).
Here we _really_ recommend that you first look at the spectrum using the QR solver in order to provide a "proper" guess for $\sigma$ and
to make sure that there are enough eigenvalues near the requested shift. The solver may take quite a long time if you request 100 eigenvalues for example,
and there are only 50 eigenvalues in the vicinity of the chosen shift.

Also note that we need the operator $[A - \sigma B]^{-1}B$, which implies that $\sigma$ can not be zero
in our case, because the matrix A can be singular (no magnetic field, for example) which removes the
guarantee that the system is properly conditioned.
