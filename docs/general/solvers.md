---
title: Solvers
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "chevron-circle-down"
last_modified_at: 2022-07-20
---

Legolas has interfaces implemented to various BLAS, LAPACK and ARPACK routines.
Below is an overview of which routines you can call, which problems are supported and how
you can configure the parfile to select the solver you want.
We have have a general eigenvalue problem of the form

$$ A\mathbf{x} = \omega B\mathbf{x} $$

where $A$ is a non-symmetric and non-Hermitian complex matrix. The $B$-matrix is always real, and in most
cases also symmetric and positive definite (depending on the physics, the Hall electron inertia term for example
breaks positive definiteness). Both matrices are block-tridiagonal, meaning they are very sparse.

<i class="fas fa-lightbulb" aria-hidden="true"></i>
**Note**: A general strategy for a thorough investigation of a certain spectrum may be as follows: first a low-resolution QR-invert run is done,
which will reveal spectral regions of interest. This can then be followed-up by a higher-resolution run using QR-invert, and/or a shift-invert Arnoldi run
near the interesting regions. Comparing the eigenvalues between both solution strategies and at different resolutions will be a good indicator of their convergence.
{: .notice--success}


## QR-invert
This is the default solver that Legolas uses, which relies on an inversion of the B-matrix to write
the eigenvalue problem in the form

$$ B^{-1}A\mathbf{x} = \omega\mathbf{x}. $$

The LAPACK routine [`dgetrf`](https://netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html)
is used to calculate the LU factorisation of $B$, followed by a call to
[`dgetri`](https://netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga56d9c860ce4ce42ded7f914fdb0683ff.html)
which uses that factorisation to invert the $B$-matrix.
Finally a call to LAPACK's [`zgeev`](https://netlib.org/lapack/explore-html/db/d55/group__complex16_g_eeigen_ga0eb4e3d75621a1ce1685064db1ac58f0.html)
is made which returns all eigenvalues and optionally the right eigenvectors.

Note that we ensure that the $B$-matrix is always nicely conditioned, such that inversion does not yield problems.
This solver can be explicitly specified in the `solvelist` through
```fortran
&solvelist
  solver = "QR-invert"
/
```
and is called by default if no `solvelist` is supplied.
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
- High memory usage (for now, rework and optimisations are under development)
- Datfiles become very large at high resolution if eigenfunctions are included, this can be mitigated by saving a subset.
- May show numerical instability at very large $Re(\omega)$ values (far into the fast sequence), this should always be visually clear in the spectrum.
  In most cases this is not a problem, since those modes will not have their eigenfunctions resolved anyway.
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>


## QZ-direct
This is a variant of the QR-invert solver, with as main difference that the $B$-matrix is not inverted
such that the eigenvalue problem is kept in its general form.
The LAPACK routine [`zggev`](https://netlib.org/lapack/explore-html/db/d55/group__complex16_g_eeigen_ga79fcce20c617429ccf985e6f123a6171.html)
is used to solve the general eigenvalue problem, returning all eigenvalues.

This solver can be specified in the `solvelist` through
```fortran
&solvelist
  solver = "QZ-direct"
/
```

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
- High memory usage (cfr. QR-invert)
- Currently no support for eigenfunctions
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>

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
  arpack_mode = "general" | "shift-invert"
  number_of_eigenvalues = 20
  which_eigenvalues = "LM" | "SM" | "LR" | "SR" | "LI" | "SI"
  maxiter = 2500
  ncv = 50
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
- `ncv`: the number of basis vectors to use during the iteration, must be at least `number_of_eigenvalues + 1`.
  Defaults to `ncv = 2 * number_of_eigenvalues`.

If the iterative solver reaches `maxiter`, only a number $j < k$ eigenvalues will be converged.
Legolas will notify you how many are converged, and you can still plot these $j$ eigenvalues and their eigenfunctions.

Note that ARPACK is better at finding large eigenvalues. We recommend using the shift-invert mode
if you want better performance for smaller eigenvalues. Ideally a combination of both is used, where
one first solves for all eigenvalues using QR-invert or the standard/general Arnoldi solver, locate
spectral regions of interest, and then follow-up with shift-invert at those locations.

### Points of note
Before going over the implemented Arnoldi-routines a few important points are highlighted here. On this page only a few key points regarding ARPACK are explained,
for a thorough overview of all its possibilities we refer to [the documentation (pdf)](http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf).

**The computational time is strongly dependent on the choice of `ncv`, `number_of_eigenvalues`, `which_eigenvalues` (and `sigma` for shift-invert)**.
The piece of text below is quite important and is taken directly from the ARPACK documentation, section 2.3.3. The parameter `nev` refers to `number_of_eigenvalues`.
> For a given `ncv`, the computational work required is proportional to `n $\cdot$ ncv$^2$` FLOPS. Setting `nev` and `ncv` for optimal performance is very much
> problem dependent. If possible, it is best to avoid setting `nev` in a way that will split clusters of eigenvalues. For example, if the five smallest eigenvalues
> are positive and on the order of $10^{-4}$ and the sixth eigenvalue is on the order of $10^{-1}$ then it is probably better to ask for `nev = 5`
> than for `nev = 3` even if the three smallest are the only ones of interest.
>
> Setting the optimal value of `ncv` relative to `nev` is not completely understood. As with the choice of `which_eigenvalues`, it depends upon the underlying approximation
> properties of the Iterative Method as well as the distribution of the eigenvalues of the $A$-matrix. As a rule of thumb, `ncv $\geq$ 2 $\cdot$ nev` is reasonable.
> There are tradeoffs due to the cost of matrix-vector products, the cost of the implicit restart mechanism and the cost of maintaining the orthogonality of the Lanczos vectors.
> If the matrix-vector product is relatively cheap (as is the case for Legolas), then a smaller value of `ncv` may lead to more matrix-vector products, but an overall decrease
> in computation time.



### General mode
Set `arpack_mode = "general"`. This solves the eigenvalue problem $AX = \omega BX$  in its general form, which technically comes down to solving the problem $R = B^{-1}AX = \omega X$.
In order to exploit matrix sparsity we never explicitly calculate $R$, but do the following instead.

In a first step the $B$-matrix is converted from its linked-list representation
to a banded matrix datastructure. Then, at every step in the iteration:
1. The matrix-vector product $u = AX$  is calculated, in which the $A$-matrix is kept in its linked-list representation. This makes the matrix-vector product efficient and fast.
2. The linear system $BR = u$ is solved for R using LAPACK's banded matrix routines and passed back to the solver.

The above steps are repeated until convergence or until `maxiter` is reached.


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
- May require trial-and-error to find an optimal set of parameters
- Fast for the largest eigenvalues (`"LM", "LR", "LI"`), slower for the others
- Only calculates spectrum extremes (largest/smallest) which are usually not modes of interest.
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>

### Shift-invert mode
Running ARPACK in shift-invert mode is usually used to enhance convergence of certain spectral regions. It relies on a transformation of the eigenvalue problem to

$$
 \Bigl(A - \sigma B\Bigr)^{-1} BX = \nu X \qquad \text{where} \qquad \nu = \dfrac{1}{\omega - \sigma}
$$

This transformation is particularly useful for finding eigenvalues near $\sigma$. It transforms small eigenvalues into large ones, making them perfectly suitable for calculation with ARPACK. Here it is **recommended** to keep `which_eigenvalues = "LM"`, since the eigenvalues $\nu_j$ of $ C \equiv \left(A - \sigma B\right)^{-1}B$ that are largest in magnitude correspond to the values $\omega_j$ of the original eigenvalue problem that are nearest to the shift $\sigma$ in absolute value.

The value of $\sigma$ can be specified by adding it to the solvelist, like so
```fortran
&solvelist
  sigma = (1.0d0, 0.05d0)
/
```
and should be a complex tuple (standard Fortran notation for complex numbers).

The eigenvalue problem is tackled as follows.

First we calculate $A - \sigma B$ and convert that to banded form. Then, at every step in the iteration:
1. The matrix-vector product $u = BX$  is calculated, in which the $B$-matrix is kept in its linked-list representation. This makes the matrix-vector product efficient and fast.
2. The linear system $\left(A - \sigma B\right)R = u$ is solved for $R$ using LAPACK's banded matrix routines and passed back to the solver.

The above steps are repeated until convergence or until `maxiter` is reached. In the final step we retransform the eigenvalues $\nu_j$ to the
eigenvalues $\omega_j$ of the original problem using

$$
\omega_j = \sigma + \dfrac{1}{\nu_j}
$$


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
- May require trial-and-error to find an optimal set of parameters
- Runtime is strongly dependent on the chosen shift $\sigma$
{% endcapture %}
<div class="notice--danger">
  {{ cons | markdownify }}
</div>
