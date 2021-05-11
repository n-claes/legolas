---
title: Regression tests
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2021-02-03
---

The "main" testing framework for the Legolas code handles regression tests, which basically compare code
output from a recent commit to previously (known) results.
These "answers" are stored in `tests/regression_tests/answers` as datfiles and logfiles, and we try to include at least
one regression test per implemented equilibrium.
The regression tests are pytest-based, meaning pytest will use Pylbo to interface with Legolas and run the
executable. Most of them run with 51 gridpoints and no eigenfunctions (although in some cases these are included), so
running a single test should take no more than 10 seconds.

## Comparison between new and stored results
Once the case is run and output files are present, pytest will proceed to do a comparison between the two.
However this is not straightforward, since the BLAS and LAPACK routines on different platforms (i.e. macOS or Linux)
may give _slightly_ different results, on top of cross-platform numerical deviations. For example, the tests run
on Linux can pass without issues, but may fail on macOS because the answers were generated on a Linux platform.

This means that we can not do a straight comparison between the eigenvalues. Including a "simple" tolerance is also
not possible, as explained below:
- There is no guarantee that all eigenvalues appear in the same order, and they hence have to be sorted. We therefore
  sort them based on the real and imaginary parts, in that order. However, an eigenvalue that is purely real may have
  a very small (e.g. $10^{-12}$) imaginary part due to numerical errors, while it is considered zero
  (or negative that small value) in the stored answers. This messes up the sorting, and if this happens a few times
  in the sequence you start comparing eigenvalues that shouldn't be compared in the first place, and the test fails.
- Say the lists are correctly sorted, then the second issue that arises is the actual comparison. For "large" eigenvalues,
  e.g. $250.1234 + 180.23i$, an error of $\pm 0.1$ (or even one) for each component is perfectly fine. However, for
  an eigenvalue of e.g. $0.00345 - 0.000152i$ even an error of $\pm 10^{-3}$ is way too big. Switching to a relative
  comparison is also not really an option, since deviations strongly depend on which equilibrium we are testing and
  you never know beforehand how "big" the relative error may be in order to be robust to small changes.

What the regression tests actually do is a figure-based _pixel comparison_. Both the eigenvalues of the tests and answers
are plotted in the imaginary plane and their values are transformed to $xy$ pixel coordinates using matplotlib's
[`transData.transform`](https://matplotlib.org/3.1.1/api/transformations.html#module-matplotlib.transforms).
We then look at a radius of 1 pixel around every point of the test case, and if there is a point from the answer
tests in this radius we flag it as fine and move on to the next point.
We've found that this method is actually quite robust, since a spectrum is quite sensitive to
even small changes. For example, modifying the value of a constant in the equations by a small number (say 0.1%)
is in most cases enough to fail the tests, since there will be at least one point that is sufficiently shifted
such that it lies outside the 1 pixel circle.

Points that are flagged as failed are added to the list, and in the end they will be plotted together with the base answers
for a visual comparison. A .png file for every failed test will be present in the `testing/regression_tests/results` folder,
with `FAILED_` prepended to the name.

## Locally running regression tests
Similar as to the [pylbo tests](../test_pylbo) you'll need pytest for this with the same plugins.
Navigate to the `tests/regression_tests` folder and execute
```bash
pytest regression.py test* -v --mpl --mpl-results-path=results
```
Due to the structure of the testing hierarchy we want to run the file `regression.py` first to properly configure
each setup, followed by all other tests in random order. The `-v` argument enables more clear testing output,
while `--mpl --mpl-results-path` enables image baseline comparison for the various multiruns. 

Regression tests are run automatically for each commit and pull request to the `master` and `develop` branches.
