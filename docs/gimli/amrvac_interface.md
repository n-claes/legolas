---
title: Interfacing with MPI-AMRVAC
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: false
last_modified_at: 2026-02-04
---


## Initialising MPI-AMRVAC with Legolas data
<!-- [Legolas file generation](../../gimli/file_generation) -->





## Numerical configurations in Legolas
Until now, all configurations in Legolas had to be defined analytically. However, the code also features a way to analyze one-dimensional numerical data. To achieve this, GIMLI is used to convert a set of numerical arrays to a Legolas-interpretable file. Legolas then imports the data from this file, interpolates to a uniform grid of user-specified resolution, and calculates the derivatives numerically.

{% capture note %}
<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
**Note:** numerical derivatives are unreliable near sharp transitions, so it is highly recommended to inspect the numerical derivatives before interpreting results.
{% endcapture %}
<div class="notice--warning">
  {{ note | markdownify }}
</div>

As a toy example, consider the Harris current sheet once more. Using NumPy we define it numerically:
```python
import numpy as np

N = 1000
x = np.linspace(-10, 10, N)
temperature = np.ones(N)
magnetic_2 = np.tanh(x)
density = np.ones(N) - magnetic_2**2 / (2. * temperature)
```
This is then passed to the GIMLI class `NumericalEquilibrium` as a dictionary,
```python
import pylbo.gimli as gl

arrays = {
    "x" : x,
    "rho0" : density,
    "T0" : temperature,
    "B02" : magnetic_2
}
equilibrium = gl.NumericalEquilibrium(arrays)
```
The accepted dictionary keys are `x`, `r`, `u1`, `rho0`, `T0`, `v01`, `v02`, `v03`, `B01`, `B02`, `B03`, and `grav`, where exactly one of the first three (`x`, `r`, `u1`) should be present, and `grav` is used for gravity. All arrays should be of equal length, and missing keys are assigned a zero array of that length. To generate the data file for Legolas, call
```python
equilibrium.to_legolas_arrays()
```
to create a `.lar` file (**L**egolas **ar**rays) in the current directory. Alternatively, the location can be specified with the `loc` argument, and `filename` is used to set the name of the file.

To run Legolas on the numerical configuration, set the equilibrium in the parfile to `equilibrium_type = "numerical"`, and specify in `paramlist` both `input_file` (default `input_file = "arrays.lar"`) and `n_input`, the resolution for the interpolation. The predefined numerical equilibrium (`smod_equil_numerical.f08`) does not include `v01`, `B01`, and `grav`, but is easily modified if necessary.

**Note:** a sixth-order accurate central difference stencil is used to calculate the numerical derivatives. Near the edges a sixth-order accurate forward and backward difference stencil is used for the left and right boundary, respectively.{: .notice--info}

{% capture note %}
<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
**Note:** the calculation of the numerical derivatives assumes an equally spaced grid. Hence, this method cannot be combined with `set_custom_grid` or `set_spacing_function`.
{% endcapture %}
<div class="notice--warning">
  {{ note | markdownify }}
</div>