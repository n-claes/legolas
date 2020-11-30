---
title: Parameter options
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_label: "Namelists"
toc_icon: "code"
last_modified_at: 2020-09-02
---

The parameter file allows for a full customisation of all Legolas variables. Lots of these variables
have default values, except for the ones you _have_ to set. These are initialised using a Fortran `ieee_quiet_nan`,
propagating the value during setup. Legolas explicitly checks for nonsense values and will warn you if something
is not right. This not only includes checks on improperly set variables but also on errors when setting the 
equilibrium: if you're dividing by zero somewhere, for example, Legolas with catch it.

**Note:** Variables that have no default value must be set either through the parfile or in an equilibrium
submodule. For the `geometry`, `x_start` and `x_end` variables we recommend setting these in a submodule.
{: .notice--info}

You can create the parfile manually using this guide, or let Pylbo take care of generating one.
More information on generating parfiles with Pylbo can be found [here](../../pylbo/generate_parfiles).

## Gridlist
This namelist includes all grid-related variables.

**Warning:** For cylindrical geometries we override the start value with `r = 0.025` if it is initialised
with zero. Having `r = 0` does not give issues per se but can introduce spurious eigenvalues: some matrix elements 
divide by $r^2$, which is not zero but corresponds to the first point in the Gaussian grid (which is very small) 
and can hence yield huge values in the A-matrix. You can force `r = 0` through `force_r0 = .true.` in the gridlist, but this is
not recommended.
{: .notice--warning}

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |    :----        |          :---:      |
| geometry     | string | geometry of the setup, `"Cartesian"` or `"cylindrical"`       | `""` (must be set)   |
| x_start      | real | start of the grid       | **NaN** (must be set)      |
| x_end        | real | end of the grid         | **NaN** (must be set)      |
| gridpoints   | int | base number of gridpoints to use |  31 |
| force_r0     | logical | forces `r=0` in cylindrical geometry | `.false.`   |
| mesh_accumulation | logical | uses mesh accumulation based on two Gaussian integrals   |   `.false.` |
| ev_1         | real | expected value of the first Gaussian integral (only used with mesh accumulation)  | 1.25  |
| sigma_1      | real | sigma of the first Gaussian integral (only used with mesh accumulation)   | 1.25  |
| ev_2         | real | expected value of the second Gaussian integral (only used with mesh accumulation)  | 1.25  |
| sigma_2      | real | sigma of the second Gaussian integral (only used with mesh accumulation)   | 1.25  |

## equilibriumlist
This namelist includes all equilibrium-related variables.

**Note:** The parameter `boundary_type` only has one allowed value for the moment (`"wall"`), additional possibilities
will be added in future versions.
{: .notice--info}

**Warning:** The parameter `remove_spurious_eigenvalues` is _really_ not recommended and should only be used as a last
resort. Keeping `r = 0.025` in cylindrical geometry usually solves all related issues you may have.
{: .notice--warning}

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |   :----        |          :---:      |
| equilibrium_type  | string | the equilibrium configuration to use, see the [various options](../../general/equilibria) | `"adiabatic_homo"`    |
| boundary_type |  string |   boundary conditions to use  | `"wall"`  |
| use_defaults  | logical | if `.true.`, uses the default values in the equilibrium submodule | `.true.` |
| remove_spurious_eigenvalues   | logical | removes the outer `nb_spurious_eigenvalues` purely real eigenvalues from the spectrum     | `.false.` |
| nb_spurious_eigenvalues | int | number of purely real eigenvalues to remove on each side of the imaginary axis  | 1 |


## physicslist
This namelist includes all physics-related variables.

**Note:** If you use the default submodule values through `use_default = .true.`, variables like `radiative_cooling` 
and `resistivity` (and related ones) are set in the submodules themselves. It's good practice to set them in the 
submodule instead of through the parfile. The parfile functionality for these variables is provided for possible
parameter studies, to prevent having to recompile the submodule over and over again.
{: .notice--info}

| Parameter    | Type | Description      | Default value     |
| :---         |  :---: |  :----        |          :---:      |
| mhd_gamma | real | the value for the ratio of specific heats gamma   | 5/3   |
| flow  | logical | if `.true.` flow is included  | `.false.`     |
| radiative_cooling | logical | if `.true.` optically thin radiative losses are included | `.false.`     |
| ncool | int | amount of points used to interpolate the cooling curve  | 4000  |
| cooling_curve |  string | cooling curve to use, options are `"jc_corona"`, `"dalgarno"`, `"ml_solar"`, `"spex"`, `"spex_dalgarno"` and `"rosner"` | `"jc_corona"` |
| external_gravity  | logical | if `.true.` external gravity is included | `.false.` |
| thermal_conduction | logical | if `.true.` thermal conduction is included   | `.false.` |
| use_fixed_tc_para | logical | if `.true.`, uses a constant value for parallel conduction    | `.false.` |
| fixed_tc_para_value  | real |  constant value to use if `use_fixed_tc_para = .true.` | 0.0   |
| use_fixed_tc_perp | logical | if `.true.`, uses a constant value for perpendicular conduction    | `.false.` |
| fixed_tc_perp_value   | real | constant value to use if `use_fixed_tc_perp = .true.` | 0.0   |
| resistivity   | logical | if `.true.` resistivity is included   | `.false.` |
| use_fixed_resistivity | logical | if `.true.`, uses a constant value for the resistivity    | `.false.` |
| fixed_eta_value   | real | constant value to use if `use_fixed_resistivity = .true.` | 0.0   |
| use_eta_dropoff   | logical | if `.true.`, smoothly drops off the resistivity profile to zero near the edges using a hyperbolic tangent profile  | `.false.` |
| dropoff_edge_dist | real | distance between grid edge and the smoothened profile center  | 0.05  |
| dropoff_width | real | sets the width over which the profile is smoothened out   | 0.1   |

## solvelist
This namelist includes all solver-related variables. For more information, see [Solvers](../../general/solvers)

| Parameter     | Type  | Description   | Default value     |
| :---          | :---: | :----         | :---:             |
| solver        | string    | which solver to use | `"QR-invert"` |
| arpack_mode   | string    | the mode for ARPACK, only used if `solver="arnoldi"` | `"standard"`  |
| number_of_eigenvalues | int   | number of eigenvalues to calculate, only used if `solver="arnoldi"`   | 10    |
| which eigenvalues | string    | which eigenvalues to calculate, only used if `solver="arnoldi"`   |   `"LM"`  |
| maxiter       | int   | the maximum number of iterations the Arnoldi solver may take, only used if `solver="arnoldi"` | 10N   |
| sigma         | complex   | sigma-value around which to do shift-invert, only for `arpack_mode="shift-invert"`  | **NaN** |

## unitslist
This namelist includes all units-related variables.

**Note:** You can set the unit normalisations in the submodule as well. You can either supply
`unit_density` or `unit_temperature`, but not both. See [Units](../../general/equations_physics/#units)
for more information.
{: .notice--info}

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |  :----        |          :---:      |
| cgs_units     | logical | if `.true.` cgs normalisations are used   | `.true.`  |
| unit_density  | real | sets the unit density | N/A (we use temperature instead) |
| unit_temperature  | real | sets the unit temperature | $10^6$ K |
| unit_magneticfield |  real | sets the unit magnetic field    | $10$ Gauss    |
| unit_length   | real | sets the unit length  | $ 10^9$ cm    |

## savelist
This namelist includes all variables related to data output.

{% capture notice-info %}
**Note:** Legolas uses a logging-based system to write information to the console, similar to Python's logging module.
The level of output is controlled through the integer `logging_level`:
* **level 0** : only critical errors are printed, everything else is suppressed
* **level 1** : only critical errors and warnings are printed
* **level 2** : prints info messages, warnings and critical errors
* **level 3** : (or higher) prints debug messages in addition to all of the above
{% endcapture %}
<div class="notice--info">
  {{ notice-info | markdownify }}
</div>

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |   :----        |          :---:      |
| write_matrices    | logical | if `.true.` the matrices are written to the datfile   |   `.false.`   |
| write_eigenfunctions  | logical | if `.true.` eigenfunctions are calculated and written to the datfile  |   `.true.`  |
| show_results  | logical | calls the python wrapper after the run is finished and plots the results, requires `pylbo_wrapper.py` to be in the same directory as the executable.  | `.true.`  |
| basename_datfile  | string | the base name of the datfile, this is prepended with the output directory and appended with `".dat"`  | `"datfile"`   |
| basename_logfile  | string | if not an empty string, saves the eigenvalues to a plain text file with name `basename_logfile.log` (mainly used for testing)  | `""`  |
| output_folder | string | the output folder, this is prepended to the datfile name  | `"output"`    |
| logging_level | int | controls the amount of console info, see the note above   | 2 |
| dry_run   |  logical | if `.true.` forces all eigenvalues to zero and does not call the solver, useful for iterating over the equilibrium implementation    | `.false.` |


## paramlist
This namelist handles constant parameters which are used in the equilibrium prescriptions, and is only read
in if `use_defaults = .false.`. This is meant to be used in parameter studies, for a comprehensive list of all
the possible parameters we refer to [this page](../../src-docs/ford/module/mod_equilibrium_params.html) in the
Legolas source code documentation. To see which parameters are used in which equilibria, take a look
[here](../equilibria/).

