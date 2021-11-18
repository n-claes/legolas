---
title: The parfile
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_label: "Namelists"
toc_icon: "code"
last_modified_at: 2021-07-26
---

The Legolas parfile allows for a full customisation of all Legolas variables. Lots of these variables
have default values, except for the ones you _have_ to set. These are initialised using a Fortran `ieee_quiet_nan`,
so if a variable is not properly set a `NaN` is propagated and flagged during the sanity checks.
Legolas explicitly inspects your equilibrium for nonsensical values and will warn you if something is not right.
Whenever the code catches something during its inspection phase it either logs a warning and continues running, or
exits with an error depending on the severity of the issue.

You can create the parfile manually using this guide, or let Pylbo take care of generating one.
More information on generating parfiles with Pylbo can be found [here](../../pylbo/generate_parfiles).

## gridlist
This namelist includes all grid-related variables.

<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
**Warning:** For cylindrical geometries we override the start value with `r = 0.025` if it is initialised
with zero. Some matrix elements scale with $1/r$ and may blow up when `r = 0`, resulting in numerical difficulties. Note that the matrix elements and equilibria
are evaluated in the _Gaussian_ grid, the first point of which may be very small but will never be exactly zero. <br><br>
Setting `force_r0 = .true.` in the gridlist forces the `r = 0` condition, but this is not recommended.
{: .notice--warning}

| Parameter    | Type     | Description | Default value     |
| :---         | :---:    |    :----  |          :---:      |
| geometry     | string   | geometry of the setup, `"Cartesian"` or `"cylindrical"`. Must be set in the parfile or equilibrium module   | - |
| x_start      | real     | start of the grid, must be set  |  - |
| x_end        | real     | end of the grid, must be set   | - |
| gridpoints   | int      | base number of gridpoints to use  |  31 |
| force_r0     | logical  | forces `r=0` in cylindrical geometry | `.false.` |
| coaxial      | logical  | use a coaxial inner boundary in cylindrical geometry | `.false.` |

## equilibriumlist
This namelist includes all equilibrium-related variables.

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |   :----        |          :---:      |
| equilibrium_type  | string | the equilibrium configuration to use, see the [various options](../../general/equilibria) | `"adiabatic_homo"`    |
| boundary_type |  string |   boundary conditions to use  | `"wall"`  |
| use_defaults  | logical | if `.true.`, uses the default values in the equilibrium submodule | `.true.` |


## physicslist
This namelist includes all physics-related variables.

| Parameter             | Type    | Description      | Default value     |
| :---                  |  :---:  |  :----        |          :---:      |
| mhd_gamma             | real    | the value for the ratio of specific heats $\gamma$   | $\frac{5}{3}$   |
| flow                  | logical | inclusion of background flow effects  | `.false.`     |
| radiative_cooling     | logical | inclusion of optically thin radiative losses | `.false.`     |
| ncool                 | int     | amount of points used to interpolate the radiative cooling tables  | 4000  |
| cooling_curve         | string  | cooling curve to use, options are `"jc_corona"`, `"dalgarno"`, `"ml_solar"`, `"spex"`, `"spex_dalgarno"` and `"rosner"` | `"jc_corona"` |
| external_gravity      | logical | inclusion of external gravity | `.false.` |
| thermal_conduction    | logical | inclusion of thermal conduction  | `.false.` |
| use_fixed_tc_para     | logical | whether to use a constant value for thermal conduction parallel to the magnetic field lines    | `.false.` |
| fixed_tc_para_value   | real    | constant value to use if `use_fixed_tc_para = .true.` | 0   |
| use_fixed_tc_perp     | logical | whether to use a constant value for thermal conduction perpendicular to the magnetic field lines    | `.false.` |
| fixed_tc_perp_value   | real    | constant value to use if `use_fixed_tc_perp = .true.` | 0   |
| resistivity           | logical | inclusion of resistivity  | `.false.` |
| use_fixed_resistivity | logical | whether to use a constant resistivity value    | `.false.` |
| fixed_eta_value       | real    | constant value to use if `use_fixed_resistivity = .true.` | 0   |
| use_eta_dropoff       | logical | if `.true.`, smoothly drops off the resistivity profile to zero near the edges using a hyperbolic tangent profile  | `.false.` |
| dropoff_edge_dist     | real    | distance between grid edge and the smoothened profile center if `use_eta_dropoff = .true.`  | 0.05  |
| dropoff_width         | real    | sets the width over which the profile is smoothened out if `use_eta_dropoff = .true.`   | 0.1   |
| viscosity             | logical | inclusion of viscosity | `.false.` |
| viscous_heating       | logical | whether to include viscous heating in the viscosity prescription | `.false.` |
| viscosity_value       | real    | constant value for the viscosity | 0 |
| incompressible        | logical | if `.true.`, uses an incompressible approximation | `.false.` |
| hall_mhd              | logical | inclusion of Hall effects | `.false.` |
| elec_intertia         | logical | whether to include the electron intertia term in Ohm's law if `hall_mhd = .true.` | `.false.` |

## solvelist
This namelist includes all solver-related variables. For more information, see [Solvers](../../general/solvers)

| Parameter     | Type  | Description   | Default value     |
| :---          | :---: | :----         | :---:             |
| solver        | string    | which solver to use | `"QR-invert"` |
| arpack_mode   | string    | the mode for ARPACK, only used if `solver="arnoldi"` | `"standard"`  |
| number_of_eigenvalues | int   | number of eigenvalues to calculate, only used if `solver="arnoldi"`   | 10    |
| which eigenvalues | string    | which eigenvalues to calculate, only used if `solver="arnoldi"`   |   `"LM"`  |
| maxiter       | int   | the maximum number of iterations the Arnoldi solver may take, only if `solver="arnoldi"` | 10N   |
| sigma         | complex   | sigma-value around which to do shift-invert, only for `arpack_mode="shift-invert"`  | - |

## unitslist
This namelist includes all units-related variables.

| Parameter    | Type | Description      | Default value     |
| :---         | :---: |  :----        |          :---:      |
| cgs_units     | logical | if `.true.` cgs normalisations are used   | `.true.`  |
| unit_density  | real | sets the unit density | N/A (we use temperature instead) |
| unit_temperature  | real | sets the unit temperature | $10^6$ K |
| unit_magneticfield |  real | sets the unit magnetic field    | $10$ Gauss    |
| unit_length   | real | sets the unit length  | $ 10^9$ cm    |
| mean_molecular_weight | real | the mean molecular weight to use | 1 |


## savelist
This namelist includes all variables related to data output.

{% capture notice-info %}
**Note:** Legolas uses a logging-based system to write information to the console, similar to Python's logging module.
The level of output is controlled through the integer `logging_level`:
* **level 0** : only critical errors are printed, everything else is suppressed. Use this for multiruns
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
| write_derived_eigenfunctions | logical | if `.true.` also calculates derived eigenfunction quantities ($\nabla \cdot \mathbf{B_1}$, $S$, $v_\parallel$, etc.) |  `.false.` |
| write_eigenfunction_subset | logical | if `.true.` only saves a part of the eigenfunctions to the datfile, based on a given radius and complex point | `.false.` |
| eigenfunction_subset_center | complex | point in the complex plane that defines the center of the circle in which to save eigenfunctions, needs `eigenfunction_subset_radius` to be set | - |
| eigenfunction_subset_radius | float | radius around `eigenfunction_subset_center`, all eigenvalues within this circle will have their eigenfunctions saved | - |
| show_results  | logical | calls the python wrapper after the run is finished and plots the results, requires `pylbo_wrapper.py` to be in the same directory as the executable.  | `.true.`  |
| basename_datfile  | string | the base name of the datfile, this is prepended with the output directory and appended with `".dat"`  | `"datfile"`   |
| basename_logfile  | string | if not an empty string, saves the eigenvalues to a plain text file with name `basename_logfile.log` (mainly used for testing)  | `""`  |
| output_folder | string | the output folder, this is prepended to the datfile name  | `"output"`    |
| logging_level | int | controls the amount of console info, see the note above   | 2 |
| dry_run   |  logical | if `.true.` forces all eigenvalues to zero and does not call the solver, useful for iterating over the equilibrium implementation    | `.false.` |


## paramlist
This namelist handles constant parameters which are used in the equilibrium prescriptions, and is only read
in if `use_defaults = .false.`. This is meant to be used in parameter studies, for a comprehensive list of all
the possible parameters we refer to [this page](../../ford/module/mod_equilibrium_params.html) in the
Legolas source code documentation. To see which parameters are used in which equilibria, take a look
[here](../equilibria/).

<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
**Note**: setting a variable with the paramlist does **NOT** automatically include it in the equilibrium configuration. For this to have any effect you have to _use_ that
variable in your submodule as well. For example, if you never set a magnetic field in your submodule and never use the `cte_B02` variable,
setting `cte_B02 = 1.0d0` in the paramlist will have no effect. Take a look at the different equilibria to see which parameters are supported.
{: .notice--warning}
