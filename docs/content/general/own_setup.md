---
title: Running your own setup
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "tools"
toc_label: "User submodule configuration"
last_modified_at: 2020-09-03
---

This page explains how to configure the custom user submodule to your needs and run your own setup.
To start navigate to any directory of your choosing and call the `setuplegolas.py` setup script.
See the [installation guide](../../getting-started/installation/#2-out-of-repository-build-recommended)
for more information.

# The parfile
When running the setup script and no parfile is found you are asked if you want to copy over a default template.
You'll need this file for the basic Legolas configuration, so please do so. Take a look at the different
[parameter options](../parameter_file/) as well. Configure this file to your needs.

# The smod_user_defined.f08 submodule
In this module you will set up your equilibrium configuration. You can either choose to copy over the default template
from the legolas `src` directory, or create your own. Both are fine, as long as the file is named `smod_user_defined.f08`
and the general template looks like this:
```fortran
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  module subroutine user_defined_eq()

    ! your setup

  end subroutine user_defined_eq
end submodule smod_user_defined
```
This is all that's needed. The `mod_equilibrium` statement between brackets refers to the parent module,
which contains the interface of the subroutines that are implemented in their respective submodules.
This also implies that you have access to _all_ global `use` statements in the `mod_equilibrium` module, without
having to explicitly import them.

In what follows all code samples go _inside_ the `user_defined_eq()` subroutine, since this is the only routine you
need to configure yourself.

## Setting the grid
You have several options here. Either you hardcode the grid parameters yourself, setting `geometry`, `x_start` and `x_end`
to fixed values. However, in some use cases it may be necessary to vary grid parameters as part of parametric studies.
We therefore provided the `allow_geometry_override` subroutine, which can be used to set default values that can be
overridden through the parfile.

**Note:** The variables `geometry`, `x_start`, `x_end` and the subroutine `call initialise_grid()` are all known
in the parent module and don't have to be imported explicitly.
{: .notice--info}

_Option 1_: full override
```fortran
call allow_geometry_override(default_geometry='cylindrical', &
                             default_x_start=0.0d0, &
                             default_x_end=1.0d0)
```
_Option 2_: partial override
```fortran
x_end = 1.0d0
call allow_geometry_override(default_geometry="Cartesian", default_x_start=0.0d0)
```
_Option 3_: hardcoded
```fortran
geometry = "Cartesian"
x_start = 0.0d0
x_end = 1.0d0
```
In the first option above you can switch between `"Cartesian"` and `cylindrical` geometries and custom grid boundaries
solely using the parfile, while in option 2 you can only modify the outer grid boundary. Option 3 does not allow for
an override at all.
After these three variables are set, you tell the code it should initialise the grid through
```fortran
call initialise_grid()
```
which allocates the arrays and initialises the grid, geometry scale factor and Gaussian grid.

## Using variables
Legolas uses a dedicated module which contains a multitude of variables that you can use in your setups.
For a comprehensive list we refer to [this page](../../src-docs/ford/module/mod_equilibrium_params.html#variable-k2).
You can use all of these in your setups, but note that these should be explicitly imported (except for `k2` and `k3`).
An example is given below.

**Note:** The variables `use_defaults`, `k2` and `k3` are known in the parent module. All other parameters you want
to use should be explicitly imported.
{: .notice--info}

```fortran
use mod_equilibrium_params, only: alpha, beta, p1

if (use_defaults) then
  k2 = 0.0d0
  k3 = 1.0d0
  alpha = 1.0d0
  beta = 5.0d0
  p1 = 1.0d0 / 2.0d0
end if
```
This example uses the variables `alpha`, `beta` and `p1` from the list and wraps them in the `use_defaults` condition.
This means that Legolas will use those values as a default if the `paramlist` namelist is not present in the parfile.
If `use_defaults` is set to `.false.` you _have_ to specify them in the `paramlist`, where you choose values for these
variables. If you forget to set a variable that you use later on its value will be `NaN`, meaning the code will catch it
during the pre-run checks and warn you.

## Setting units
**Note:** All public variables/subroutine defined in
[`mod_units`](../../src-docs/ford/module/mod_units.html#variable-unit_length) are known to the parent module and
do not have to be imported explicitly in the submodule.
{: .notice--info}

Specifying the unit normalisations can either be done through the parfile (the `unitslist`) or in the submodule
directly. This is only relevant if you include radiative cooling, thermal conduction or (temperature-dependent)
resistivity, since these units are used to re-dimensionalise the variables for the corresponding calculations and/or
table lookups. To set the units you can do, for example
```fortran
use_cgs = .true.
call set_normalisations(new_unit_density=1.6727d-15, &
                        new_unit_magneticfield=22.5d0, &
                        new_unit_length=1.0d10)
```
which sets values in cgs units, if you set `use_cgs = .false.` these are interpreted as SI units. Instead of specifying
`new_unit_density` you can choose a reference temperature `new_unit_temperature` instead,
see [units](../equations_physics/#units) for more information.

## Including additional physics
**Note:** You can set everything here either in the parfile or in the submodule. In case of the latter,
the logicals `resistivity`, `flow`, `radiative_cooling`, `external_gravity` and `thermal_conduction` are
known in the parent module. Variables like `use_fixed_resistivity`, `fixed_eta_value` and related ones on the other hand
must be imported explicitly.
{: .notice--info}

Including additional physics in your setup is quite straightforward: you set the corresponding logical to `.true.` and
set up the additional variables.

### Resistivity
In the case of resistivity you have two options after setting `resistivity = .true.`:
- A constant value over the entire grid (with a possible drop-off near the edges,
  see the [physicslist](../parameter_file/#physicslist) namelist). Set `use_fixed_resistivity = .true.` and
  specify `fixed_eta_value`.
- A realistic temperature-dependent resistivity profile based on the Spitzer resistivity. Here it suffices
  to only specify `resistivity = .true.`, since `use_fixed_resistivity` is `.false.` by default.
  It's perhaps a good idea to specify unit normalisations in this case, if the default ones are not sufficient to your
  needs.

### Thermal conduction
Set `thermal_conduction = .true.`. This is analogous to resistivity, in the sense that you can either specify a
constant value for the parallel and perpendicular components (separately), or use a fully realistic profile that
depends on density, temperature and magnetic field. The constant values are controlled through the variables
`use_fixed_tc_para` and `fixed_tc_para_value` for the parallel thermal conduction component, and through
`use_fixed_tc_perp` and `fixed_tc_perp_value` for the perpendicular component. If you use the realistic profile
you should set unit normalisations.

### Radiative cooling
Set `radiative_cooling = .true.` and specify the cooling curve, see the [physicslist](../parameter_file/#physicslist)
for possible options. For radiative cooling we only have a temperature-dependent profile due to the use of realistic
cooling tables, hence you should set unit normalisations when you use this.

### Flow
Set `flow = .true.` and specify the various velocity components in the corresponding fields. See below.

### External gravity
Set `external_gravity = .true.` and specify the gravity component in the corresponding field. See below.

**Tip:** The variable `g` in the `mod_equilibrium_params` module is usually used if you want to set a constant gravity
value.
{: .notice--success}

## Defining the equilibrium state

**Please note** <i class="fas fa-hard-hat"></i> <i class="fas fa-hammer"></i>  
We are planning a (minor) refactor of some attributes in the (near) future, so keep a look out for changes.
{: .notice--danger}

**Tip:** to make sure your implementation is in order you can do a dry run first (`dry_run = .true.` in the parfile)
and plot the equilibrium arrays using Pylbo. Once you're sure that everything is correct, remove the `dry_run` line from
the parfile and you're good to go.
{: .notice--success}

The final thing that should be done is specifying the actual equilibrium state. The `mod_equilibrium` parent module
contains various different types (we call them _fields_) each containing multiple attributes. For example,
the variable `rho_field` contains (you guessed it) density-related variables, accessible through
`rho_field % rho0` for the regular density and `rho_field % d_rho0_dr` for the density derivative.
Analogous for the other fields. The advantage of this is that everything is contained to its particular type,
meaning they can be easily passed to other subroutines, and adding a new attribute is as easy as appending it to the
corresponding type. For a full list of the various fields and attributes see
[mod_types](../../src-docs/ford/module/mod_types.html#type-density_type).

You only have to specify the `rho_field`, `v_field`, `B_field`, `T_field` and `grav_field` (the names are self-explanatory);
if a field is not set all attributes are assumed to be zero. For example, if you don't include flow or gravity, you
don't have to set those fields explicitly.

**Warning:** _technically_ you have access to `eta_field`, `rc_field` and `kappa_field` in the submodule, corresponding
to the resistivity, radiative cooling and thermal conduction types. These are set automatically by Legolas
and should _not_ be set manually.
{: .notice--warning}

Since equilibrium expressions can be quite complicated, it is best to set everything in a loop where you iterate
over the _Gaussian_ grid, NOT over the _base_ grid. The base grid is only used to set up the Gaussian grid,
which is where the equilibrium quantities are evaluated in anyway. All field attributes are 1D-arrays,
so you can simply iterate over them when setting the equilibrium.
In the example below we include flow and a radius-dependent external gravity profile.

**Note:** as indicated in the warning above you should not set `eta_field` manually. There is one exception through,
which concerns the second magnetic field derivatives. These are only used when resistivity is included (hence why they are
a part of `eta_field` and not of `B_field`) and should be set explicitly if `resistivity = .true.`.
You can do this through `eta_field % dd_B02_dr` and `eta_field % dd_B03_dr`.
{: .notice--info}

```fortran
! specify variables to use later on
real(dp)    :: x
integer     :: i

! set the grid variables
geometry = "Cartesian"
x_start = 0.0d0
x_end = 2.0d0
call initialise_grid()

! loop over the gaussian points
do i = 1, gauss_gridpts
  ! evaluate equilibria in the Gaussian grid
  x = grid_gauss(i)

  ! set equilibria + derivatives
  rho_field % rho0(i) = 5.0d0 - x**2
  rho_field % d_rho0_dr(i) = -2.0d0 * x
  T_field % T0(i) = 2.0d0 ! constant, so no need to set d_T0_dr
  v_field % v02(i) = x
  v_field % d_v02_dr(i) = 1.0d0
  B_field % B02(i) = sin(x)
  B_field % d_B02_dr(i) = cos(x)
  B_field % B03(i) = cos(x)
  B_field % d_B03_dr(i) = -sin(x)
  B_field % B0(i) = sqrt(B_field % B02(i)**2 + B_field % B03(i)**2)
  ! no resistivity, so no second B_field derivatives needed
  grav_field % grav(i) = 1.0d0 / x**2
end do
```
The kind `dp` above denotes the intrinsic double precision `real64` from Fortran `iso_fortran_env`, defined in
`mod_global_variables`. The variables `dp`, `gauss_gridpts` and `grid_gauss` are known in the parent module and
don't have to be imported.

If you prefer to work with pressure instead of temperature that is also possible.
You can specify a pressure and density profile, and set temperature as such:
```fortran
real(dp) :: x
! pressure arrays
real(dp) :: px(gauss_gridpts), dpx(gauss_gridpts)
integer  :: i

do i = 1, gauss_gridpts
  x = grid_gauss(i)
  rho_field % rho0(i) = 5.0d0 - x**2
  rho_field % d_rho0_dr(i) = -2.0d0 * x
  px(i) = 0.5d0 * x
  dpx(i) = 0.5d0
  ! ideal gas law to set temperature, T' = (p'*rho - rho'*p) / rho**2
  T_field % T0(i) = px(i) / rho_field % rho0(i)
  T_field % d_T0_dr(i) = (  dpx(i) * rho_field % rho0(i) &
                            - rho_field % d_rho0_dr(i) * px(i)  ) &
                         / rho_field % rho0(i)**2
end do
```

For more examples you can look at implementation of the various descendants of the
[`mod_equilibrium`](../../src-docs/ford/module/mod_equilibrium.html) parent module.

**Note:** all equilibria should satisfy the [equilibrium conditions](../equations_physics/#equilibrium-requirements).
The examples here are only for illustration purposes and do not satisfy those conditions. In principle Legolas
will run if your equilibrium does not check out, but you _will_ be warned.
{: .notice--info}
