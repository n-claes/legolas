---
title: Running your own setup
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "tools"
toc_label: "User submodule configuration"
last_modified_at: 2021-07-26
---

This page explains how to configure the custom user submodule to your needs. We assume that you are already familiar with
how to set up a legolas directory and how to run legolas in a separate directory. If this is not the case see [how to run your first problem](../../getting-started/running).
In what follows we also assume that both a `smod_user_defined.f08` file and parfile are present in the current working directory. You can copy both of these over when running
the setup script.

## The submodule itself
The submodule `smod_user_defined.f08` is used to set up your own equilibrium configuration. The general template looks like this:
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
which contains the interface of the subroutines that are implemented in their respective submodules. This also implies that you have access to all public variables/routines/functions
that are either defined or used in the [`mod_equilibium`](../../ford/module/mod_equilibrium.html) module through host association, without having to `use` them explicitly.

In what follows all code samples go _inside_ the `user_defined_eq()` subroutine, since this is the only routine you
need to configure yourself.

## Setting the grid
Here you have 2 options: let Legolas take care of everything _or_ provide a custom grid yourself.

**Note:** The variables `geometry`, `x_start`, `x_end` and the subroutines `initialise_grid` and `allow_geometry_override` are all known
in the parent module and don't have to be used explicitly.
{: .notice--info}
### Using a default grid
In this case you supply the geometry, start and end of the grid, and Legolas will create a uniform grid for you. Here you can either explicitly set
the `geometry`, `x_start` and `x_end` variables to some values, or you use the `allow_geometry_override` subroutine:
```fortran
call allow_geometry_override(default_geometry="Cartesian", default_x_start=0.0d0, default_x_end=1.0d0)
call initialise_grid()
```
This allows you to vary the grid properties through the parfile, without having to recompile your module. Legolas will default to the values
given here if nothing is specified in the parfile.
The call to `initialise_grid` sets up the (uniform) grid, Gaussian grid and scale factors.


### Using a custom grid
The second possibility is defining a grid yourself. There are no restrictions on the custom grid (besides it needing to be strictly increasing),
meaning you have full control over start, end and grid spacing. This can be particularly useful if you want mesh accumulation near a certain region to better resolve some
sharp gradients instead of increasing the base resolution. Simply define an array of size `gridpts`, fill it and supply it to the initialisation routine:
```fortran
real(dp) :: my_grid(gridpts)
integer :: i

x_start = 0.0d0
x_end = 2.0d0
do i = 1, gridpts
  ! fill the grid however you want
end do
call initialise_grid(custom_grid=my_grid)
```
The `dp` here is a global parameter, corresponding to the double precision value `real64` from `iso_fortran_env`. The initialisation routine takes care of
setting up the Gaussian grid automatically, using the custom grid you provided. You can even mix this with [parfile variables](../../general/parameter_file/#paramlist)
to change the grid between different runs. You can take a look at the pre-implemented
[flux tube equilibria](../../ford/sourcefile/smod_equil_coronal_flux_tube.f08.html) for inspiration, which use a separate grid spacing in 3 different regions.

## Using variables
Legolas uses a dedicated module which contains a multitude of variables that you can use in your setups.
For a comprehensive list we refer to [this page](../../ford/module/mod_equilibrium_params.html#variable-k2).
You can use all of these in your setups, but note that these should be explicitly imported (except for `k2` and `k3`).
An example is given below.

**Note:** The variables `use_defaults`, `k2` and `k3` are known in the parent module. All other parameters should be explicitly used.
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
[`mod_units`](../../ford/module/mod_units.html#variable-unit_length) are known to the parent module and
do not have to be used explicitly in the submodule.
{: .notice--info}

Specifying the unit normalisations can either be done through the parfile (the `unitslist`) or in the submodule
directly. This is only relevant if you include radiative cooling, thermal conduction or (temperature-dependent)
resistivity, since these units are used to re-dimensionalise the variables for the corresponding calculations and/or
table lookups. To set the units you can do
```fortran
use_cgs = .true.
call set_normalisations( &
  new_unit_density=1.6727d-15, &
  new_unit_magneticfield=22.5d0, &
  new_unit_length=1.0d10, &
  new_mean_molecular_weight=1.0d0 &  ! this one is optional
)
```
which sets values in cgs units, if you set `use_cgs = .false.` these are interpreted as SI units. Instead of specifying
`new_unit_density` you can choose a reference temperature `new_unit_temperature` instead. Setting the mean molecular weight is optional,
the default value is 1.

## Including additional physics
**Note:** Don't forget to enable additional physics in the parfile's `physicslist`. Another option is enabling them in the submodule itself.
{: .notice--info}
### Resistivity
In the case of resistivity you have two options:
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

If you only set `thermal_conduction = .true.` this will include both parallel and perpendicular conduction based on the
Spitzer conductivity. If you want to omit perpendicular thermal conduction, set it to zero using the fixed value flags.

### Radiative cooling
Set `radiative_cooling = .true.` and specify the cooling curve, see the [physicslist](../parameter_file/#physicslist)
for possible options. For radiative cooling we only have a temperature-dependent profile due to the use of realistic
cooling tables, hence you should set unit normalisations when you use this.

### Flow
Set `flow = .true.` and specify the various velocity components in the corresponding fields. See below.

### External gravity
Set `external_gravity = .true.` and specify the gravity component in the corresponding field. See below.

## Defining the equilibrium configuration

<i class="fa fa-file-contract" aria-hidden="true"></i>
**Disclaimer**: due to the finite element representation used by Legolas the equilibrium profiles (or their derivatives) do not have to be
continuous. Do note though that smooth, continuous profiles give the best results, and you should try to achieve this as much as possible. Instead of a discontinuity
you can try to connect different regions with a smooth sine or hyperbolic tangent profile. This is not at all required and Legolas will run just fine with strong jumps or gradients,
but know that in that case you should double check if the results make sense and ensure enough gridpoints to resolve everything.
{: .notice--danger}

The final thing that should be done is specifying the actual equilibrium state. The `mod_equilibrium` parent module
contains various types (we call them _fields_) each containing multiple attributes. For example,
the variable `rho_field` contains (you guessed it) density-related variables, accessible through
`rho_field % rho0` for the regular density and `rho_field % d_rho0_dr` for the density derivative.
Analogous for the other fields. All attributes are arrays with the same size as the Gaussian grid. For a full list of the various fields and attributes see
[mod_types](../../ford/module/mod_types.html#type-density_type).

You only have to specify the `rho_field`, `v_field`, `B_field`, `T_field` and `grav_field` (the names are self-explanatory);
if a field is not set all attributes are assumed to be zero. For example, if you don't include flow or gravity, you
don't have to set those fields explicitly.

Note that you have to use the **Gaussian** grid and NOT the **base** grid!
The base grid is only used to set up the Gaussian grid, and only the latter is used to define the equilibrium in. When setting the equilibrium attributes
you can either use Numpy-like array operations using `gaussian_grid`, or you iterate over every gridpoint individually. Here we give the resistive tearing
modes as an example, you can take a look at the various pre-implemented equilibria `legolas/equilibria` folder for more inspiration.

<i class="fa fa-exclamation-triangle" aria-hidden="true"></i>
**Note:** all equilibria you implement have to satisfy force balance, and in the case of non-adiabatic effects also thermal balance.
Legolas will warn you if one of the conditions is not satisfied. Note that you can run your setup if this is not the case, but know that there is a high probability that the results
don't make much sense. See [this page](../../physics/balance_reqs/#equilibrium-requirements) for more information.
{: .notice--warning}

```fortran
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  module subroutine user_defined_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: alpha, beta, cte_rho0

    geometry = "Cartesian"
    x_start = -0.5d0
    x_end = 0.5d0
    call initialise_grid()

    if (use_defaults) then
      resistivity = .true.
      use_fixed_resistivity = .true.
      fixed_eta_value = 1.0d-4

      k2 = 0.49d0
      k3 = 0.0d0

      alpha = 4.73884d0
      beta  = 0.15d0
      cte_rho0 = 1.0d0
    end if

    rho_field % rho0 = cte_rho0
    B_field % B02 = sin(alpha * grid_gauss)
    B_field % B03 = cos(alpha * grid_gauss)
    B_field % B0 = 1.0d0  ! sqrt(B02**2 + B03**2)
    T_field % T0 = 0.5d0 * beta

    B_field % d_B02_dr = alpha * cos(alpha * grid_gauss)
    B_field % d_B03_dr = -alpha * sin(alpha * grid_gauss)
    eta_field % dd_B02_dr = -alpha**2 * sin(alpha * grid_gauss)
    eta_field % dd_B03_dr = -alpha**2 * cos(alpha * grid_gauss)
  end subroutine user_defined_eq
end submodule smod_user_defined
```

<i class="fas fa-lightbulb" aria-hidden="true"></i>
**Tip:** to make sure your implementation is in order you can do a dry run first (`dry_run = .true.` in the `savelist`)
and plot the equilibrium arrays using Pylbo. Once you're sure that everything is correct, remove the `dry_run` line from
the parfile and you're good to go.
{: .notice--info}