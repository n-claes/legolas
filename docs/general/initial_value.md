---
title: Initial-value solver
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "chevron-circle-down"
last_modified_at: 2026-03-09
---

Legolas includes an initial-value solver that integrates perturbations forward in time directly from the FEM matrices, without requiring a separate code.
Given a set of initial perturbation profiles, it evolves the system

$$B \frac{d\mathbf{x}}{dt} = -iA\mathbf{x}$$

using the implicit theta-method (implicit midpoint at the default `alpha = 0.5`), and saves snapshots of the solution at regular intervals alongside the standard eigenvalue output.

This page assumes familiarity with how to set up and run Legolas. If not, see [running your first problem](../../getting-started/running) and [implementing a custom setup](../../general/own_setup).

## Supported physics types

The initial-value solver currently supports:

| Physics type | Perturbed components |
|---|---|
| `isothermal-1d` | $\rho_1$, $v_1$ |
| `hd-1d` | $\rho_1$, $v_1$, $T_1$ |
| `hd` | $\rho_1$, $v_1$, $v_2$, $v_3$, $T_1$ |

## Configuration

IVP mode is configured through the `ivplist` namelist in the parfile:

```fortran
&ivplist
  enabled         = .true.
  alpha           = 0.5    ! implicitness (0 = forward Euler, 1 = backward Euler, 0.5 = implicit midpoint)
  t_end           = 10.0   ! end time
  n_steps         = 1000   ! number of time steps
  n_snapshots     = 100    ! number of snapshots to save
  snapshot_stride = 10     ! save every n-th step (overrides n_snapshots if set)
/
```

Snapshots are written to the datfile automatically when `enabled = .true.`.

## Specifying initial conditions

Initial conditions are set in the `user_defined_eq` subroutine of your `smod_user_defined.f08` file.
An `initial_conditions` object is passed in alongside the standard `settings`, `grid`, `background`, and `physics` objects.
You set perturbation profiles using type-bound setter routines — each takes a function (and optionally its derivative) matching the signature `f(x) result(y)` where `x` and `y` are both `real(dp)` arrays of the same size.

The physics type must be set in the parfile (not in the Fortran submodule):
```fortran
&physicslist
  physics_type = "isothermal-1d"
/
```

A minimal example with a Gaussian density perturbation:

```fortran
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  module subroutine user_defined_eq(settings, grid, background, physics, initial_conditions)
    type(settings_t), intent(inout)           :: settings
    type(grid_t), intent(inout)               :: grid
    type(background_t), intent(inout)         :: background
    type(physics_t), intent(inout)            :: physics
    type(initial_conditions_t), intent(inout) :: initial_conditions

    ! --- equilibrium ---
    background%density%rho0   => uniform_density
    background%temperature%T0 => uniform_temperature

    ! --- initial conditions ---
    call initial_conditions%set_ic_density_funcs(rho_func=gaussian_rho)

  end subroutine user_defined_eq


  pure function uniform_density(x) result(rho)
    real(dp), intent(in) :: x(:)
    real(dp) :: rho(size(x))
    rho = 1.0_dp
  end function uniform_density

  pure function uniform_temperature(x) result(T)
    real(dp), intent(in) :: x(:)
    real(dp) :: T(size(x))
    T = 1.0_dp
  end function uniform_temperature

  pure function gaussian_rho(x) result(rho1)
    real(dp), intent(in) :: x(:)
    real(dp) :: rho1(size(x))
    real(dp), parameter :: x0 = 0.5_dp, sigma = 0.05_dp
    rho1 = exp(-((x - x0) / sigma)**2)
  end function gaussian_rho

end submodule smod_user_defined
```

The setter routines available on the `initial_conditions` object are:

| Subroutine | Component |
|---|---|
| `set_ic_density_funcs(rho_func [, drho_func])` | $\rho_1$ |
| `set_ic_velocity_1_funcs(v01_func, dv01_func)` | $v_1$ |
| `set_ic_velocity_2_funcs(v02_func, dv02_func)` | $v_2$ |
| `set_ic_velocity_3_funcs(v03_func, dv03_func)` | $v_3$ |
| `set_ic_temperature_funcs(T_func [, dT_func])` | $T_1$ |

Components not set default to zero. Components not present in the chosen physics type are silently ignored.
It is not necessary to set the density or temperature derivatives when using the default basis functions, therefore
arguments are optional.

## Post-processing with Pylbo

Once Legolas has run, load the datfile and retrieve the snapshots:

```python
import pylbo
import matplotlib.pyplot as plt

ds = pylbo.load("output/my_datfile.dat")

if ds.has_iv_snapshots:
    ivp = ds.get_iv_snapshots()
```

The returned `IVPSolution` object provides several methods for inspection:

### Space-time heatmap

```python
fig, ax = plt.subplots()
ivp.plot_space_time_heatmap("rho", ax=ax)
plt.show()
```

### Spatial profiles at selected snapshots

```python
fig, ax = plt.subplots()
ivp.plot_spatial_slices("rho", snap_indices=[0, 25, 50, 99], ax=ax)
plt.show()
```

### Accessing raw data

```python
# shape: (n_snapshots, n_points)
rho_data = ivp.get_component("rho")

# physical times at each snapshot
print(ivp.times)

# spatial coordinate
print(ivp.x_domain)
```
