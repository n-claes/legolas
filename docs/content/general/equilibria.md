---
title: Implemented equilibria
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2020-09-03
---

For a complete list of implemented equilibria we refer to the [Legolas docs](../../src-docs/ford/lists/modules.html).
Look at the documentation of every submodule for information on default parameters, which one you can vary and where
it comes from. Below we give a small overview of which equilibria are implemented together with the included 
physical effects, and links to the source docs.

**Pure adiabatic**
- [Homogeneous medium](../../src-docs/ford/module/smod_equil_adiabatic_homo.html): _Cartesian_
  ```fortran 
  equilibrium_type="adiabatic_homo"
  ```
- [Flux tube under coronal conditions](../../src-docs/ford/module/smod_equil_coronal_flux_tube.html): _cylindrical_
  ```fortran 
  equilibrium_type="coronal_flux_tube"
  ```
- [Flux tube under photospheric conditions](../../src-docs/ford/module/smod_equil_photospheric_flux_tube.html): _cylindrical_
  ```fortran 
  equilibrium_type="photospheric_flux_tube"
  ```
- [Constant axial current](../../src-docs/ford/module/smod_equil_constant_current.html): _cylindrical_
  ```fortran 
  equilibrium_type="constant_current_tokamak"
  ```

**Adiabatic + flow**
- [Flow driven instabilities](../../src-docs/ford/module/smod_equil_flow_driven_instabilities.html): _Cartesian_
    - [Kelvin-Helmholtz HD instabilities](../../src-docs/ford/module/smod_equil_KHI.html)
      ```fortran 
      equilibrium_type="constant_current_tokamak"
      ```
- [Internal kink modes in a force-free magnetic field](../../src-docs/ford/module/smod_equil_internal_kink_instability.html): _cylindrical_
  ```fortran 
  equilibrium_type="internal_kink"
  ```
- [Kelvin-Helmholtz and current-driven instabilities](../../src-docs/ford/module/smod_equil_kelvin_helmholtz_cd.html): _cylindrical_
  ```fortran 
  equilibrium_type="kelvin_helmholtz_cd"
  ```
- [Rotating plasma cylinder](../../src-docs/ford/module/smod_equil_rotating_plasma_cylinder.html): _cylindrical_
  ```fortran 
  equilibrium_type="rotating_plasma_cylinder"
  ```
- [Rayleigh-Taylor instabilities in a rotating theta-pinch](../../src-docs/ford/module/smod_equil_RTI_theta_pinch.html): _cylindrical_
  ```fortran 
  equilibrium_type="RTI_theta_pinch"
  ```
- [Suydam cluster modes](../../src-docs/ford/module/smod_equil_suydam_cluster.html): _cylindrical_
  ```fortran 
  equilibrium_type="suydam_cluster"
  ```
    
**Adiabatic + external gravity**
- [Gravito-acoustic HD modes](../../src-docs/ford/module/smod_equil_gravito_acoustic.html): _Cartesian_
  ```fortran 
  equilibrium_type="gravito_acoustic"
  ```
- [Gravito-acoustic MHD modes](../../src-docs/ford/module/smod_equil_gravito_mhd.html): _Cartesian_
  ```fortran 
  equilibrium_type="gravito_mhd"
  ```
- [Interchange modes](../../src-docs/ford/module/smod_equil_interchange_modes.html): _Cartesian_
  ```fortran 
  equilibrium_type="interchange_modes"
  ```

**Adiabatic + external gravity + flow**
- [Flow driven instabilities](../../src-docs/ford/module/smod_equil_flow_driven_instabilities.html): _Cartesian_
    - [Rayleigh-Taylor instabilities](../../src-docs/ford/module/smod_equil_RTI.html)
      ```fortran 
      equilibrium_type="rayleigh_taylor"
      ```
    - [Kelvin-Helmholtz + Rayleigh-Taylor instabilities](../../src-docs/ford/module/smod_equil_RTI_KHI.html)
      ```fortran 
      equilibrium_type="RTI_KHI"
      ```
- [Magneto-rotational instabilities](../../src-docs/ford/module/smod_equil_MRI.html): _cylindrical accretion disk_
  ```fortran 
  equilibrium_type="MRI_accretion"
  ```

**Pure resistive**
- [Resistive homogeneous medium](../../src-docs/ford/module/smod_equil_resistive_homo.html): _Cartesian_
  ```fortran 
  equilibrium_type="resistive_homo"
  ```
- [Resistive tearing modes](../../src-docs/ford/module/smod_equil_resistive_tearing.html): _Cartesian_
  ```fortran 
  equilibrium_type="resistive_tearing"
  ```
- [Resonant absorption](../../src-docs/ford/module/smod_equil_resonant_absorption.html): _Cartesian_
  ```fortran 
  equilibrium_type="resonant_absorption"
  ```

**Resistive + flow**
- [Resistive tearing modes with flow](../../src-docs/ford/module/smod_equil_resistive_tearing_flow.html): _Cartesian_
  ```fortran 
  equilibrium_type="resistive_tearing_flow"
  ```

**Non-adiabatic**
- [Discrete Alfvén waves](../../src-docs/ford/module/smod_equil_discrete_alfven.html): _cylindrical_
  ```fortran 
  equilibrium_type="discrete_alfven"
  ```
- [Gold-Hoyle](../../src-docs/ford/module/smod_equil_gold_hoyle.html): _cylindrical_
  ```fortran 
  equilibrium_type="gold_hoyle"
  ```
- [Magnetothermal instabilities](../../src-docs/ford/module/smod_equil_magnetothermal_instabilities.html):  _cylindrical_
  ```fortran 
  equilibrium_type="magnetothermal_instabilities"
  ```