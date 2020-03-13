!
! SUBMODULE: smod_equil_resistive_homo
!
! DESCRIPTION:
!> Submodule defining a homogeneous medium in Cartesian geometry with a
!! constant resistivity value. From Advanced Magnetohydrodynamics by Goedbloed,
!! Keppens and Poedts, page 156.
submodule (mod_equilibrium) smod_equil_resistive_homo
  implicit none

contains

  module subroutine resistive_homo_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: beta

    call allow_geometry_override(default_geometry="Cartesian", default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 0.001d0

    if (use_defaults) then
      k2 = 0.0d0
      k3 = 1.0d0
      beta = 0.25d0
    end if

    !! Equilibrium
    rho_field % rho0 = 1.0d0
    B_field % B02    = 0.0d0
    B_field % B03    = 1.0d0
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    ! n = 1, kB = 1, mu0 = 1 for T0
    T_field % T0     = beta * (B_field % B0)**2 / (2.0d0)

  end subroutine resistive_homo_eq

end submodule smod_equil_resistive_homo
