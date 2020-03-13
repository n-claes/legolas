!
! SUBMODULE: smod_equil_adiabatic_homo
!
! DESCRIPTION:
!> Submodule defining an adiabatic homogeneous medium in Cartesian geometry.
submodule (mod_equilibrium) smod_equil_adiabatic_homo
  implicit none

contains

  module subroutine adiabatic_homo_eq()

    call allow_geometry_override(default_geometry="Cartesian", default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      k2 = 0
      k3 = dpi
    end if

    rho_field % rho0 = 1.0d0
    T_field % T0     = 1.0d0
    B_field % B02    = 0.0d0
    B_field % B03    = 1.0d0
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03)**2)

  end subroutine adiabatic_homo_eq

end submodule smod_equil_adiabatic_homo
