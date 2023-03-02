! =============================================================================
!> This submodule defines a simple, adiabatic homogeneous medium in Cartesian
!! geometry. The geometry can be overridden using the parfile.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = \(\pi\)
!! - <tt>cte_rho0</tt> = 1 : used to set the density value.
!! - <tt>cte_T0</tt> = 1 : used to set the temperature value.
!! - <tt>cte_B02</tt> = 0 : used to set the By value.
!! - <tt>cte_B03</tt> = 1 : used to set the Bz value.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_adiabatic_homo
  implicit none

contains

  !> Sets the equilibrium.
  module procedure adiabatic_homo_eq
    use mod_equilibrium_params, only: cte_rho0, cte_T0, cte_B02, cte_B03

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      k2 = 0
      k3 = dpi

      cte_rho0 = 1.0d0
      cte_T0 = 1.0d0
      cte_B02 = 0.0d0
      cte_B03 = 1.0d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    rho_field % rho0 = cte_rho0
    T_field % T0     = cte_T0
    B_field % B02    = cte_B02
    B_field % B03    = cte_B03
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03)**2)
  end procedure adiabatic_homo_eq

end submodule smod_equil_adiabatic_homo
