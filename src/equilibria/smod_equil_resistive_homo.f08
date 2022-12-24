! =============================================================================
!> This submodule defines a simple, homogeneous medium in Cartesian
!! geometry with a constant resistivity value. The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from section 14.3, p. 550 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : used to set the density value.
!! - <tt>cte_B02</tt> = 0 : used to set the By value.
!! - <tt>cte_B03</tt> = 1 : used to set the Bz value.
!! - <tt>beta</tt> = 0.25 : used to constrain the temperature value.
!! - fixed resistivity value of 0.001
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_resistive_homo
  implicit none

contains

  !> Sets the equilibrium.
  module procedure resistive_homo_eq
    use mod_equilibrium_params, only: beta, cte_rho0, cte_B02, cte_B03

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_resistivity(fixed_resistivity_value=0.001_dp)

      k2 = 0.0d0
      k3 = 1.0d0
      beta = 0.25d0
      cte_rho0 = 1.0d0
      cte_B02 = 0.0d0
      cte_B03 = 1.0d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    rho_field % rho0 = cte_rho0
    B_field % B02    = cte_B02
    B_field % B03    = cte_B03
    B_field % B0     = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    T_field % T0     = beta * (B_field % B0)**2 / (2.0d0)
  end procedure resistive_homo_eq

end submodule smod_equil_resistive_homo
