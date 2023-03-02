! =============================================================================
!> This submodule defines an equilibrium in Cartesian geometry with a constant
!! resistivity value. Parameters are taken in such a way as to allow for
!! resistive tearing modes. The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from section 14.3, p. 551 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0.49
!! - <tt>k3</tt> = 0
!! - <tt>cte_rho0</tt> = 1 : used to set the density value.
!! - <tt>alpha</tt> = 4.73884 : parameter in the magnetic field prescription.
!! - <tt>beta</tt> = 0.25 : used to constrain the temperature value.
!! - fixed resistivity value of 0.0001
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_resistive_tearing
  implicit none

contains

  !> Sets the equilibrium
  module procedure resistive_tearing_modes_eq
    use mod_equilibrium_params, only: alpha, beta, cte_rho0

    real(dp)              :: x
    integer               :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(-0.5_dp, 0.5_dp)
      call settings%physics%enable_resistivity(fixed_resistivity_value=0.0001_dp)

      k2 = 0.49d0
      k3 = 0.0d0

      alpha = 4.73884d0
      beta  = 0.15d0
      cte_rho0 = 1.0d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      B_field % B02(i)    = sin(alpha * x)
      B_field % B03(i)    = cos(alpha * x)
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = beta * (B_field % B0(i))**2 / (2.0d0)

      B_field % d_B02_dr(i) = alpha * cos(alpha * x)
      B_field % d_B03_dr(i) = -alpha * sin(alpha * x)
      ! No d_T0_dr needed, as B0**2 is independent of r

      eta_field % dd_B02_dr(i) = -alpha**2 * sin(alpha * x)
      eta_field % dd_B03_dr(i) = -alpha**2 * cos(alpha * x)
    end do
  end procedure resistive_tearing_modes_eq

end submodule smod_equil_resistive_tearing
