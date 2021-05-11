! =============================================================================
!> This submodule defines a magnetic flux tube embedded in a uniform magnetic
!! environment. In this case the flux tube is under photospheric conditions
!! \( c_{Ae} < c_s < c_{se} < c_A \) where the subscript e denotes the outer region.
!! More specifically the equilibrium is defined as
!! \( c_{Ae} = c_s/2, c_{se} = 3c_s/2, c_A = 2c_s \). The geometry can be overridden
!! in the parfile, and is cylindrical by default for \( r \in [0, 10] \).
!!
!! This equilibrium is taken from chapter 6, fig. 6.5 in
!! _Roberts, Bernard (2019). MHD Waves in the Solar Atmosphere.
!! Cambridge University Press._ [DOI](https://doi.org/10.1017/9781108613774).
!! @note For best results, it is recommended to enable mesh accumulation. @endnote
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 2
!! - <tt>cte_rho0</tt> = 1 : density value for the inner tube.
!! - <tt>cte_p0</tt> = 1 : pressure value for the inner tube.
!! - <tt>r0</tt> = 1 : radius of the inner tube.
!!
!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_coronal_flux_tube
submodule(mod_equilibrium) smod_equil_photospheric_flux_tube
  implicit none

contains

  module subroutine photospheric_flux_tube_eq()
    use mod_global_variables, only: dp_LIMIT, gamma
    use mod_equilibrium_params, only: cte_rho0, cte_p0, r0

    real(dp)  :: r, rho_e, p_e, B_0, B_e
    integer   :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=0.0d0, default_x_end=10.0d0)

    if (use_defaults) then
      cte_rho0 = 1.0d0
      cte_p0 = 1.0d0
      r0 = 1.0d0

      k2 = 0.0d0
      k3 = 2.0d0
    end if

    call initialise_grid()

    rho_e = 8.0d0 * (2.0d0 * gamma + 1.0d0) * cte_rho0 / (gamma + 18.0d0)
    p_e = 18.0d0 * (2.0d0 * gamma + 1.0d0) * cte_p0 / (gamma + 18.0d0)
    B_0 = 2.0d0 * sqrt(gamma * cte_p0)
    B_e = sqrt(2.0d0 * gamma * cte_p0 * (2.0d0 * gamma + 1.0d0) / (gamma + 18.0d0))

    if (r0 > x_end) then
      call log_message("equilibrium: inner cylinder radius r0 > x_end", level='error')
    else if (r0 < x_start) then
      call log_message("equilibrium: inner cylinder radius r0 < x_start", level='error')
    end if

    ! check pressure balance
    if (abs(cte_p0 + 0.5d0 * B_0**2 - p_e - 0.5d0 * B_e**2) > dp_LIMIT) then
      error stop "equilibrium: total pressure balance not satisfied"
    end if

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      if (r > r0) then
        rho_field % rho0(i) = rho_e
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_e
        T_field % T0(i) = p_e / rho_e
      else
        rho_field % rho0(i) = cte_rho0
        B_field % B02(i) = 0.0d0
        B_field % B03(i) = B_0
        T_field % T0(i) = cte_p0 / cte_rho0
      end if
      B_field % B0(i)  = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
    end do

  end subroutine photospheric_flux_tube_eq

end submodule smod_equil_photospheric_flux_tube
