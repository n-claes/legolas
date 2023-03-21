! =============================================================================
!> This submodule defines a cylindrical equilibrium, resembling a
!! rotating plasma cylinder. The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from
!! _Nijboer, R. J., Holst, B., Poedts, S., & Goedbloed, J. P. (1997).
!!  Calculating magnetohydrodynamic flow spectra. Computer physics communications, 106(1-2), 39-52_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 0
!! - <tt>cte_rho0</tt> = 1 : used to set the density value.
!! - <tt>cte_p0</tt> = 0.1 : used to set the pressure.
!! - <tt>p1</tt> = 8 : sets the constant a21.
!! - <tt>p2</tt> = 0 : sets the constant a22.
!! - <tt>p3</tt> = 0 : sets the constant a3.
!! - <tt>p4</tt> = 1 : sets the constant b21.
!! - <tt>p5</tt> = 0 : sets the constant b22.
!! - <tt>p6</tt> = 0 : sets the constant b3.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_rotating_plasma_cylinder
  use mod_equilibrium_params, only: cte_rho0, cte_p0, p1, p2, p3, p4, p5, p6
  implicit none

  real(dp)    :: a21, a22, a3, b21, b22, b3

contains

  module procedure rotating_plasma_cyl_eq
    real(dp)    :: r
    integer     :: i

    if (settings%equilibrium%use_defaults) then  ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_flow()

      a21 = 8.0_dp
      a22 = 0.0_dp
      a3  = 0.0_dp
      b21 = 1.0_dp
      b22 = 0.0_dp
      b3  = 0.0_dp
      cte_rho0 = 1.0
      cte_p0  = 0.1_dp

      k2  = 1.0_dp
      k3  = 0.0_dp
    else  ! LCOV_EXCL_STOP
      a21 = p1
      a22 = p2
      a3 = p3
      b21 = p4
      b22 = p5
      b3 = p6
    end if
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      rho_field % rho0(i) = rho0()
      v_field % v02(i) = v02(r)
      v_field % v03(i) = v03()
      B_field % B02(i) = B02(r)
      B_field % B03(i) = B03()
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = T0(r)

      B_field % d_B02_dr(i) = dB02(r)
      v_field % d_v02_dr(i) = dv02(r)
      T_field % d_T0_dr(i)  = dT0(r)
    end do
  end procedure rotating_plasma_cyl_eq


  real(dp) function rho0()
    rho0 = cte_rho0
  end function rho0

  real(dp) function T0(r)
    real(dp), intent(in) :: r
    T0 = (1.0_dp / rho0()) * ( &
      cte_p0 + 0.5_dp * (a21**2 - 2.0_dp * b21**2) * r**2 &
      + (2.0_dp / 3.0_dp) * (a21 * a22 - b21 * b22) * r**3 &
      + (1.0_dp / 4.0_dp) * (a22**2 - b22**2) * r**4 &
    )
  end function T0

  real(dp) function dT0(r)
    real(dp), intent(in) :: r
    dT0 = (1.0_dp / rho0()) * ( &
      (a21**2 - 2.0_dp * b21**2) * r &
      + 2.0_dp * (a21 * a22 - b21 * b22) * r**2 &
      + (a22**2 - b22**2) * r**3 &
    )
  end function dT0

  real(dp) function v02(r)
    real(dp), intent(in) :: r
    v02 = a21 * r + a22 * r**2
  end function v02

  real(dp) function dv02(r)
    real(dp), intent(in) :: r
    dv02 = a21 + 2.0_dp * a22 * r
  end function dv02

  real(dp) function v03()
    v03 = a3
  end function v03

  real(dp) function B02(r)
    real(dp), intent(in) :: r
    B02 = b21*r + b22 * r**2
  end function B02

  real(dp) function dB02(r)
    real(dp), intent(in) :: r
    dB02 = b21 + 2.0_dp * b22 * r
  end function dB02

  real(dp) function B03()
    B03 = b3
  end function B03

end submodule smod_equil_rotating_plasma_cylinder
