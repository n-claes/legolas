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
  implicit none

contains

  module subroutine rotating_plasma_cyl_eq()
    use mod_equilibrium_params, only: cte_rho0, cte_p0, p1, p2, p3, p4, p5, p6

    real(dp)    :: a21, a22, a3, b21, b22, b3
    real(dp)    :: r
    integer     :: i

    call allow_geometry_override( &
      default_geometry="cylindrical", default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then  ! LCOV_EXCL_START
      flow = .true.

      a21 = 8.0d0
      a22 = 0.0d0
      a3  = 0.0d0
      b21 = 1.0d0
      b22 = 0.0d0
      b3  = 0.0d0
      cte_rho0 = 1.0
      cte_p0  = 0.1d0

      k2  = 1.0d0
      k3  = 0.0d0
    else  ! LCOV_EXCL_STOP
      a21 = p1
      a22 = p2
      a3 = p3
      b21 = p4
      b22 = p5
      b3 = p6
    end if

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      v_field % v02(i) = a21*r + a22*r**2
      v_field % v03(i) = a3
      B_field % B02(i) = b21*r + b22*r**2
      B_field % B03(i) = b3
      B_field % B0(i)  = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)  = (1.0d0 / (rho_field % rho0(i))) &
        * (cte_p0 + 0.5d0 * (a21**2 - 2.0d0*b21**2)*r**2 &
        + (2.0d0/3.0d0)*(a21*a22 - b21*b22)*r**3 &
        + (1.0d0/4.0d0)*(a22**2 - b22**2)*r**4)

      B_field % d_B02_dr(i) = b21 + 2.0d0*b22*r
      B_field % d_B03_dr(i) = 0.0d0
      v_field % d_v02_dr(i) = a21 + 2.0d0*a22*r
      v_field % d_v03_dr(i) = 0.0d0
      T_field % d_T0_dr(i)  = (1.0d0 / (rho_field % rho0(i))) * ( &
        (a21**2 - 2.0d0*b21**2)*r &
        + 2.0d0*(a21*a22 - b21*b22)*r**2 &
        + (a22**2 - b22**2)*r**3 &
      )
    end do
  end subroutine rotating_plasma_cyl_eq

end submodule smod_equil_rotating_plasma_cylinder
