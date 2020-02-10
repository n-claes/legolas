!
! SUBMODULE: smod_equil_kelvin_helmholtz_cd
!
! DESCRIPTION:
!> Submodule defining an unperturbed magnetised jet model in cylindrical geometry.
!! Obtained from Baty & Keppens, Astrophys. J 580 (2002).
submodule (mod_equilibrium) smod_equil_kelvin_helmholtz_cd
  implicit none

contains

  module subroutine kh_cd_instability_eq()
    use mod_equilibrium_params, only: V, cte_p0, Bz0, rc, Bth0, rj

    real(dp)    :: r, a
    integer     :: i

    ! Jet radius, other parameters in function of rj
    rj    = 1.0d0

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 2.0d0 * rj
    call initialise_grid() ! Initialise grid

    flow = .true.

    if (use_defaults) then
      V = 1.63d0
      cte_p0 = 1.0d0
      Bz0 = 0.25d0

      ! UNI
      !rc    = 1.0d0
      !Bth0  = 0.0d0

      ! HEL1
      !rc    = 2.0d0
      !Bth0  = 0.4d0 * (rc**2+rj**2) / (rj*rc)

      ! HEL2
      rc    = 0.5d0
      Bth0  = 0.4d0 * (rc**2+rj**2) / (rj*rc)

      k2  = -1.0d0
    end if

    a = 0.1d0 * rj
    k3  = dpi / rj

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      !! Equilibrium
      rho_field % rho0(i) = 1.0d0
      v_field % v03(i)    = (V/2.0d0) * tanh((rj-r)/a)
      B_field % B02(i)    = Bth0 * r*rc / (rc**2 + r**2 )
      B_field % B03(i)    = Bz0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = cte_p0 / (rho_field % rho0(i)) - (Bth0**2/(2.0d0*(rho_field % rho0(i)))) &
                              * (1 - rc**4/(rc**2+r**2)**2)

      !! Derivatives
      B_field % d_B02_dr(i) = Bth0 * rc * (rc**2-r**2) / (r**2+rc**2)**2
      v_field % d_v03_dr(i) = - (V/(2.0d0*a)) / cosh((rj-r)/a)**2
      T_field % d_T0_dr(i)  = - (2.0d0*Bth0**2/(rho_field % rho0(i))) * rc**4*r / (r**2+rc**2)**3
    end do

  end subroutine kh_cd_instability_eq

end submodule smod_equil_kelvin_helmholtz_cd
