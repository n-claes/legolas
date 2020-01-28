!
! SUBMODULE: smod_equil_magneto_rotational_instability
! 
! DESCRIPTION:
!> Submodule defining an equilibrium for a magneto-rotational instability in 
!! cylindrical geometry.
!! Obtained from Magnetohydrodynamics (2019), Fig. 13.17
!! Also appears in Goedbloed, Phys. Plasmas 25, 032110 (2018), Fig. 13
submodule (mod_equilibrium) smod_equil_magneto_rotational_instability
  implicit none
  
contains 
  
  module subroutine magneto_rotational_eq()
    real(dp)      :: r, v1, p1, B21, B31
    real(dp)      :: p0(gauss_gridpts), d_p0_dr(gauss_gridpts)
    integer       :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()

    flow = .true.

    !! Parameters
    v1  = 1.0d0
    p1  = 0.01d0
    B21 = 0.01d0
    B31 = 0.01d0

    k2  = 1.0d0
    k3  = 70.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      p0(i)       = p1 / sqrt(r)**5
      d_p0_dr(i)  = -(5.0d0/2.0d0) * p1 / sqrt(r)**7

      !! Equilibrium
      rho_field % rho0(i) = r**(-3.0d0/2.0d0)
      T_field % T0(i) = p0(i) / (rho_field % rho0(i))
      v_field % v02(i) = v1 / sqrt(r)
      B_field % B02(i) = B21 / r**(5.0d0/4.0d0)
      B_field % B03(i) = B31 / r**(5.0d0/4.0d0)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      !! Derivatives
      rho_field % d_rho0_dr(i) = -(3.0d0/2.0d0) * r**(-5.0d0/2.0d0)
      T_field % d_T0_dr(i) = (d_p0_dr(i) * (rho_field % rho0(i)) &
                              - (rho_field % d_rho0_dr(i)) * p0(i)) / (rho_field % rho0(i))**2
      v_field % d_v02_dr(i) = -0.5d0 * v1 / sqrt(r)**3
      B_field % d_B02_dr(i) = -(5.0d0/4.0d0) * B21 / r**(9.0d0/4.0d0)
      B_field % d_B03_dr(i)   = -(5.0d0/4.0d0) * B31 / r**(9.0d0/4.0d0)

      ! TODO: eta?
      eta_field % dd_B02_dr(i)  = (45.0d0/16.0d0) * B21 / r**(13.0d0/4.0d0)
      eta_field % dd_B03_dr(i)  = (45.0d0/16.0d0) * B31 / r**(13.0d0/4.0d0)
    end do
  end subroutine magneto_rotational_eq
  
end submodule smod_equil_magneto_rotational_instability
