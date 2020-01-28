!
! SUBMODULE: smod_equil_rotating_plasma_cylinder
! 
! DESCRIPTION:
!> Submodule defining a rotating plasma cylinder in cylindrical geometry.
!! From Nijboer et al., J Plasma Phys 58(1) (1997).
submodule (mod_equilibrium) smod_equil_rotating_plasma_cylinder
  implicit none
  
contains 
  
  module subroutine rotating_plasma_cyl_eq()
    real(dp)    :: a21, a22, a3, b21, b22, b3, p0
    real(dp)    :: r
    integer     :: i

    geometry = 'cylindrical'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid() ! Initialise grid

    flow = .true.

    !! Parameters
    a21 = 8.0d0
    a22 = 0.0d0
    a3  = 0.0d0
    b21 = 1.0d0
    b22 = 0.0d0
    b3  = 0.0d0
    p0  = 0.1d0

    k2  = 1.0d0
    k3  = 0.0d0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      !! Equilibrium
      rho_field % rho0(i) = 1.0d0
      v_field % v02(i) = a21*r + a22*r**2
      v_field % v03(i) = a3
      B_field % B02(i) = b21*r + b22*r**2
      B_field % B03(i) = b3
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = (1.0d0 / (rho_field % rho0(i))) &
                        * (p0 + 0.5d0 * (a21**2 - 2.0d0*b21**2)*r**2 &
                            + (2.0d0/3.0d0)*(a21*a22 - b21*b22)*r**3 &
                            + (1.0d0/4.0d0)*(a22**2 - b22**2)*r**4)
      !! Derivatives
      B_field % d_B02_dr(i) = b21 + 2.0d0*b22*r
      B_field % d_B03_dr(i) = 0.0d0
      v_field % d_v02_dr(i) = a21 + 2.0d0*a22*r
      v_field % d_v03_dr(i) = 0.0d0
      T_field % d_T0_dr(i) = (1.0d0 / (rho_field % rho0(i))) * ( &
                                (a21**2 - 2.0d0*b21**2)*r + 2.0d0*(a21*a22 - b21*b22)*r**2 &
                                + (a22**2 - b22**2)*r**3   )
      !TODO eta?
      eta_field % dd_B02_dr = 2.0d0*b22
    end do
  
  end subroutine rotating_plasma_cyl_eq
  
end submodule smod_equil_rotating_plasma_cylinder
