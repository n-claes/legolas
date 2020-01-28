!
! SUBMODULE: smod_equil_interface_modes
! 
! DESCRIPTION:
!> Submodule defining ... TODO
submodule (mod_equilibrium) smod_equil_interface_modes
  implicit none
  
contains 
  
  module subroutine interface_modes_eq()
    real(dp)      :: x, B0, B_e, rho0, rho_e, T0, T_e
    integer       :: i
    
    geometry = 'Cartesian'
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    !! Parameters
    rho_e = 1.0d0
    T_e   = 1.0d0
    B_e   = 1.0d0
    rho0  = 0.0d0
    T0    = 0.0d0
    B0    = sqrt(2*rho_e*T_e+B_e**2)

    k2  = 0.0d0
    k3  = 1.0d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      !! Equilibrium
      if (x > 0) then
        rho_field % rho0(i) = rho_e
        T_field % T0(i) = T_e
        B_field % B03(i)= B_e
      else
        rho_field % rho0(i)= rho0
        T_field % T0(i) = T0
        B_field % B03(i)= B0
      end if
    end do
  end subroutine interface_modes_eq
  
end submodule smod_equil_interface_modes
