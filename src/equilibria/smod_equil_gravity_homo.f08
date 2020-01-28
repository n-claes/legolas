!
! SUBMODULE: smod_equil_gravity_homo
! 
! DESCRIPTION:
!> Submodule defining a homogeneous medium in Cartesian geometry with a 
!! constant gravity term included.
submodule (mod_equilibrium) smod_equil_gravity_homo
  implicit none
  
contains 
  
  module subroutine gravity_homo_eq()
    real(dp)  :: x, g
    integer   :: i

    geometry = 'Cartesian'
    x_start = 0.0d0
    x_end   = 1.0d0
    call initialise_grid()
    
    external_gravity = .true.

    k2 = dpi
    k3 = dpi

    !! Parameters
    B_field % B02 = 1.0d0
    B_field % B03 = 1.0d0
    B_field % B0 = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    T_field % T0 = 1.0d0
    grav_field % grav = 0.5d0
    
    ! Not homogeneous, but consistent with Nijboer equilibrium equation (7)
    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      g = grav_field % grav(i)
      
      rho_field % rho0(i) = exp(-g * x)
      rho_field % d_rho0_dr(i) = -g * exp(-g * x)
    end do
    
  end subroutine gravity_homo_eq
  
end submodule smod_equil_gravity_homo
