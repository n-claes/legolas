!
! SUBMODULE: smod_equil_interface_modes
!
! DESCRIPTION:
!> Submodule defining an equilibrium of a two-layer plasma.
!! Obtained from B. Roberts, MHD Waves in the Solar Atmosphere (2019), Ch. 4
submodule (mod_equilibrium) smod_equil_interface_modes
  implicit none

contains

  module subroutine interface_modes_eq()
    use mod_global_variables, only: gridpts, dp_LIMIT

    real(dp)      :: x, B0, B_e, rho0, rho_e, T0, T_e
    integer       :: i

    geometry = 'Cartesian'
    ! Override values from par file
    x_start = -0.5d0
    x_end   = 0.5d0
    call initialise_grid()

    !! Parameters
    rho_e = 10.0d0
    T_e   = 10.0d0
    B_e   = 10.0d0
    rho0  = 1.0d0
    T0    = 1.0d0
    B0    = sqrt(2.0d0*rho_e*T_e + B_e**2 - 2.0d0*rho0*T0)

    if (abs(rho0*T0 + 0.5d0*B0**2 - rho_e*T_e - 0.5*B_e**2) > dp_LIMIT) then
      stop "Total pressure balance is not satisfied."
    end if

    k2  = 0.0d0
    k3  = 1.0d0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      !! Equilibrium
      if (x > 0) then
        rho_field % rho0(i) = rho_e
        T_field % T0(i)     = T_e
      else
        rho_field % rho0(i) = rho0
        T_field % T0(i)     = T0
      end if

      B_field % B03(i)      = 0.5d0 * (B_e + B0) + 0.5d0 * (B_e - B0) * tanh(gridpts*x)
      B_field % d_B03_dr(i) = 0.5d0 * (B_e - B0) * gridpts / cosh(gridpts*x)**2
    end do

  end subroutine interface_modes_eq

end submodule smod_equil_interface_modes
