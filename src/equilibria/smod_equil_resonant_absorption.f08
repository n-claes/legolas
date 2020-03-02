!
! SUBMODULE: smod_equil_resonant_absorption
!
! DESCRIPTION:
!> Submodule defining an inhomogeneous medium in Cartesian geometry with a constant
!! resistivity value. Two (constant) density profiles are defined which are connected
!! by a sine profile, representing the interface inbetween. This density profile allows
!! for resonant absorption (hence quasi-modes).
!! From Van Doorsselaere and Poedts, Plasma Physics and Controlled Fusion 49.3 (2007), 261.
submodule (mod_equilibrium) smod_equil_resonant_absorption
  implicit none

contains

  module subroutine resonant_absorption_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: p1, p2, r0

    real(dp)  :: x, s, l, rho_left, rho_right
    integer   :: i

    geometry = "Cartesian"
    x_start = 0.0d0
    x_end = 1.0d0
    call initialise_grid()

    resistivity = .true.
    use_fixed_resistivity = .true.
    fixed_eta_value = 10.0d0**(-3.2d0)

    if (use_defaults) then
      k2 = 1.0d0
      k3 = 0.25d0

      rho_left = 9.0d0
      rho_right = 1.0d0 ! so zeta = rho_left/rho_right = 9
      l = 0.94d0
    else
      rho_left = p1
      rho_right = p2
      l = r0
    end if

    s = 0.5d0

    !! Equilibrium
    B_field % B02 = 0.0d0
    B_field % B03 = 1.0d0
    B_field % B0 = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    T_field % T0 = 0.0d0  ! zero pressure

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      if (x >= 0.0d0 .and. x < s - 0.5d0*l) then
        rho_field % rho0(i) = rho_left
      else if (s - 0.5d0*l <= x .and. x <= s + 0.5d0*l) then
        rho_field % rho0(i) = 0.5d0*rho_left * ( (1 + rho_right/rho_left) - (1 - rho_right/rho_left)*sin(dpi*(x - s)/l) )
      else
        rho_field % rho0(i) = rho_right
      end if
    end do

  end subroutine resonant_absorption_eq

end submodule smod_equil_resonant_absorption
