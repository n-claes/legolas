! =============================================================================
!> This submodule defines an inhomogeneous medium in Cartesian geometry with
!! a constant resistivity value. Two (constant) density profiles are defined
!! which are connected by a sine profile, representing the interface inbetween.
!! This density profile allows for resonant absorption, the geometry
!! can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Van Doorsselaere, T., & Poedts, S. (2007).
!!  Modifications to the resistive MHD spectrum due to changes in the equilibrium.
!!  Plasma Physics and Controlled Fusion, 49(3), 261_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 0.05
!! - <tt>p1</tt> = 0.9 : used to set the left density value.
!! - <tt>p2</tt> = 0.1 : used to set the right density value.
!! - <tt>r0</tt> = 0.2 : used to set the width of the interface.
!! - <tt>cte_B02</tt> = 0 : used to set the By value.
!! - <tt>cte_B03</tt> = 1 : used to set the Bz value.
!! - <tt>cte_T0</tt> = 0 : used to set the temperature. Zero by default.
!! - fixed eta value of \(10^{-3.2}\)
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_resonant_absorption
  implicit none

contains

  !> Sets the equilibrium
  module subroutine resonant_absorption_eq()
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_equilibrium_params, only: p1, p2, r0, cte_T0, cte_B02, cte_B03

    real(dp)  :: x, s, r0, rho_left, rho_right, zeta
    integer   :: i

    call allow_geometry_override(default_geometry='Cartesian', default_x_start=0.0d0, default_x_end=1.0d0)
    call initialise_grid()

    if (use_defaults) then
      resistivity = .true.
      use_fixed_resistivity = .true.
      fixed_eta_value = 10.0d0**(-3.2d0)

      k2 = 1.0d0
      k3 = 0.05d0

      rho_left = 0.9d0
      rho_right = 0.1d0
      r0 = 0.2d0
      cte_T0 = 0.0d0
      cte_B02 = 0.0d0
      cte_B03 = 1.0d0
    else
      rho_left = p1
      rho_right = p2
    end if

    s = 0.5d0 * (x_start + x_end)
    zeta = rho_left / rho_right

    B_field % B02 = cte_B02
    B_field % B03 = cte_B03
    B_field % B0 = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    T_field % T0 = cte_T0

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      if (x >= x_start .and. x < s - 0.5d0*r0) then
        rho_field % rho0(i) = rho_left
      else if (s - 0.5d0*r0 <= x .and. x <= s + 0.5d0*r0) then
        rho_field % rho0(i) = 0.5d0 * rho_left * (1.0d0 + 1.0d0 / zeta - (1.0d0 - 1.0d0 / zeta) * sin(dpi * (x - s) / r0))
        rho_field % d_rho0_dr(i) = dpi * rho_left * (1.0d0 / zeta - 1.0d0) * cos(dpi * (x - s) / r0) / (2.0d0 * r0)
      else
        rho_field % rho0(i) = rho_right
      end if
    end do
  end subroutine resonant_absorption_eq

end submodule smod_equil_resonant_absorption
