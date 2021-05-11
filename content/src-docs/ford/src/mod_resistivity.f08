! =============================================================================
!> Module containing resistivity-related routines, calculates
!! and sets the resistivity values based on the equilibrium configuration.
module mod_resistivity
  use mod_global_variables, only: dp, gauss_gridpts
  use mod_physical_constants, only: dpi, Z_ion, coulomb_log
  implicit none

  private

  public :: set_resistivity_values

contains


  !> This routines sets all resistivity values in \p eta_field,
  !! and calls all other relevant subroutines defined in this module.
  subroutine set_resistivity_values(T_field, eta_field)
    use mod_types, only: temperature_type, resistivity_type
    use mod_global_variables, only: use_eta_dropoff

    !> the type containing the temperature attributes
    type(temperature_type), intent(in)    :: T_field
    !> the type containing the resistivity attributes
    type(resistivity_type), intent(inout) :: eta_field

    call get_eta(T_field % T0, eta_field % eta)
    call get_deta_dT(T_field % T0, eta_field % d_eta_dT)
    if (use_eta_dropoff) then
      call set_eta_dropoff(eta_field)
    end if
  end subroutine set_resistivity_values


  !> Calculates the resistivity. Returns either the full Spitzer resistivity
  !! based on the equilibrium parameters, or a fixed resistivity value if
  !! specified in the global variables module. The unit resistivity is also set in this routine.
  !! If a fixed resistivity is used, eta is assumed to be normalised and
  !! the unit resistivity remains unity.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the resistivity is normalised on exit.
  subroutine get_eta(T0_eq, eta)
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_units, only: unit_temperature, set_unit_resistivity, unit_resistivity

    !> equilibrium temperature
    real(dp), intent(in)    :: T0_eq(gauss_gridpts)
    !> resistivity values, normalised on exit
    real(dp), intent(inout) :: eta(gauss_gridpts)

    real(dp)                :: ec, me, e0, kB, eta_1MK
    real(dp)                :: T0_dimfull(gauss_gridpts)

    if (use_fixed_resistivity) then
      eta = fixed_eta_value
      return
    end if

    call get_constants(ec, me, e0, kB)
    T0_dimfull = T0_eq * unit_temperature
    eta = (4.0d0 / 3.0d0) * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
          / ((4 * dpi * e0)**2 * (kB * T0_dimfull)**(3.0d0 / 2.0d0))
    ! Set the unit resistivity, such that the normalised resistivity
    ! at 1 MK equals approximately 0.1. This can be done, as a unit current
    ! can be chosen at random.
    eta_1MK = (4.0d0 / 3.0d0) * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
          / ((4 * dpi * e0)**2 * (kB * 1.0d6)**(3.0d0 / 2.0d0))
    call set_unit_resistivity(eta_1MK / 0.1d0)
    ! Renormalise
    eta = eta / unit_resistivity
  end subroutine get_eta


  !> Calculates the derivative of the resistivity.
  !! Returns the derivative of the full Spitzer resistivity with respect
  !! to the equilibrium temperature, if a fixed resistivity was set instead
  !! this routine returns zero.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the derivative of the resistivity is normalised on exit.
  subroutine get_deta_dT(T0_eq, deta_dT)
    use mod_global_variables, only: use_fixed_resistivity
    use mod_units, only: unit_temperature, unit_deta_dT

    !> equilibrium temperature
    real(dp), intent(in)    :: T0_eq(gauss_gridpts)
    !> derivative of resistivity with respect to temperature, normalised on exit
    real(dp), intent(out)   :: deta_dT(gauss_gridpts)

    real(dp)                :: ec, me, e0, kB
    real(dp)                :: T0_dimfull(gauss_gridpts)

    if (use_fixed_resistivity) then
      deta_dT = 0.0d0
      return
    end if

    call get_constants(ec, me, e0, kB)
    T0_dimfull = T0_eq * unit_temperature
    deta_dT = -2.0d0 * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
              / ((4 * dpi * e0)**2 * kB**(3.0d0 / 2.0d0) * T0_dimfull**(5.0d0 / 2.0d0))
    deta_dT = deta_dT / unit_deta_dT
  end subroutine get_deta_dT


  !> Sets a hyperbolic tangent profile for the resistivity so it goes smoothly to zero
  !! near the edges. The location and width of the dropoff profile can be controlled
  !! through <tt>dropoff_edge_dist</tt> and <tt>dropoff_width</tt>.
  !! @warning Throws an error if:
  !!
  !! - a dropoff profile is requested but the resistivity is not constant. @endwarning
  subroutine set_eta_dropoff(eta_field)
    use mod_types, only: resistivity_type
    use mod_logging, only: log_message, char_log
    use mod_grid, only: grid_gauss
    use mod_physical_constants, only: dpi
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value, &
                                    dropoff_edge_dist, dropoff_width

    !> the type containing the resistivity attributes
    type(resistivity_type), intent(inout) :: eta_field

    real(dp)  :: sleft, sright, width, etaval
    real(dp)  :: x, shift, stretch
    integer   :: i

    if (.not. use_fixed_resistivity) then
      call log_message('eta dropoff only possible with a fixed resistivity value', level='error')
    end if

    width = dropoff_width
    etaval = fixed_eta_value
    sleft = grid_gauss(1) + 0.5d0 * width + dropoff_edge_dist
    sright = grid_gauss(gauss_gridpts) - dropoff_edge_dist - 0.5d0 * width

    shift = etaval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
    stretch = etaval / (tanh(dpi) - tanh(-dpi))

    call log_message('setting eta-dropoff profile', level='info')
    write(char_log, '(f20.2)') dropoff_width
    call log_message('dropoff width = ' // adjustl(trim(char_log)), level='info')
    write(char_log, '(f20.2)') dropoff_edge_dist
    call log_message('distance from edge = ' // adjustl(trim(char_log)), level='info')

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      if (sleft - 0.5d0*width <= x .and. x <= sleft + 0.5d0*width) then
        eta_field % eta(i) = shift + stretch * tanh(2.0d0 * dpi * (x - sleft) / width)
        eta_field % d_eta_dr(i) = (2.0d0 * dpi * stretch / width) / cosh(2.0d0 * dpi * (x - sleft) / width)**2
      else if (sleft + 0.5d0*width < x .and. x < sright - 0.5d0*width) then
        eta_field % eta(i) = etaval
        eta_field % d_eta_dT(i) = 0.0d0
      else if (sright - 0.5d0*width <= x .and. x <= sright + 0.5d0*width) then
        eta_field % eta(i) = shift + stretch * tanh(2.0d0 * dpi * (sright - x) / width)
        eta_field % d_eta_dr(i) = (-2.0d0 * dpi * stretch / width) / cosh(2.0d0 * dpi * (sright - x) / width)**2
      else
        eta_field % eta(i) = 0.0d0
        eta_field % d_eta_dr(i) = 0.0d0
      end if
    end do
  end subroutine set_eta_dropoff


  !> Retrieves resistivity constants.
  !! Returns all physical constants used to calculate the Spitzer resistivity,
  !! either in cgs (default) or SI depending on the unit system chosen.
  subroutine get_constants(ec, me, e0, kB)
    use mod_global_variables, only: cgs_units
    use mod_physical_constants, only: ec_cgs, me_cgs, kB_cgs, ec_si, me_si, e0_si, kB_si

    !> electric charge
    real(dp), intent(out) :: ec
    !> electron mass
    real(dp), intent(out) :: me
    !> permittivity of free space
    real(dp), intent(out) :: e0
    !> Boltzmann constant
    real(dp), intent(out) :: kB

    if (cgs_units) then
      ec = ec_cgs
      me = me_cgs
      e0 = 1.0d0
      kB = kB_cgs
    else
      ec = ec_si
      me = me_si
      e0 = e0_si
      kB = kB_si
    end if
  end subroutine get_constants

end module mod_resistivity
