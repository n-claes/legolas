! =============================================================================
!> Module containing resistivity-related routines, calculates
!! and sets the resistivity values based on the equilibrium configuration.
module mod_resistivity
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi, Z_ion, coulomb_log
  use mod_settings, only: settings_t
  implicit none

  private

  public :: set_resistivity_values

contains


  !> This routines sets all resistivity values in \p eta_field,
  !! and calls all other relevant subroutines defined in this module.
  subroutine set_resistivity_values(settings, T_field, eta_field)
    use mod_types, only: temperature_type, resistivity_type

    type(settings_t), intent(in) :: settings
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)    :: T_field
    !> the type containing the resistivity attributes
    type(resistivity_type), intent(inout) :: eta_field

    call get_eta(settings, T_field % T0, eta_field % eta)
    call get_deta_dT(settings, T_field % T0, eta_field % d_eta_dT)
    if (settings%physics%resistivity%use_dropoff) then
      call set_eta_dropoff(settings, eta_field)
    end if
  end subroutine set_resistivity_values


  !> Calculates the resistivity. Returns either the full Spitzer resistivity
  !! based on the equilibrium parameters, or a fixed resistivity value if
  !! specified in the global variables module. The unit resistivity is also set in this routine.
  !! If a fixed resistivity is used, eta is assumed to be normalised and
  !! the unit resistivity remains unity.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the resistivity is normalised on exit.
  subroutine get_eta(settings, T0_eq, eta)
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in) :: T0_eq(:)
    !> resistivity values, normalised on exit
    real(dp), intent(inout) :: eta(:)

    real(dp) :: ec, me, kB
    real(dp) :: T0_dimfull(size(T0_eq))

    if (settings%physics%resistivity%has_fixed_resistivity()) then
      eta = settings%physics%resistivity%get_fixed_resistivity()
      return
    end if

    call get_constants(ec, me, kB)
    T0_dimfull = T0_eq * settings%units%get_unit_temperature()
    eta = ( &
      (4.0d0 / 3.0d0) * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
      / (kB * T0_dimfull)**(3.0d0 / 2.0d0) &
    ) / settings%units%get_unit_resistivity()
  end subroutine get_eta


  !> Calculates the derivative of the resistivity.
  !! Returns the derivative of the full Spitzer resistivity with respect
  !! to the equilibrium temperature, if a fixed resistivity was set instead
  !! this routine returns zero.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the derivative of the resistivity is normalised on exit.
  subroutine get_deta_dT(settings, T0_eq, deta_dT)
    type(settings_t), intent(in) :: settings
    !> equilibrium temperature
    real(dp), intent(in) :: T0_eq(:)
    !> derivative of resistivity with respect to temperature, normalised on exit
    real(dp), intent(out) :: deta_dT(:)

    real(dp) :: ec, me, kB
    real(dp) :: unit_temperature, unit_deta_dT
    real(dp) :: T0_dimfull(size(T0_eq))

    if (settings%physics%resistivity%has_fixed_resistivity()) then
      deta_dT = 0.0d0
      return
    end if

    call get_constants(ec, me, kB)
    unit_temperature = settings%units%get_unit_temperature()
    T0_dimfull = T0_eq * unit_temperature
    unit_deta_dT = settings%units%get_unit_resistivity() / unit_temperature
    deta_dT = ( &
      -2.0d0 * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
      / (kB**(3.0d0 / 2.0d0) * T0_dimfull**(5.0d0 / 2.0d0)) &
    ) / unit_deta_dT
  end subroutine get_deta_dT


  !> Sets a hyperbolic tangent profile for the resistivity so it goes smoothly to zero
  !! near the edges. The location and width of the dropoff profile can be controlled
  !! through <tt>dropoff_edge_dist</tt> and <tt>dropoff_width</tt>.
  !! @warning Throws an error if:
  !!
  !! - a dropoff profile is requested but the resistivity is not constant. @endwarning
  subroutine set_eta_dropoff(settings, eta_field)
    use mod_types, only: resistivity_type
    use mod_logging, only: logger, str
    use mod_grid, only: grid_gauss
    use mod_physical_constants, only: dpi

    type(settings_t), intent(in) :: settings
    !> the type containing the resistivity attributes
    type(resistivity_type), intent(inout) :: eta_field

    real(dp) :: sleft, sright, width, etaval, edge_dist
    real(dp) :: x, shift, stretch
    integer :: i, gauss_gridpts

    if (.not. settings%physics%resistivity%has_fixed_resistivity()) then
      call logger%error("eta dropoff only possible with a fixed resistivity value")
      return
    end if

    gauss_gridpts = settings%grid%get_gauss_gridpts()
    width = settings%physics%dropoff_width
    edge_dist = settings%physics%dropoff_edge_dist
    etaval = settings%physics%resistivity%get_fixed_resistivity()
    sleft = grid_gauss(1) + 0.5d0 * width + edge_dist
    sright = grid_gauss(gauss_gridpts) - edge_dist - 0.5d0 * width

    shift = etaval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
    stretch = etaval / (tanh(dpi) - tanh(-dpi))

    call logger%info("setting eta-dropoff profile")
    call logger%info("dropoff width = " // str(width))
    call logger%info("distance from edge = " // str(edge_dist))

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
  !! Returns all physical constants used to calculate the Spitzer resistivity.
  subroutine get_constants(ec, me, kB)
    use mod_physical_constants, only: ec_cgs, me_cgs, kB_cgs

    !> electric charge
    real(dp), intent(out) :: ec
    !> electron mass
    real(dp), intent(out) :: me
    !> Boltzmann constant
    real(dp), intent(out) :: kB

    ec = ec_cgs
    me = me_cgs
    kB = kB_cgs
  end subroutine get_constants

end module mod_resistivity
