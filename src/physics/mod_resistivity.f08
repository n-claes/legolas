! =============================================================================
!> @brief   Module containing resistivity-related routines.
!! @details This module is responsible for calculating and setting the resistivity
!!          values based on the equilibrium configuration.
module mod_resistivity
  use mod_global_variables, only: dp, gauss_gridpts
  use mod_physical_constants, only: dpi, Z_ion, coulomb_log
  implicit none

  private

  public :: set_resistivity_values

contains


  !> @brief   Main routine to set the resistivity.
  !! @brief   This routines sets all resistivity values in \p eta_field,
  !!          and calls all other relevant subroutines defined in this module.
  !! @param[in] T_field   the type containing the temperature attributes
  !! @param[in, out] eta_field  the type containing the resistivity attributes
  subroutine set_resistivity_values(T_field, eta_field)
    use mod_types, only: temperature_type, resistivity_type

    type(temperature_type), intent(in)    :: T_field
    type(resistivity_type), intent(inout) :: eta_field

    call get_eta(T_field % T0, eta_field % eta)
    call get_deta_dT(T_field % T0, eta_field % d_eta_dT)
  end subroutine set_resistivity_values


  !> @brief   Calculates the resistivity.
  !! @details Returns either the full Spitzer resistivity based on the equilibrium
  !!          parameters, or a fixed resistivity value if specified in the global
  !!          variables module. The unit resistivity is also set in this routine.
  !!          If a fixed resistivity is used, eta is assumed to be normalised and
  !!          the unit resistivity remains unity.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the resistivity is normalised on exit.
  !! @param[in]  T0_eq  equilibrium temperature
  !! @param[out] eta    resistivity values
  subroutine get_eta(T0_eq, eta)
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_units, only: unit_temperature, set_unit_resistivity, unit_resistivity

    real(dp), intent(in)    :: T0_eq(gauss_gridpts)
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
    !! Set the unit resistivity, such that the normalised resistivity
    !! at 1 MK equals approximately 0.1. This can be done, as a unit current
    !! can be chosen at random.
    eta_1MK = (4.0d0 / 3.0d0) * sqrt(2.0d0 * dpi) * Z_ion * ec**2 * sqrt(me) * coulomb_log &
          / ((4 * dpi * e0)**2 * (kB * 1.0d6)**(3.0d0 / 2.0d0))
    call set_unit_resistivity(eta_1MK / 0.1d0)
    !! Renormalise
    eta = eta / unit_resistivity
  end subroutine get_eta


  !> @brief   Calculates the derivative of the resistivity.
  !! @details Returns the derivative of the full Spitzer resistivity with respect
  !!          to the equilibrium temperature. If a fixed resistivity was set instead
  !!          this routine returns zero.
  !! @note    The temperature should be normalised when calling this routine,
  !!          the derivative of the resistivity is normalised on exit.
  !! @param[in]  T0_eq    normalised equilibrium temperature
  !! @param[out] deta_dT  normalised derivative of resistivity with respect to temperature
  subroutine get_deta_dT(T0_eq, deta_dT)
    use mod_global_variables, only: use_fixed_resistivity
    use mod_units, only: unit_temperature, unit_deta_dT

    real(dp), intent(in)    :: T0_eq(gauss_gridpts)
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


  !> @brief   Retrieves resistivity constants.
  !! @details Returns all physical constants used to calculate the Spitzer resistivity,
  !!          either in cgs (default) or SI depending on the unit system chosen.
  !! @param[out] ec   electric charge
  !! @param[out] me   electron mass
  !! @param[out] e0   permittivity of free space
  !! @param[out] kB   Boltzmann constant
  subroutine get_constants(ec, me, e0, kB)
    use mod_global_variables, only: cgs_units
    use mod_physical_constants, only: ec_cgs, me_cgs, kB_cgs, ec_si, me_si, e0_si, kB_si

    real(dp), intent(out) :: ec, me, e0, kB

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
