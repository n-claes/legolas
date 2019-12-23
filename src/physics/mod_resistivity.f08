!
! MODULE: mod_resistivity
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to calculate the resistivity contributions.
!
module mod_resistivity
  use mod_global_variables, only: dp, gauss_gridpts
  use mod_physical_constants, only: dpi, Z_ion, coulomb_log
  implicit none

  private
  
  public :: get_eta
  public :: get_deta_dT

contains

  !> Calculates the Spitzer resistivity based on the equilibrium temperatures.
  !! @param[in]  T0   Array with the equilibrium temperatures, in K
  !! @param[out] eta  Array containing the Spitzer resistivity as a function of T
  subroutine get_eta(T0, eta)
    use mod_global_variables, only: use_fixed_resistivity, fixed_eta_value
    use mod_physical_constants, only: unit_temperature, set_unit_resistivity, unit_resistivity
    
    real(dp), intent(in)    :: T0(gauss_gridpts)
    real(dp), intent(out)   :: eta(gauss_gridpts)

    real(dp)                :: ec, me, e0, kB, eta_1MK
    real(dp)                :: T0_denorm(gauss_gridpts)

    if (use_fixed_resistivity) then
      eta(:) = fixed_eta_value
      return
    end if

    !! Denormalise variables for calculation
    T0_denorm = T0 * unit_temperature

    call get_constants(ec, me, e0, kB)

    eta = dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
          ((4 * dpi * e0)**2 * (kB * T0_denorm)**1.5d0)

    !! Set the unit resistivity, such that the normalised resistivity
    !! at 1 MK equals approximately 0.1. This can be done, as a unit current
    !! can be chosen at random.
    eta_1MK = dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
              ((4 * dpi * e0)**2 * (kB * 1.0d6)**1.5d0)
    call set_unit_resistivity(eta_1MK / 0.1d0)

    !! Renormalise
    eta = eta / unit_resistivity

  end subroutine get_eta

  !> Calculates the derivative of the Spitzer resistivity with respect to
  !! temperature.
  !! @param[in]  T0      Array with the equilibrium temperatures, in K
  !! @param[out] deta_dT Array containing the derivative of the Spitzer
  !!                     resistivity with respect to temperature
  subroutine get_deta_dT(T0, deta_dT)
    use mod_global_variables, only: use_fixed_resistivity
    use mod_physical_constants, only: unit_temperature, unit_deta_dT
    
    real(dp), intent(in)    :: T0(gauss_gridpts)
    real(dp), intent(out)   :: deta_dT(gauss_gridpts)

    real(dp)                :: ec, me, e0, kB
    real(dp)                :: T0_denorm(gauss_gridpts)

    if (use_fixed_resistivity) then
      deta_dT(:) = 0.0d0
    end if

    !! Denormalise variables for calculation
    T0_denorm = T0 * unit_temperature

    call get_constants(ec, me, e0, kB)

    deta_dT = -1.5d0 * dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
              ((4 * dpi * e0)**2 * kB**1.5d0 * T0_denorm**2.5d0)
    !! Renormalise
    deta_dT = deta_dT / unit_deta_dT

  end subroutine get_deta_dT

  !> Gets the physical constants relevant for the Spitzer resistivity.
  !! Retrieves the values in SI or cgs units.
  !! @param[out] ec   The electric charge in SI or cgs units
  !! @param[out] me   The electron mass in SI or cgs units
  !! @param[out] e0   The permittivity of free space in SI or cgs units
  !! @param[out] kB   The Boltzmann constant in SI or cgs units
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
