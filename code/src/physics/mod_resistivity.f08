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
  use mod_global_variables
  use mod_physical_constants
  implicit none

  public

  private :: get_constants


contains

  !> Calculates the Spitzer resistivity based on the equilibrium temperatures.
  !! @param[in]  T0   Array with the equilibrium temperatures, in K
  !! @param[out] eta  Array containing the Spitzer resistivity as a function of T
  subroutine get_eta(T0, eta)
    real(dp), intent(in)    :: T0(4*gridpts)
    real(dp), intent(out)   :: eta(4*gridpts)

    real(dp)                :: ec, me, e0, kB, eta_1MK

    call get_constants(ec, me, e0, kB)

    eta = dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
          ((4 * dpi * e0)**2 * (kB * T0)**1.5d0)

    !! Set the unit resistivity, such that the normalised resistivity
    !! at 1 MK equals approximately 0.1. This can be done, as a unit current
    !! can be chosen at random.
    eta_1MK = dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
              ((4 * dpi * e0)**2 * (kB * 1.0d6)**1.5d0)
    call set_unit_resistivity(eta_1MK / 0.1d0)

  end subroutine get_eta

  !> Calculates the derivative of the Spitzer resistivity with respect to
  !! temperature.
  !! @param[in]  T0      Array with the equilibrium temperatures, in K
  !! @param[out] deta_dT Array containing the derivative of the Spitzer
  !!                     resistivity with respect to temperature
  subroutine get_deta_dT(T0, deta_dT)
    real(dp), intent(in)    :: T0(4*gridpts)
    real(dp), intent(out)   :: deta_dT(4*gridpts)

    real(dp)                :: ec, me, e0, kB

    call get_constants(ec, me, e0, kB)

    deta_dT = -1.5d0 * dpi * Z_ion * ec**2 * sqrt(me) * log(coulomb_log) / &
              ((4 * dpi * e0)**2 * kB**1.5d0 * T0**2.5d0)

  end subroutine get_deta_dT

  !> Gets the physical constants relevant for the Spitzer resistivity.
  !! Retrieves the values in SI or cgs units.
  !! @param[out] ec   The electric charge in SI or cgs units
  !! @param[out] me   The electron mass in SI or cgs units
  !! @param[out] e0   The permittivity of free space in SI or cgs units
  !! @param[out] kB   The Boltzmann constant in SI or cgs units
  subroutine get_constants(ec, me, e0, kB)
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
