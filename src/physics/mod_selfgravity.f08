module mod_selfgravity
  use mod_global_variables, only: dp
  implicit none

  private

  public :: get_gravity_prefactor

contains

  real(dp) function get_gravity_prefactor()
    use mod_physical_constants, only: dpi, bigG_cgs
    use mod_units, only: unit_length, unit_mass, unit_time

    get_gravity_prefactor = ( &
      4.0d0 * dpi * bigG_cgs / (unit_length**3 / (unit_mass * unit_time**2 )) &
    )
  end function get_gravity_prefactor

end module mod_selfgravity