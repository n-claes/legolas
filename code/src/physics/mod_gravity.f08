module mod_gravity
  use mod_global_variables
  implicit none

  !> Gravitational acceleration
  real(dp)                     :: g

  !> Surface gravity of Earth in m/s^2 (SI)
  real(dp), parameter, private :: g_earth_si = 9.81d0
  !> Surface gravity of Earth in cm/s^2 (cgs)
  real(dp), parameter, private :: g_earth_cgs = 9.81d3
  !> Surface gravity of the Sun in m/s^2 (SI)
  real(dp), parameter, private :: g_solar_si = 274d0
  !> Surface gravity of the Sun in cm/s^2 (cgs)
  real(dp), parameter, private :: g_solar_cgs = 274d3


contains

  subroutine initialise_gravity()

    ! Obtain gravity type
    select case(gravity_type)

    case('earth')
      write(*, *) "Using Earth-like external gravity."
      if (cgs_units) then
        g = g_earth_cgs
      else
        g = g_earth_si
      end if

    case('solar')
      write(*, *) "Using solar-like external gravity."
      if (cgs_units) then
        g = g_solar_cgs
      else
        g = g_solar_si
      end if

    case default
      write(*, *) "Gravity type not defined correctly."
      write(*, *) "Currently set on:   ", gravity_type
      stop
    end select

  end subroutine initialise_gravity






end module mod_gravity
