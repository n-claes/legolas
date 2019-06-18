module mod_gravity
  use mod_global_variables
  implicit none

  private

  !> Gravitational acceleration
  real(dp)            :: g

  !> Surface gravity of Earth in m/s^2 (SI)
  real(dp), parameter :: g_earth_si = 9.81d0
  !> Surface gravity of Earth in cm/s^2 (cgs)
  real(dp), parameter :: g_earth_cgs = 9.81d3
  !> Surface gravity of the Sun in m/s^2 (SI)
  real(dp), parameter :: g_solar_si = 274d0
  !> Surface gravity of the Sun in cm/s^2 (cgs)
  real(dp), parameter :: g_solar_cgs = 274d3

  public  :: g
  public  :: initialise_gravity

contains

  subroutine initialise_gravity()
    if (external_gravity) then
      call set_gravity()
    else
      g = 0.0d0
    end if
  end subroutine initialise_gravity

  subroutine set_gravity()
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
  end subroutine set_gravity

end module mod_gravity
