!
! MODULE: mod_gravity
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to set the external gravity.
!
module mod_gravity
  use mod_global_variables
  implicit none

  public

  !> Gravitational acceleration
  real(dp), allocatable       :: grav_eq(:)

contains

  !> Initialises the gravitational acceleration parameter grav, set to
  !! zero if gravity is not included.
  subroutine initialise_gravity()
    allocate(grav_eq(gauss_gridpts))

    grav_eq = 0.0d0
  end subroutine initialise_gravity

end module mod_gravity
