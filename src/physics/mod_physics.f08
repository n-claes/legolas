module mod_physics
  use mod_resistivity, only: resistivity_t, new_resistivity
  implicit none

  private

  type, public :: physics_t
    type(resistivity_t) :: resistivity
  contains
    procedure, public :: delete
  end type physics_t

  public :: new_physics

contains

  function new_physics() result(physics)
    type(physics_t) :: physics
    physics%resistivity = new_resistivity()
  end function new_physics


  pure subroutine delete(this)
    class(physics_t), intent(inout) :: this
    call this%resistivity%delete()
  end subroutine delete

end module mod_physics
