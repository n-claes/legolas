module mod_physics
  use mod_physics_utils, only: physics_i
  use mod_resistivity, only: resistivity_t, new_resistivity
  implicit none

  private

  type, public :: physics_t
    type(resistivity_t) :: resistivity
  contains
    procedure, public :: set_resistivity_funcs
    procedure, public :: delete
  end type physics_t

  public :: new_physics

contains

  function new_physics() result(physics)
    type(physics_t) :: physics
    physics%resistivity = new_resistivity()
  end function new_physics


  subroutine set_resistivity_funcs(this, eta_func, detadT_func, detadr_func)
    class(physics_t), intent(inout) :: this
    procedure(physics_i), optional :: eta_func
    procedure(physics_i), optional :: detadT_func
    procedure(physics_i), optional :: detadr_func

    if (present(eta_func)) this%resistivity%eta => eta_func
    if (present(detadT_func)) this%resistivity%detadT => detadT_func
    if (present(detadr_func)) this%resistivity%detadr => detadr_func
  end subroutine set_resistivity_funcs


  pure subroutine delete(this)
    class(physics_t), intent(inout) :: this
    call this%resistivity%delete()
  end subroutine delete

end module mod_physics
