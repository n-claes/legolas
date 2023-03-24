module mod_physics
  use mod_physics_utils, only: physics_i
  use mod_resistivity, only: resistivity_t, new_resistivity
  use mod_gravity, only: gravity_t, new_gravity
  implicit none

  private

  type, public :: physics_t
    type(resistivity_t) :: resistivity
    type(gravity_t) :: gravity
  contains
    procedure, public :: set_resistivity_funcs
    procedure, public :: set_gravity_funcs
    procedure, public :: delete
  end type physics_t

  public :: new_physics

contains

  function new_physics() result(physics)
    type(physics_t) :: physics
    physics%resistivity = new_resistivity()
    physics%gravity = new_gravity()
  end function new_physics


  subroutine set_resistivity_funcs(this, eta_func, detadT_func, detadr_func)
    class(physics_t), intent(inout) :: this
    procedure(physics_i) :: eta_func
    procedure(physics_i), optional :: detadT_func
    procedure(physics_i), optional :: detadr_func

    this%resistivity%eta => eta_func
    if (present(detadT_func)) this%resistivity%detadT => detadT_func
    if (present(detadr_func)) this%resistivity%detadr => detadr_func
  end subroutine set_resistivity_funcs


  subroutine set_gravity_funcs(this, g0_func)
    class(physics_t), intent(inout) :: this
    procedure(physics_i) :: g0_func
    this%gravity%g0 => g0_func
  end subroutine set_gravity_funcs


  pure subroutine delete(this)
    class(physics_t), intent(inout) :: this
    call this%resistivity%delete()
    call this%gravity%delete()
  end subroutine delete

end module mod_physics
