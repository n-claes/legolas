module mod_physics
  use mod_physics_utils, only: physics_i
  use mod_resistivity, only: resistivity_t, new_resistivity
  use mod_gravity, only: gravity_t, new_gravity
  use mod_hall, only: hall_t, new_hall
  use mod_thermal_conduction, only: conduction_t, new_conduction
  use mod_radiative_cooling, only: cooling_t, new_cooling
  implicit none

  private

  type, public :: physics_t
    type(resistivity_t) :: resistivity
    type(gravity_t) :: gravity
    type(hall_t) :: hall
    type(conduction_t) :: conduction
    type(cooling_t) :: cooling
  contains
    procedure, public :: set_resistivity_funcs
    procedure, public :: set_gravity_funcs
    procedure, public :: set_parallel_conduction_funcs
    procedure, public :: set_perpendicular_conduction_funcs
    procedure, public :: delete
  end type physics_t

  public :: new_physics

contains

  function new_physics() result(physics)
    type(physics_t) :: physics
    physics%resistivity = new_resistivity()
    physics%gravity = new_gravity()
    physics%hall = new_hall()
    physics%conduction = new_conduction()
    physics%cooling = new_cooling()
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


  subroutine set_parallel_conduction_funcs( &
    this, tcpara_func, dtcparadT_func, dtcparadr_func &
  )
    class(physics_t), intent(inout) :: this
    procedure(physics_i) :: tcpara_func
    procedure(physics_i), optional :: dtcparadT_func
    procedure(physics_i), optional :: dtcparadr_func

    this%conduction%tcpara => tcpara_func
    if (present(dtcparadT_func)) this%conduction%dtcparadT => dtcparadT_func
    if (present(dtcparadr_func)) this%conduction%dtcparadr => dtcparadr_func
  end subroutine set_parallel_conduction_funcs


  subroutine set_perpendicular_conduction_funcs( &
    this, &
    tcperp_func, &
    dtcperpdT_func, &
    dtcperpdrho_func, &
    dtcperpdB2_func, &
    dtcperpdr_func &
  )
    class(physics_t), intent(inout) :: this
    procedure(physics_i) :: tcperp_func
    procedure(physics_i), optional :: dtcperpdT_func
    procedure(physics_i), optional :: dtcperpdrho_func
    procedure(physics_i), optional :: dtcperpdB2_func
    procedure(physics_i), optional :: dtcperpdr_func

    this%conduction%tcperp => tcperp_func
    if (present(dtcperpdT_func)) this%conduction%dtcperpdT => dtcperpdT_func
    if (present(dtcperpdrho_func)) this%conduction%dtcperpdrho => dtcperpdrho_func
    if (present(dtcperpdB2_func)) this%conduction%dtcperpdB2 => dtcperpdB2_func
    if (present(dtcperpdr_func)) this%conduction%dtcperpdr => dtcperpdr_func
  end subroutine set_perpendicular_conduction_funcs


  pure subroutine delete(this)
    class(physics_t), intent(inout) :: this
    call this%resistivity%delete()
    call this%gravity%delete()
    call this%hall%delete()
    call this%conduction%delete()
    call this%cooling%delete()
  end subroutine delete

end module mod_physics
