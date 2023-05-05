module mod_conduction_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: conduction_settings_t
    logical, private :: para_conduction
    logical, private :: fixed_tc_para
    real(dp), private :: fixed_tc_para_value
    logical, private :: perp_conduction
    logical, private :: fixed_tc_perp
    real(dp), private :: fixed_tc_perp_value

  contains

    procedure, public :: disable
    procedure, public :: is_enabled
    procedure, public :: enable_para_conduction
    procedure, public :: has_parallel_conduction
    procedure, public :: set_fixed_tc_para
    procedure, public :: get_fixed_tc_para
    procedure, public :: has_fixed_tc_para
    procedure, public :: enable_perp_conduction
    procedure, public :: has_perpendicular_conduction
    procedure, public :: set_fixed_tc_perp
    procedure, public :: get_fixed_tc_perp
    procedure, public :: has_fixed_tc_perp
  end type conduction_settings_t

  public :: new_conduction_settings

contains

  pure function new_conduction_settings() result(conduction)
    type(conduction_settings_t) :: conduction

    conduction%para_conduction = .false.
    conduction%fixed_tc_para = .false.
    conduction%fixed_tc_para_value = 0.0_dp
    conduction%perp_conduction = .false.
    conduction%fixed_tc_perp = .false.
    conduction%fixed_tc_perp_value = 0.0_dp
  end function new_conduction_settings


  pure logical function is_enabled(this)
    class(conduction_settings_t), intent(in) :: this
    is_enabled = this%has_parallel_conduction() .or. this%has_perpendicular_conduction()
  end function is_enabled


  pure subroutine disable(this)
    class(conduction_settings_t), intent(inout) :: this
    this%para_conduction = .false.
    this%fixed_tc_para = .false.
    this%perp_conduction = .false.
    this%fixed_tc_perp = .false.
  end subroutine disable


  pure subroutine enable_para_conduction(this)
    class(conduction_settings_t), intent(inout) :: this
    this%para_conduction = .true.
  end subroutine enable_para_conduction


  pure logical function has_parallel_conduction(this)
    class(conduction_settings_t), intent(in) :: this
    has_parallel_conduction = this%para_conduction
  end function has_parallel_conduction


  pure subroutine set_fixed_tc_para(this, tc_para)
    class(conduction_settings_t), intent(inout) :: this
    real(dp), intent(in) :: tc_para
    this%fixed_tc_para_value = tc_para
    this%fixed_tc_para = .true.
    this%para_conduction = .true.
  end subroutine set_fixed_tc_para


  pure real(dp) function get_fixed_tc_para(this)
    class(conduction_settings_t), intent(in) :: this
    get_fixed_tc_para = this%fixed_tc_para_value
  end function get_fixed_tc_para


  pure logical function has_fixed_tc_para(this)
    class(conduction_settings_t), intent(in) :: this
    has_fixed_tc_para = this%is_enabled() .and. this%fixed_tc_para
  end function has_fixed_tc_para


  pure subroutine enable_perp_conduction(this)
    class(conduction_settings_t), intent(inout) :: this
    this%perp_conduction = .true.
  end subroutine enable_perp_conduction


  pure logical function has_perpendicular_conduction(this)
    class(conduction_settings_t), intent(in) :: this
    has_perpendicular_conduction = this%perp_conduction
  end function has_perpendicular_conduction


  pure subroutine set_fixed_tc_perp(this, tc_perp)
    class(conduction_settings_t), intent(inout) :: this
    real(dp), intent(in) :: tc_perp
    this%fixed_tc_perp_value = tc_perp
    this%fixed_tc_perp = .true.
    this%perp_conduction = .true.
  end subroutine set_fixed_tc_perp


  pure real(dp) function get_fixed_tc_perp(this)
    class(conduction_settings_t), intent(in) :: this
    get_fixed_tc_perp = this%fixed_tc_perp_value
  end function get_fixed_tc_perp


  pure logical function has_fixed_tc_perp(this)
    class(conduction_settings_t), intent(in) :: this
    has_fixed_tc_perp = this%is_enabled() .and. this%fixed_tc_perp
  end function has_fixed_tc_perp

end module mod_conduction_settings
