module mod_bg_magnetic
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: bg_magnetic_t
    procedure(real(dp)), pointer, nopass :: B01
    procedure(real(dp)), pointer, nopass :: B02
    procedure(real(dp)), pointer, nopass :: dB02
    procedure(real(dp)), pointer, nopass :: ddB02
    procedure(real(dp)), pointer, nopass :: B03
    procedure(real(dp)), pointer, nopass :: dB03
    procedure(real(dp)), pointer, nopass :: ddB03
  contains
    procedure, public :: get_B0
    procedure, public :: get_dB0
    procedure, public :: delete
  end type bg_magnetic_t

  public :: new_bg_magnetic

contains

  function new_bg_magnetic(default_func) result(bg_magnetic)
    procedure(real(dp)) :: default_func
    type(bg_magnetic_t) :: bg_magnetic
    bg_magnetic%B01 => default_func
    bg_magnetic%B02 => default_func
    bg_magnetic%dB02 => default_func
    bg_magnetic%ddB02 => default_func
    bg_magnetic%B03 => default_func
    bg_magnetic%dB03 => default_func
    bg_magnetic%ddB03 => default_func
  end function new_bg_magnetic


  impure elemental real(dp) function get_B0(this, x)
    class(bg_magnetic_t), intent(in) :: this
    real(dp), intent(in) :: x
    get_B0 = sqrt(this%B01(x)**2 + this%B02(x)**2 + this%B03(x)**2)
  end function get_B0


  impure elemental real(dp) function get_dB0(this, x)
    class(bg_magnetic_t), intent(in) :: this
    real(dp), intent(in) :: x
    get_dB0 = (this%B02(x) * this%dB02(x) + this%B03(x) * this%dB03(x)) / this%get_B0(x)
  end function get_dB0


  pure subroutine delete(this)
    class(bg_magnetic_t), intent(inout) :: this
    nullify(this%B01)
    nullify(this%B02)
    nullify(this%dB02)
    nullify(this%ddB02)
    nullify(this%B03)
    nullify(this%dB03)
    nullify(this%ddB03)
  end subroutine delete

end module mod_bg_magnetic
