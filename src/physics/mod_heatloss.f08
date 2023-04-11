module mod_heatloss
  use mod_global_variables, only: dp
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_radiative_cooling, only: cooling_t
  implicit none
  private

  type(settings_t), pointer :: settings => null()
  type(background_t), pointer :: background => null()
  type(cooling_t), pointer :: cooling =>  null()

  type, public :: heatloss_t
    procedure(real(dp)), pointer, nopass :: L0
    procedure(real(dp)), pointer, nopass :: dLdT
    procedure(real(dp)), pointer, nopass :: dLdrho

  contains
    procedure, public :: delete
  end type heatloss_t

  public :: new_heatloss

contains

  function new_heatloss(settings_tgt, background_tgt, cooling_tgt) result(heatloss)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(cooling_t), target, intent(in) :: cooling_tgt
    type(heatloss_t) :: heatloss

    settings => settings_tgt
    background => background_tgt
    cooling => cooling_tgt

    heatloss%L0 => get_L0
    heatloss%dLdT => get_dLdT
    heatloss%dLdrho => get_dLdrho
  end function new_heatloss


  real(dp) function get_L0(x)
    real(dp), intent(in) :: x
    get_L0 = 0.0_dp
  end function get_L0


  real(dp) function get_dLdT(x)
    real(dp), intent(in) :: x
    get_dLdT = background%density%rho0(x) * cooling%dlambdadT(x)
  end function get_dLdT


  real(dp) function get_dLdrho(x)
    real(dp), intent(in) :: x
    get_dLdrho = cooling%lambdaT(x)
  end function get_dLdrho


  subroutine delete(this)
    class(heatloss_t), intent(inout) :: this
    nullify(settings)
    nullify(background)
    nullify(cooling)
    nullify(this%L0)
    nullify(this%dLdT)
    nullify(this%dLdrho)
  end subroutine delete

end module mod_heatloss
