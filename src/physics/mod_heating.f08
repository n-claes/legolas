module mod_heating
  use mod_global_variables, only: dp
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_function_utils, only: zero_func
  implicit none

  private

  type(settings_t), pointer :: settings => null()
  type(background_t), pointer :: background => null()

  type, public :: heating_t
    procedure(real(dp)), pointer, nopass :: H
    procedure(real(dp)), pointer, nopass :: dHdT
    procedure(real(dp)), pointer, nopass :: dHdrho

  contains
    procedure, public :: delete
  end type heating_t

  public :: new_heating

contains

  function new_heating(settings_tgt, background_tgt) result(heating)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(heating_t) :: heating
    settings => settings_tgt
    background => background_tgt

    heating%H => zero_func
    heating%dHdT => zero_func
    heating%dHdrho => zero_func
  end function new_heating

  subroutine delete(this)
    class(heating_t), intent(inout) :: this
    nullify(settings)
    nullify(background)
    nullify(this%H)
    nullify(this%dHdT)
    nullify(this%dHdrho)
  end subroutine delete

end module mod_heating
