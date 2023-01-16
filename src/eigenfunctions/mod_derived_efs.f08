module mod_derived_efs
  use mod_global_variables, only: dp, str_len_arr
  use mod_settings, only: settings_t
  use mod_base_efs, only: base_ef_t
  implicit none

  private

  type, extends(base_ef_t), public :: derived_ef_t
  contains

    procedure, public :: initialise

  end type derived_ef_t


contains

  pure subroutine initialise(this, name, ef_grid_size, nb_efs)
    class(derived_ef_t), intent(inout) :: this
    character(str_len_arr), intent(in) :: name
    integer, intent(in) :: ef_grid_size
    integer, intent(in) :: nb_efs

    call this%base_ef_t%initialise(name, ef_grid_size, nb_efs)
  end subroutine initialise


end module mod_derived_efs
