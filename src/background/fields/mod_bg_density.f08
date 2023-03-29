module mod_bg_density
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: bg_density_t
    procedure(real(dp)), pointer, nopass :: rho0
    procedure(real(dp)), pointer, nopass :: drho0
  contains
    procedure :: delete
  end type bg_density_t

  public :: new_bg_density

contains

  function new_bg_density(default_func) result(bg_density)
    procedure(real(dp)) :: default_func
    type(bg_density_t) :: bg_density
    bg_density%rho0 => default_func
    bg_density%drho0 => default_func
  end function new_bg_density


  pure subroutine delete(this)
    class(bg_density_t), intent(inout) :: this
    nullify(this%rho0)
    nullify(this%drho0)
  end subroutine delete

end module mod_bg_density
