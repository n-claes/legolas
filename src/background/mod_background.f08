module mod_background
  use mod_bg_profiles, only: zero
  use mod_bg_density, only: bg_density_t, new_bg_density
  implicit none

  private

  type, public :: background_t
    type(bg_density_t) :: density
  end type background_t

  public :: new_background

contains

  function new_background() result(background)
    type(background_t) :: background

    background%density = new_bg_density(default_func=zero)
  end function new_background

end module mod_background
