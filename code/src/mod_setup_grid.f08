module mod_setup_grid
  implicit none

contains

  subroutine initialise_grid(grid))
    use mod_global_variables
    
    real, intent(out)            :: grid(integral_gridpts)

    dx = (x_end - x_start) / integral_gridpts
    do i = 1, integral_gridpts
      grid(i) = x_start + (i - 1)*dx

  end subroutine




end module mod_setup_grid
