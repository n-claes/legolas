module mod_setup_grid
  implicit none

contains

  subroutine initialise_grid(grid)
    use mod_global_variables

    real, intent(out)            :: grid(integral_gridpts)
    integer                      :: i
    double precision             :: dx

    ! minus one here to include x_end
    dx = (x_end - x_start) / (integral_gridpts-1)
    do i = 1, integral_gridpts
      grid(i) = x_start + (i - 1)*dx
    end do

  end subroutine initialise_grid




end module mod_setup_grid
