module mod_setup_grid
  implicit none

contains

  subroutine initialise_grid(grid)
    use mod_global_variables

    double precision, intent(inout)  :: grid(gridpts)
    integer                          :: i
    double precision                 :: dx


    ! minus one here to include x_end
    dx = (x_end - x_start) / (integral_gridpts-1)
    do i = 1, gridpts
      grid(i) = x_start + (i - 1)*dx
    end do

    call accumulate_mesh(grid)

  end subroutine initialise_grid

  !> Subroutine to re-grid the mesh to a non-uniform spacing.
  !! This is based on two Gaussian curves with known widths (from
  !! mod_global_variables); the grid is accumulated near each Gaussian maximum,
  !! using equidistribution based on the integral under the curve defined
  !! by the function gaussian().
  subroutine accumulate_mesh(grid)
    use mod_global_variables

    double precision, intent(inout)  :: grid(gridpts)
    integer                      :: i, integral_gridpts_x5
    double precision             :: dx, xi, gauss_x, bgf, fact
    double precision             :: x_sum, int_eval_1, int_eval_1_inv

    bgf  = 0.3d0 !background
    fact = 0.1d0

    ! first evaluation of integral
    integral_gridpts_x5 = 5*integral_gridpts
    dx = (grid(gridpts) - grid(1)) / float(integral_gridpts_x5)
    xi = grid(1)
    x_sum = 0.0d0

    do i = 1, integral_gridpts_x5
      gauss_x = gaussian(xi, bgf, fact)
      x_sum   = x_sum + gauss_x * dx
      xi      = xi + dx
    end do

    int_eval_1 = x_sum
    int_eval_1_inv = 1.0d0 / x_sum



    ! second evaluation




    return

  end subroutine accumulate_mesh

  !> Function to calculate a Gaussian curve based on known widths and
  !! expected values (from mod_global_variables). The Gaussian is evaluated
  !! in x.
  function gaussian(x, bgf, fact) result(f_gauss)
    use mod_global_variables

    double precision, intent(in)    :: x, bgf, fact
    double precision                :: f_gauss
    double precision                :: gauss_1, gauss_2, norm_1, norm_2

    norm_1 = 1.0d0 / (sigma_1 * sqrt(2*dpi))
    norm_2 = 1.0d0 / (sigma_2 * sqrt(2*dpi))

    gauss_1 = norm_1 * exp(-0.5d0 * ((x - ev_1) / sigma_1)**2)
    gauss_2 = norm_2 * exp(-0.5d0 * ((x - ev_2) / sigma_2)**2)

    f_gauss = bgf + (1.0d0 - bgf) * (gauss_1 + fact*gauss_2) / fact

  end function gaussian


  subroutine grid_clean()
    return
  end subroutine grid_clean




end module mod_setup_grid
