module mod_global_variables
  implicit none
  public

  !> Grid-related parameters
  character(:), allocatable         :: geometry
  double precision                  :: x_start, x_end
  integer                           :: gridpts, matrix_gridpts
  integer                           :: integral_gridpts
  double precision                  :: ev_1, ev_2, sigma_1, sigma_2

  real, parameter                   :: dpi = 3.141592653589793238462643383279

  logical, save                     :: mesh_accumulation

  !! Equilibrium-related parameters
  !> Number of Gaussian points
  integer, parameter                   :: n_gauss = 4
  !> Gaussian nodes in the interval [-1, 1]
  double precision, dimension(n_gauss) :: gaussian_nodes = &
                                              (/ -0.861136311594053, &
                                                 -0.339981043584856, &
                                                  0.339981043584856, &
                                                  0.861136311594053  /)
  !> Gaussian weights in the interval [-1, 1]
  double precision, dimension(n_gauss) :: gaussian_weights = &
                                              (/ 0.347854845137454, &
                                                 0.652145154862546, &
                                                 0.652145154862546, &
                                                 0.347854845137454  /)

contains

  subroutine init_variables()
    geometry = "Cartesian"
    x_start  = 0.0d0
    x_end    = 1.0d0
    gridpts  = 11
    matrix_gridpts = 16 * gridpts
    integral_gridpts = gridpts - 1

    mesh_accumulation = .false.

    !> expected values Gaussian
    ev_1 = 7.25d0
    ev_2 = 0.75d0
    !> standard deviations Gaussian
    sigma_1 = 1.0d0
    sigma_2 = 2.0d0


  end subroutine init_variables

  subroutine variables_clean()
    deallocate(geometry)
  end subroutine variables_clean

end module mod_global_variables
