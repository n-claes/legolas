module mod_global_variables
  use, intrinsic :: iso_fortran_env
  implicit none
  public

  !! Fortran-2008 type standards
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: qp = real128

  !> Default length for strings
  integer, parameter :: str_len = 125

  !! Physical parameters
  real(dp), parameter       :: gamma = 5.0d0/3.0d0
  real(dp), parameter       :: gamma_1 = 1.0d0 - gamma
  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)                  :: k2
  !> Wavenumber in z-direction (Cartesian) or z-direction (cylindrical)
  real(dp)                  :: k3
  ! Radiative cooling
  logical, save             :: radiative_cooling
  integer, parameter        :: ncool = 4000
  character(len=str_len)    :: cooling_curve
  ! External gravity
  logical, save             :: external_gravity
  ! Thermal conduction
  logical, save             :: thermal_conduction

  !> Grid-related parameters
  character(:), allocatable         :: geometry
  real(dp)                          :: x_start, x_end
  integer                           :: gridpts, matrix_gridpts
  integer                           :: integral_gridpts
  real(dp)                          :: ev_1, ev_2, sigma_1, sigma_2

  real(dp), parameter               :: dpi = 3.141592653589793238462643383279

  logical, save                     :: mesh_accumulation

  !! Equilibrium-related parameters
  !> Number of Gaussian points
  integer, parameter                :: n_gauss = 4
  !> Gaussian nodes in the interval [-1, 1]
  real(dp), dimension(n_gauss) :: gaussian_nodes = &
                                      (/ -0.861136311594053, &
                                         -0.339981043584856, &
                                          0.339981043584856, &
                                          0.861136311594053  /)
  !> Gaussian weights in the interval [-1, 1]
  real(dp), dimension(n_gauss) :: gaussian_weights = &
                                      (/ 0.347854845137454, &
                                         0.652145154862546, &
                                         0.652145154862546, &
                                         0.347854845137454  /)

  !! Block-related parameters
  !> Number of equations
  integer, parameter        :: nb_eqs = 8
  !> Dimension of one finite element integral block, e.g. A(1,1).
  !  See page 190 of advanced MHD
  integer, parameter        :: dim_integralblock = 2
  !> Dimension of one quadblock subblock, 16x16
  !  One subblock has elements A(1,1) -> A(8,8), where every element A(i,i)
  !  is in itself a block of dimension 2x2 (dim_integralblock above)
  integer, parameter        :: dim_subblock = nb_eqs * dim_integralblock
  !> Dimension of the quadblock, containing 4 subblocks
  !  The quadblock is the matrix that is shifted on the main diagonal
  integer, parameter        :: dim_quadblock = 2*dim_subblock

contains

  subroutine initialise_variables()
    geometry = "Cartesian"
    x_start  = 0.0d0
    x_end    = 1.0d0
    gridpts  = 11
    matrix_gridpts = 16 * gridpts
    integral_gridpts = gridpts - 1

    mesh_accumulation  = .false.
    external_gravity   = .true.
    radiative_cooling  = .false.
    cooling_curve = "SPEX_DM"
    thermal_conduction = .true.

    !> expected values Gaussian
    ev_1 = 7.25d0
    ev_2 = 0.75d0
    !> standard deviations Gaussian
    sigma_1 = 1.0d0
    sigma_2 = 2.0d0


  end subroutine initialise_variables

  subroutine variables_clean()
    deallocate(geometry)
  end subroutine variables_clean

end module mod_global_variables
