!
! MODULE: mod_global_variables
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module containing all global variables.
!
module mod_global_variables
  use, intrinsic :: iso_fortran_env
  implicit none

  public

  !! Fortran-2008 type standards
  !> Single-precision value
  integer, parameter :: sp = real32
  !> Double-precision value
  integer, parameter :: dp = real64
  !> Quadruple-precision value
  integer, parameter :: qp = real128
  !> Default length for strings
  integer, parameter :: str_len = 125

  !> Values smaller than this are forced to zero
  real(dp), parameter :: dp_LIMIT = 1.0d-12

  !! Physics parameters
  !> Complex number
  complex(dp), parameter    :: ic = (0.0d0, 1.0d0)
  !> Complex real
  complex(dp), parameter    :: ir = (1.0d0, 0.0d0)
  !> Boolean for cgs units
  logical, save             :: cgs_units
  !> Ratio of specific heats gamma
  real(dp), protected       :: gamma
  !> Variable for (gamma - 1)
  real(dp), protected       :: gamma_1
  !> Boolean to enable flow
  logical, save             :: flow
  !> Boolean to enable radiative cooling
  logical, save             :: radiative_cooling
  !> Number of points to interpolate the radiative cooling curve
  integer                   :: ncool
  !> Name of the cooling curve to use
  character(len=str_len)    :: cooling_curve
  !> Boolean for external gravity
  logical, save             :: external_gravity
  !> Boolean to enable thermal conduction
  logical, save             :: thermal_conduction
  !> Boolean to enable resistivity
  logical, save             :: resistivity
  !> Boolean to set a fixed resistivity instead of a temperature dependence
  logical, save             :: use_fixed_resistivity
  !> Sets the fixed value for the resistivity
  real(dp)                  :: fixed_eta_value
  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)                  :: k2
  !> Wavenumber in z-direction (Cartesian) or z-direction (cylindrical)
  real(dp)                  :: k3

  !! Grid-related parameters
  !> Geometry of the problem: 'Cartesian' or 'cylindrical'
  character(len=str_len)    :: geometry
  !> Starting value for the x-array
  real(dp)                  :: x_start
  !> Final value of the x-array
  real(dp)                  :: x_end
  !> Amount of gridpoints of the initial array
  integer, protected        :: gridpts
  !> Amount of gridpoint of the gaussian array
  integer, protected        :: gauss_gridpts
  !> Gridpoints of matrices A and B, equal to 16 * gridpts
  integer, protected        :: matrix_gridpts
  !> Amount of gridpoints of an eigenfunction array
  integer, protected        :: ef_gridpts

  !! Mesh-accumulation parameters
  !> Boolean to enable mesh accumulation
  logical, save             :: mesh_accumulation
  !> Expected value 1 for the Gaussian function
  real(dp)                  :: ev_1
  !> Expected value 2 for the Gaussian function
  real(dp)                  :: ev_2
  !> Standard deviation 1 for the Gaussian function
  real(dp)                  :: sigma_1
  !> Standard deviation 2 for the Gaussian function
  real(dp)                  :: sigma_2

  !! Equilibrium-related parameters
  !> Number of Gaussian points
  integer, parameter           :: n_gauss = 4
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

  !> Type of equilibrium to set up
  character(len=str_len)       :: equilibrium_type
  !> Type of boundary to use
  character(len=str_len)       :: boundary_type

  !! Block-related parameters
  !> Number of equations
  integer, parameter        :: nb_eqs = 8
  !> Dimension of one finite element integral block, e.g. A(1,1).
  !! See page 190 of advanced MHD.
  integer, parameter        :: dim_integralblock = 2
  !> Dimension of one quadblock subblock, 16x16
  !! One subblock has elements A(1,1) -> A(8,8), where every element A(i,i)
  !! is in itself a block of dimension 2x2 (dim_integralblock above) due to
  !! quadratic and cubic finite elements.
  integer, parameter        :: dim_subblock = nb_eqs * dim_integralblock
  !> Dimension of the quadblock, containing 4 subblocks
  !! The quadblock is the matrix that is shifted on the main diagonal
  integer, parameter        :: dim_quadblock = 2*dim_subblock

  !! Data IO-related parameters (savelist)
  !> Write matrices A and B to file when finished
  logical, save             :: write_matrices
  !> Write eigenvectors to file when finished
  logical, save             :: write_eigenvectors
  !> Write eigenfunctions to file when finished
  logical, save             :: write_eigenfunctions
  !> Call python script when finishing to plot results
  logical, save             :: show_results
  !> Plot matrices A and B when finished
  logical, save             :: show_matrices
  !> Plot eigenfunctions when finished
  logical, save             :: show_eigenfunctions

  !! Related parameters for filenames (filelist)
  !> Name for the eigenvalues file
  character(len=str_len)    :: savename_eigenvalues
  !> Name for the eigenfunction grid file
  character(len=str_len)    :: savename_efgrid
  !> Name for the matrix files
  character(len=str_len)    :: savename_matrix
  !> Name for the eigenvector files
  character(len=str_len)    :: savename_eigenvectors
  !> Name for the eigenfunction files
  character(len=str_len)    :: savename_eigenfunctions
  !> Name for the configuration file
  character(len=str_len)    :: savename_config

contains

  !> Subroutine to set gamma and (gamma - 1)
  !! @param[in] gamma_in  The ratio of specific heats gamma
  subroutine set_gamma(gamma_in)
    real(dp), intent(in)    :: gamma_in

    gamma   = gamma_in
    gamma_1 = gamma - 1.0d0
  end subroutine set_gamma

  !> Subroutine to set the gridpoints and dependent gridpoints for
  !! more control at initialisation.
  !! @param[in] gridpts_in  Gridpoints for the initial array
  subroutine set_gridpts(gridpts_in)
    integer, intent(in) :: gridpts_in

    gridpts        = gridpts_in
    gauss_gridpts  = 4*(gridpts - 1)
    matrix_gridpts = 16 * gridpts
    ef_gridpts     = 2*gridpts - 1
  end subroutine set_gridpts

  !> Subroutine to override the matrix gridpoints. This is solely used for
  !! testing purposes, as the matrix gridpoints are controlled by 'gridpts'.
  !! @param[in] gridpts_in  Dimension of the matrix
  subroutine set_matrix_gridpts(gridpts_in)
    integer, intent(in) :: gridpts_in

    matrix_gridpts = gridpts_in
  end subroutine set_matrix_gridpts

end module mod_global_variables
