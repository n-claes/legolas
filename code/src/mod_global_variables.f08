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

  !! Physics parameters
  !> Complex number
  complex(dp), parameter    :: ic = (0.0d0, 1.0d0)
  !> Ratio of specific heats gamma
  real(dp), parameter       :: gamma = 5.0d0/3.0d0
  !> Variable for (gamma - 1)
  real(dp), parameter       :: gamma_1 = 1.0d0 - gamma
  !> Boolean to enable flow
  logical, save             :: flow
  !> Boolean to enable radiative cooling
  logical, save             :: radiative_cooling
  !> Number of points to interpolate the radiative cooling curve
  integer, parameter        :: ncool = 4000
  !> Name of the cooling curve to use
  character(len=str_len)    :: cooling_curve
  !> Boolean to enable external gravity
  logical, save             :: external_gravity
  !> Strength of the gravitational acceleration, can be 'earth' or 'solar'
  character(len=str_len)    :: gravity_type
  !> Boolean to enable thermal conduction
  logical, save             :: thermal_conduction
  !> Boolean to enable resistivity
  logical, save             :: resistivity
  !> Boolean for cgs units
  logical, save             :: cgs_units = .true.

  !! Grid-related parameters
  !> Geometry of the problem: 'Cartesian' or 'cylindrical'
  character(:), allocatable :: geometry
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
  !> Type of equilibrium to set up
  character(:), allocatable :: equilibrium_type

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

contains

  !> Initialises all global variables.
  subroutine initialise_variables()
    use mod_physical_constants

    ! Normalisations
    call set_unit_length(1.0d9)         ! cm
    call set_unit_numberdensity(1.0d9)  ! cm*3
    call set_unit_temperature(1.0d6)    ! K
    call set_normalisations(cgs_units)

    ! Grid related parameters
    geometry = "cylindrical"
    x_start  = 0.0d0
    x_end    = 1.0d0

    call set_gridpts(11)

    ! Mesh accumulation parameters
    mesh_accumulation  = .false.
    !> expected values Gaussian
    ev_1 = 7.25d0
    ev_2 = 0.75d0
    !> standard deviations Gaussian
    sigma_1 = 1.0d0
    sigma_2 = 2.0d0

    ! Physics related parameters
    flow               = .true.
    external_gravity   = .true.
    gravity_type       = "solar"
    radiative_cooling  = .true.
    cooling_curve      = "JCcorona"
    thermal_conduction = .true.
    resistivity        = .true.

    ! Equilibrium related parameters
    equilibrium_type = "adiabatic homogeneous"

  end subroutine initialise_variables

  !> Subroutine to set the gridpoints and dependent gridpoints for
  !! more control at initialisation.
  !! @param[in] gridpts_in  Gridpoints for the initial array
  subroutine set_gridpts(gridpts_in)
    integer, intent(in) :: gridpts_in

    gridpts        = gridpts_in
    gauss_gridpts  = 4*(gridpts - 1)
    matrix_gridpts = 16 * gridpts
  end subroutine set_gridpts

  !> Subroutine to override the matrix gridpoints. This is solely used for
  !! testing purposes, as the matrix gridpoints are controlled by 'gridpts'.
  !! @param[in] gridpts_in  Dimension of the matrix
  subroutine set_matrix_gridpts(gridpts_in)
    integer, intent(in) :: gridpts_in

    matrix_gridpts = gridpts_in
  end subroutine set_matrix_gridpts

  !> Deallocates the variables in this module.
  subroutine variables_clean()
    deallocate(geometry)
    deallocate(equilibrium_type)
  end subroutine variables_clean

end module mod_global_variables
