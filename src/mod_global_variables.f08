! =============================================================================
!> All global variables are defined in this module, and can be accessed
!! throughout the entire code. Most variables are first initialised
!! to a default value, and are properly set either in the parfile
!! or in an equilibrium submodule.
module mod_global_variables
  use, intrinsic :: iso_fortran_env
  implicit none

  public

  !> single-precision value
  integer, parameter :: sp = real32
  !> double-precision value
  integer, parameter :: dp = real64
  !> quadruple-precision value
  integer, parameter :: qp = real128
  !> default length for strings
  integer, parameter :: str_len = 125
  !> default length for strings in arrays
  integer, parameter :: str_len_arr = 16

  !> values smaller than this are forced to zero
  real(dp), parameter :: dp_LIMIT = 1.0d-12
  !> NaN value (ieee_quiet_nan)
  real(dp), protected :: NaN

  !> complex number i
  complex(dp), parameter    :: ic = (0.0d0, 1.0d0)
  !> complex real
  complex(dp), parameter    :: ir = (1.0d0, 0.0d0)
  !> boolean for cgs units, defaults to <tt>True</tt>
  logical, save             :: cgs_units
  !> ratio of specific heats gamma, defaults to 5/3
  real(dp), protected       :: gamma
  !> variable for (gamma - 1)
  real(dp), protected       :: gamma_1
  !> boolean for incompressible approximation, defaults to <tt>False</tt>
  logical, save             :: incompressible
  !> boolean for flow, defaults to <tt>False</tt>
  logical, save             :: flow
  !> boolean for radiative cooling, defaults to <tt>False</tt>
  logical, save             :: radiative_cooling
  !> number of points to interpolate the radiative cooling curve, defaults to 4000
  integer                   :: ncool
  !> name of the cooling curve to use, defaults to \p jc_corona
  character(len=str_len)    :: cooling_curve
  !> boolean for external gravity, defaults to <tt>False</tt>
  logical, save             :: external_gravity
  !> boolean for thermal conduction, defaults to <tt>False</tt>
  logical, save             :: thermal_conduction
  !> boolean to set a fixed value for parallel conduction, defaults to <tt>False</tt>
  logical, save             :: use_fixed_tc_para
  !> defines the fixed value for parallel conduction, defaults to 0
  real(dp)                  :: fixed_tc_para_value
  !> boolean to set a fixed value for perpendicular conduction, defaults to <tt>False</tt>
  logical, save             :: use_fixed_tc_perp
  !> defines the fixed value for perpendicular conduction, defaults to 0
  real(dp)                  :: fixed_tc_perp_value
  !> boolean for resistivity, defaults to <tt>False</tt>
  logical, save             :: resistivity
  !> boolean to set a fixed value for the resistivity, defaults to <tt>False</tt>
  logical, save             :: use_fixed_resistivity
  !> defines the fixed value for the resistivity, defaults to 0
  real(dp)                  :: fixed_eta_value
  !> boolean for a resistivity profile that goes to zero on edges, defaults to <tt>False</tt>
  logical, save             :: use_eta_dropoff
  !> distance between the grid edge and start of dropoff, defaults to <tt>0.05</tt>
  real(dp)                  :: dropoff_edge_dist
  !> width of the dropoff region, defaults to <tt>0.1</tt>
  real(dp)                  :: dropoff_width
  !> boolean for viscosity, defaults to <tt>False</tt>
  logical, save             :: viscosity
  !> boolean to include viscoous heating, defaults to <tt>False</tt>
  logical, save             :: viscous_heating
  !> defines the fixed value for either the dynamic or kinematic viscosity, defaults to 0
  real(dp)                  :: viscosity_value
  !> boolean to use Hall MHD, defaults to <tt>False</tt>
  logical, save             :: hall_mhd
  !> boolean to use dropoff profile for Hall parameter, defaults to <tt>False </tt>
  logical, save             :: hall_dropoff
  !> boolean to use electron pressure in Ohm's law, defaults to <tt>False</tt>
  logical, save             :: elec_pressure
  !> boolean to use electron inertia in Ohm's law, defaults to <tt>False</tt>
  logical, save             :: elec_inertia
  !> boolean to use dropoff profile for inertia parameter, defaults to <tt>False </tt>
  logical, save             :: inertia_dropoff
  !> defines the geometry of the problem, defaults depend on chosen equilibrium
  character(len=str_len)    :: geometry
  !> defines the presence of a coaxial inner boundary for a cylindrical geometry,
  !! defaults to <tt>False</tt>
  logical, save             :: coaxial
  !> start value of the base grid, defaults depend on chosen equilibrium
  real(dp)                  :: x_start
  !> end value of the base grid, defaults depend on chosen equilibrium
  real(dp)                  :: x_end
  !> number of gridpoints in the base grid, defaults to 31
  integer, protected        :: gridpts
  !> boolean to force r=0 in cylindrical geometry, defaults to <tt>False</tt>
  logical, save             :: force_r0
  !> number of gridpoints in the gaussian grid, automatically set by \p gridpts
  integer, protected        :: gauss_gridpts
  !> size of the A and B matrices, automatically set by \p gridpts
  integer, protected        :: matrix_gridpts
  !> size of a single eigenfunction array, automatically set by \p gridpts
  integer, protected        :: ef_gridpts

  !> number of Gaussian nodes
  integer, parameter           :: n_gauss = 4
  !> values for the Gaussian nodes in [-1, 1]
  real(dp), dimension(n_gauss) :: gaussian_nodes = &
                                      (/ -0.861136311594053, &
                                         -0.339981043584856, &
                                          0.339981043584856, &
                                          0.861136311594053  /)
  !> weights for the Gaussian nodes in [-1, 1]
  real(dp), dimension(n_gauss) :: gaussian_weights = &
                                      (/ 0.347854845137454, &
                                         0.652145154862546, &
                                         0.652145154862546, &
                                         0.347854845137454  /)

  !> name of the equilibrium to set up, determines the submodule, defaults to \p "adiabatic_homo"
  character(len=str_len)       :: equilibrium_type
  !> type of boundary conditions, defaults to \p "wall"
  character(len=str_len)       :: boundary_type
  !> use default values for parameters in the chosen submodule, defaults to <tt>True</tt>
  logical, save                :: use_defaults
  !> boolean for spurious eigenvalue removal, defaults to <tt>False</tt>
  logical, save                :: remove_spurious_eigenvalues
  !> amount of eigenvalues to remove on each side of the imaginary axis, defaults to 1
  integer                      :: nb_spurious_eigenvalues

  !> total number of equations
  integer, parameter        :: nb_eqs = 8
  !> dimension of one finite element integral block, e.g. A(1, 2)
  integer, parameter        :: dim_integralblock = 2
  !> dimension of one subblock, 4 of these define a quadblock
  integer, parameter        :: dim_subblock = nb_eqs * dim_integralblock
  !> dimension of one quadblock, this is the block shifted along the main diagonal
  integer, parameter        :: dim_quadblock = 2*dim_subblock

  !> boolean to write both matrices to the datfile, defaults to <tt>False</tt>
  logical, save             :: write_matrices
  !> boolean to write the eigenfunctions to the datfile, defaults to <tt>True</tt>
  logical, save             :: write_eigenfunctions
  !> boolean to call the Python wrapper and plot the results, defaults to <tt>True</tt>
  logical, save             :: show_results
  !> base name for the datfile, defaults to \p "datfile"
  character(len=str_len)    :: basename_datfile
  !> base name for the logfile, defaults to ""
  character(len=str_len)    :: basename_logfile
  !> path to the output folder, defaults to "."
  character(len=str_len)    :: output_folder
  !> sets the logging level, defaults to 2 (errors, warnings and info)
  integer                   :: logging_level
  !> boolean for doing a dry run, defaults to <tt>False</tt>
  logical, save             :: dry_run

  !> method to solve the eigenvalue problem, defaults to <tt>"QR-invert"</tt>
  character(len=str_len)    :: solver
  !> selected mode for ARPACK, defaults to <tt>"standard"</tt>
  character(len=str_len)    :: arpack_mode
  !> number of eigenvalues to calculate with ARPACK, defaults to 100
  integer                   :: number_of_eigenvalues
  !> which eigenvalues to calculate, defaults to <tt>"LM"<tt> (largest magnitude)
  character(len=2)          :: which_eigenvalues
  !> maximum number of Arnoldi iterations
  integer                   :: maxiter
  !> sigma value, only used in shift-invert mode
  complex(dp)               :: sigma

contains


  !> Initialises the global variables in this module.
  !! All variables in this module are first set to their default values.
  !! These are either regular values or NaN, the latter in case variables
  !! must be explicitly set in the parfile or equilibrium submodule.
  !! @note  The variables <tt>geometry</tt>, <tt>x_start</tt> and
  !!        <tt>x_end</tt> are explicitly set to NaN.
  subroutine initialise_globals()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    NaN = ieee_value(NaN, ieee_quiet_nan)

    !! physics variables
    cgs_units = .true.
    gamma = 5.0d0/3.0d0
    incompressible = .false.
    flow = .false.
    radiative_cooling = .false.
    ncool = 4000
    cooling_curve = 'jc_corona'
    external_gravity = .false.
    thermal_conduction = .false.
    use_fixed_tc_para = .false.
    fixed_tc_para_value = 0.0d0
    use_fixed_tc_perp = .false.
    fixed_tc_perp_value = 0.0d0
    resistivity = .false.
    use_fixed_resistivity = .false.
    fixed_eta_value = 0.0d0
    use_eta_dropoff = .false.
    dropoff_edge_dist = 0.05d0
    dropoff_width = 0.1d0
    viscosity = .false.
    viscous_heating = .false.
    viscosity_value = 0.0d0
    hall_mhd = .false.
    hall_dropoff = .false.
    elec_pressure = .false.
    elec_inertia = .false.
    inertia_dropoff = .false.

    !! grid variables
    ! do not initialise these three so they MUST be set in submodules/parfile
    geometry = ""
    coaxial = .false.
    x_start = NaN
    x_end = NaN
    gridpts = 31
    force_r0 = .false.
    call set_gridpts(gridpts)

    !! equilibrium variables
    equilibrium_type = 'adiabatic_homo'
    boundary_type = 'wall'
    use_defaults = .true.
    remove_spurious_eigenvalues = .false.
    nb_spurious_eigenvalues = 1

    !! post-processing parameters
    write_matrices = .false.
    write_eigenfunctions = .true.
    show_results = .true.
    logging_level = 2
    dry_run = .false.

    !! file-saving variables
    basename_datfile = "datfile"
    basename_logfile = ""
    output_folder = "output"

    !! solution method variables
    solver = "QR-invert"
    arpack_mode = "standard"
    number_of_eigenvalues = 100
    which_eigenvalues = "LM"
    ! this defaults to 10*N, but we technically don't know N yet
    maxiter = 0
  end subroutine initialise_globals


  !> Sets the ratio of specific heats gamma and its corresponding
  !! value gamma - 1.
  subroutine set_gamma(gamma_in)
    !> the value for gamma
    real(dp), intent(in)    :: gamma_in

    if (incompressible) then
      gamma = 1.0d8
    else
      gamma = gamma_in
    end if
    gamma_1 = gamma - 1.0d0
  end subroutine set_gamma


  !> Sets all gridpoint-related variables: sets the base number of gridpoints,
  !! the gridpoints of the Gaussian grid, matrix sizes and size of the
  !! eigenfunction arrays.
  subroutine set_gridpts(gridpts_in)
    !> amount of gridpoints for the base grid
    integer, intent(in) :: gridpts_in

    gridpts        = gridpts_in
    gauss_gridpts  = 4*(gridpts - 1)
    matrix_gridpts = 16 * gridpts
    ef_gridpts     = 2*gridpts - 1
  end subroutine set_gridpts

end module mod_global_variables
