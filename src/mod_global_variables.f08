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
  integer, parameter :: str_len = 500
  !> default length for strings in arrays
  integer, parameter :: str_len_arr = 16

  !> values smaller than this are forced to zero, bit higher than machine precision
  real(dp), parameter :: dp_LIMIT = 5.0d-15
  !> NaN value (ieee_quiet_nan)
  real(dp), protected :: NaN

  !> complex number i
  complex(dp), parameter    :: ic = (0.0d0, 1.0d0)
  !> complex real
  complex(dp), parameter    :: ir = (1.0d0, 0.0d0)
  !> boolean for cgs units, defaults to <tt>True</tt>
  logical, save             :: cgs_units
  !> distance between the grid edge and start of dropoff, defaults to <tt>0.05</tt>
  real(dp)                  :: dropoff_edge_dist
  !> width of the dropoff region, defaults to <tt>0.1</tt>
  real(dp)                  :: dropoff_width
  !> boolean to use Hall MHD, defaults to <tt>False</tt>
  logical, save             :: hall_mhd
  !> boolean to use substitution for Hall elements (as presented in the paper), defaults to <tt>True</tt>
  logical, save             :: hall_substitution
  !> boolean to use dropoff profile for Hall parameter, defaults to <tt>False</tt>
  logical, save             :: hall_dropoff
  !> boolean to use electron inertia in Ohm's law, defaults to <tt>False</tt>
  logical, save             :: elec_inertia
  !> boolean to use dropoff profile for inertia parameter, defaults to <tt>False</tt>
  logical, save             :: inertia_dropoff
  !> fraction of number of electrons to number of all particles (between 0 and 1), defaults to <tt>0.5</tt>
  real(dp)                  :: electron_fraction

  !> number of Gaussian nodes
  integer, parameter           :: n_gauss = 4
  !> values for the Gaussian nodes in [-1, 1]
  real(dp), parameter :: gaussian_nodes(n_gauss) = [ &
    -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053 &
  ]
  !> weights for the Gaussian nodes in [-1, 1]
  real(dp), parameter :: gaussian_weights(n_gauss) = [ &
    0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454 &
  ]

  !> sets the logging level, defaults to 2 (errors, warnings and info)
  integer                   :: logging_level

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
    dropoff_edge_dist = 0.05d0
    dropoff_width = 0.1d0
    hall_mhd = .false.
    hall_substitution = .true.
    hall_dropoff = .false.
    elec_inertia = .false.
    inertia_dropoff = .false.
    electron_fraction = 0.5d0

    !! post-processing parameters
    logging_level = 2
  end subroutine initialise_globals

end module mod_global_variables
