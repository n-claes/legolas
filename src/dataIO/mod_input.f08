module mod_input
  use mod_global_variables
  use mod_equilibrium_params
  implicit none

  private

  integer       :: unit_par = 101
  logical, save :: parfile_present

  public :: read_parfile
  public :: get_parfile

contains

  !> Reads in the supplied parfile and sets all global variables accordingly.
  !! @param[in] parfile   The name of the parfile
  subroutine read_parfile(parfile)
    use mod_physical_constants, only: set_unit_length, set_unit_numberdensity, &
                                      set_unit_temperature, set_normalisations

    character(len=*), intent(in)  :: parfile

    real(dp)    :: mhd_gamma
    ! integer     :: big      !< incommpressible limit test
    real(dp)    :: unit_length, unit_numberdensity, unit_temperature, &
                   unit_velocity
    integer     :: gridpoints

    namelist /gridlist/ geometry, x_start, x_end, gridpoints
    namelist /meshlist/ mesh_accumulation, ev_1, ev_2, sigma_1, sigma_2
    namelist /physicslist/ mhd_gamma, flow, radiative_cooling, ncool, &
                           cooling_curve, external_gravity, &
                           thermal_conduction, resistivity, &
                           use_fixed_resistivity, fixed_eta_value
    namelist /unitslist/ cgs_units, unit_length, unit_numberdensity, &
                         unit_temperature, unit_velocity
    namelist /equilibriumlist/ equilibrium_type, boundary_type, use_defaults
    namelist /savelist/ write_matrices, write_eigenvectors, &
                        write_eigenfunctions, show_results, show_matrices, &
                        show_eigenfunctions
    namelist /filelist/ savename_eigenvalues, savename_efgrid, &
                        savename_matrix, savename_eigenvectors, &
                        savename_eigenfunctions, savename_config
    namelist /paramlist/ k2, k3, cte_rho0, cte_T0, cte_B02, cte_B03, cte_v02, cte_v03, &
                         cte_p0, p1, p2, p3, p4, p5, p6, p7, p8, &
                         alpha, beta, delta, theta, tau, lambda, nu, &
                         r0, rc, rj, Bth0, Bz0, V, j0, g

    parfile_present = .true.
    if (parfile == "") then
      parfile_present = .false.
    end if

    !! Set defaults
    !> Gridlist defaults
    geometry = "Cartesian"          !< geometry of the problem
    x_start  = 0.0d0                !< start of the grid
    x_end    = 1.0d0                !< end of the grid
    gridpoints = 31                 !< gridpoints for regular grid

    !> Meshlist defaults
    mesh_accumulation = .false.
    ev_1 = 1.25d0                   !< expected value gaussian 1
    ev_2 = 1.25d0                   !< expected value gaussian 2
    sigma_1 = 1.0d0                 !< standard deviation gaussian 1
    sigma_2 = 2.0d0                 !< standard deviation gaussian 2

    !> Physicslist defaults
    mhd_gamma = 5.0d0 / 3.0d0       !< ratio of specific heats
    ! mhd_gamma = real(huge(big))     !< incompressible limit test
    flow  = .false.                 !< use flow
    radiative_cooling = .false.     !< use radiative cooling
    ncool = 4000                    !< points for cooling curve interpolation
    cooling_curve = "JCcorona"      !< radiative cooling curve to use
    external_gravity = .false.      !< use external gravity
    thermal_conduction = .false.    !< use thermal conduction
    resistivity = .false.           !< use resistivity
    use_fixed_resistivity = .false. !< use fixed resistivity
    fixed_eta_value = 0.0d0         !< value for fixed resistivity

    !> Unitslist defaults
    cgs_units = .true.              !< use cgs units
    unit_length = 1.0d9             !< length unit in cm (cgs) or m (mks)
    unit_numberdensity = 1.0d9      !< cm**-3 (cgs) or m**-3 (mks)
    unit_temperature = 1.0d6        !< temperature unit in kelvin
    unit_velocity = 0.0d0           !< velocity unit in cm/s (cgs) or m/s (mks)

    !> Equilibriumlist defaults
    equilibrium_type = "Adiabatic homogeneous"  !< precoded equilibrium to use
    boundary_type = 'wall'                      !< type of boundary condition
    use_defaults = .true.                       !< use defaults for equilibrium parameters

    !> Savelist defaults
    write_matrices = .false.        !< write matrices A and B when finished
    write_eigenvectors = .false.    !< writes eigenvectors to file
    write_eigenfunctions = .false.  !< writes eigenfunctions to file
    show_results = .true.           !< plot spectrum when finished
    show_matrices = .false.         !< plot matrices A and B when finished
    show_eigenfunctions = .false.   !< plots the eigenfunctions when finished

    !> Filelist defaults
    savename_eigenvalues = "eigenvalues"
    savename_efgrid = "ef_grid"
    savename_matrix = "matrix"
    savename_eigenvectors = "eigenvectors"
    savename_eigenfunctions = "eigenfunctions"
    savename_config = "configuration"

    !! Initialise equilibrium parameters to nan. These are controlled
    !! using the par file and/or the equilibrium submodules.
    call init_equilibrium_params()

    !! Read parfile, if supplied
    if (parfile_present) then
      open(unit_par, file=trim(parfile), status='old')
      !! Start reading namelists, rewind so they can appear out of order
            rewind(unit_par)
            read(unit_par, gridlist, end=1001)

      1001  rewind(unit_par)
            read(unit_par, meshlist, end=1002)

      1002  rewind(unit_par)
            read(unit_par, physicslist, end=1003)

      1003  rewind(unit_par)
            read(unit_par, unitslist, end=1004)

      1004  rewind(unit_par)
            read(unit_par, equilibriumlist, end=1005)

      1005  rewind(unit_par)
            read(unit_par, savelist, end=1006)

      1006  rewind(unit_par)
            read(unit_par, filelist, end=1007)

      1007  rewind(unit_par)
            read(unit_par, paramlist, end=1008)

      1008  close(unit_par)
    end if

    !> Set up grid and normalisations
    call set_gridpts(gridpoints)
    call set_gamma(mhd_gamma)

    call set_unit_length(unit_length)
    call set_unit_numberdensity(unit_numberdensity)
    call set_unit_temperature(unit_temperature)
    call set_normalisations(cgs_units)

  end subroutine read_parfile


  !> Retrieves the parfile passed as command line argument.
  !! @param[out] filename_par   The name of the parfile
  subroutine get_parfile(filename_par)
    character(len=str_len), intent(out) :: filename_par

    integer                             :: num_args, i
    character(len=20), allocatable      :: args(:)
    logical                             :: file_exists


    num_args = command_argument_count()
    allocate(args(num_args))
    do i = 1, num_args
      call get_command_argument(i, args(i))
    end do

    filename_par = ""

    do i = 1, num_args, 2
      select case(args(i))
      case('-i')
        filename_par = args(i+1)
      case default
        write(*, *) "Unable to read in command line arguments."
        stop
      end select
    end do

    if (filename_par == "") then
      write(*, *) "No parfile supplied, using default configuration."
      return
    end if

    !! Check if supplied file exists
    inquire(file=trim(filename_par), exist=file_exists)
    if (.not. file_exists) then
      write(*, *) "Parfile not found: ", trim(filename_par)
      stop
    end if
  end subroutine get_parfile

end module mod_input
