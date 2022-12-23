! =============================================================================
!> Module to handle parfile reading.
!! Contains subroutines to retrieve the parfile based on the commandline arguments
!! and to read the parfile, setting the global variables.
module mod_input
  use mod_global_variables, only: dp, str_len, NaN
  use mod_logging, only: log_message
  use mod_settings, only: settings_t
  implicit none

  private

  integer, parameter :: unit_par = 101

  public :: read_parfile
  public :: get_parfile

contains


  subroutine read_parfile(parfile, settings)
    character(len=*), intent(in) :: parfile
    type(settings_t), intent(inout) :: settings
    integer :: iostat


    if (parfile == "") then
      call log_message("no parfile supplied!", level="error")
      return
    end if

    open(unit_par, file=trim(parfile), status="old")
      do while (iostat == 0)
        ! rewind after reading so namelists can appear out of order
        call read_gridlist(unit_par, settings, iostat)
        rewind(unit_par)
        call read_savelist(unit_par, settings, iostat)
        rewind(unit_par)
        call read_solvelist(unit_par, settings, iostat)
        rewind(unit_par)
        call read_physicslist(unit_par, settings, iostat)
        rewind(unit_par)
        call read_equilibriumlist(unit_par, settings, iostat)
        rewind(unit_par)
        call read_paramlist(unit_par, iostat)
        rewind(unit_par)
        call read_unitlist(unit_par, settings, iostat)
        rewind(unit_par)
        exit
      end do
    close(unit_par)

    if (settings%solvers%get_solver() == "none") call settings%io%set_all_io_to_false()
  end subroutine read_parfile


  subroutine read_gridlist(unit, settings, iostat)
    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    character(len=str_len) :: geometry
    integer :: gridpoints
    real(dp) :: x_start, x_end
    logical :: coaxial
    logical :: force_r0

    namelist /gridlist/ &
      geometry, gridpoints, x_start, x_end, coaxial, force_r0

    ! defaults
    geometry = ""
    gridpoints = 0
    x_start = 0.0_dp
    x_end = 1.0_dp
    coaxial = .false.
    force_r0 = .false.

    read(unit, nml=gridlist, iostat=iostat)

    call settings%grid%set_geometry(geometry)
    call settings%grid%set_gridpts(gridpoints)
    call settings%grid%set_grid_boundaries(grid_start=x_start, grid_end=x_end)
    settings%grid%coaxial = coaxial
    settings%grid%force_r0 = force_r0
  end subroutine read_gridlist


  subroutine read_savelist(unit, settings, iostat)
    use mod_global_variables, only: logging_level

    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    logical :: write_matrices, write_eigenvectors, write_residuals
    logical :: write_eigenfunctions, write_derived_eigenfunctions
    logical :: write_eigenfunction_subset
    logical :: show_results
    real(dp) :: eigenfunction_subset_radius
    complex(dp) :: eigenfunction_subset_center
    character(len=str_len) :: basename_datfile, output_folder

    namelist /savelist/ &
      write_matrices, write_eigenvectors, write_residuals, write_eigenfunctions, &
      write_derived_eigenfunctions, write_eigenfunction_subset, &
      show_results, basename_datfile, output_folder, logging_level, &
      eigenfunction_subset_radius, eigenfunction_subset_center

    ! defaults
    write_matrices = .false.
    write_eigenvectors = .false.
    write_residuals = .false.
    write_eigenfunctions = .true.
    write_derived_eigenfunctions = .false.
    write_eigenfunction_subset = .false.
    show_results = .true.
    basename_datfile = "datfile"
    output_folder = "output"
    logging_level = 2
    eigenfunction_subset_radius = NaN
    eigenfunction_subset_center = cmplx(NaN, NaN, kind=dp)

    read(unit, nml=savelist, iostat=iostat)

    settings%io%write_matrices = write_matrices
    settings%io%write_eigenvectors = write_eigenvectors
    settings%io%write_residuals = write_residuals
    settings%io%write_eigenfunctions = write_eigenfunctions
    settings%io%write_derived_eigenfunctions = write_derived_eigenfunctions
    settings%io%write_ef_subset = write_eigenfunction_subset
    if (write_eigenfunction_subset) call check_eigenfunction_subset_params( &
      center=eigenfunction_subset_center, radius=eigenfunction_subset_radius &
    )
    settings%io%ef_subset_radius = eigenfunction_subset_radius
    settings%io%ef_subset_center = eigenfunction_subset_center
    settings%io%show_results = show_results
    call settings%io%set_basename_datfile(basename_datfile)
    call settings%io%set_output_folder(output_folder)
  end subroutine read_savelist


  subroutine read_solvelist(unit, settings, iostat)
    use mod_global_variables, only: dp_LIMIT

    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    character(len=str_len) :: solver, arpack_mode
    character(len=2) :: which_eigenvalues
    integer :: number_of_eigenvalues, maxiter, ncv
    real(dp) :: tolerance
    complex(dp) :: sigma

    namelist /solvelist/ &
      solver, arpack_mode, which_eigenvalues, number_of_eigenvalues, &
      maxiter, ncv, tolerance, sigma

    ! defaults
    solver = "QR-cholesky"
    arpack_mode = "standard"
    which_eigenvalues = "LM"
    number_of_eigenvalues = 10
    maxiter = 0
    ncv = 0
    tolerance = dp_LIMIT
    sigma = (0.0_dp, 0.0_dp)

    read(unit, nml=solvelist, iostat=iostat)

    call settings%solvers%set_solver(solver)
    call settings%solvers%set_arpack_mode(arpack_mode)
    settings%solvers%number_of_eigenvalues = number_of_eigenvalues
    settings%solvers%which_eigenvalues = which_eigenvalues
    settings%solvers%maxiter = maxiter
    settings%solvers%ncv = ncv
    settings%solvers%tolerance = tolerance
    settings%solvers%sigma = sigma
  end subroutine read_solvelist


  subroutine read_physicslist(unit, settings, iostat)
    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    logical :: flow, incompressible, radiative_cooling, external_gravity, &
      parallel_conduction, perpendicular_conduction, use_fixed_tc_para, &
      use_fixed_tc_perp, resistivity, use_fixed_resistivity, use_eta_dropoff, &
      viscosity, viscous_heating, hall_mhd, hall_substitution, hall_dropoff, &
      elec_inertia, inertia_dropoff
    integer :: ncool
    character(len=str_len) :: cooling_curve
    real(dp) :: mhd_gamma
    real(dp) :: fixed_tc_para_value, fixed_tc_perp_value, fixed_resistivity_value
    real(dp) :: viscosity_value
    real(dp) :: electron_fraction
    real(dp) :: dropoff_edge_dist, dropoff_width

    namelist /physicslist/ &
      mhd_gamma, flow, incompressible, radiative_cooling, external_gravity, parallel_conduction, perpendicular_conduction, use_fixed_tc_para, &
      use_fixed_tc_perp, resistivity, use_fixed_resistivity, use_eta_dropoff, &
      viscosity, viscous_heating, hall_mhd, hall_substitution, hall_dropoff, &
      elec_inertia, inertia_dropoff, ncool, cooling_curve, fixed_tc_para_value, &
      fixed_tc_perp_value, fixed_resistivity_value, dropoff_edge_dist, dropoff_width, &
      viscosity_value, electron_fraction, dropoff_edge_dist, dropoff_width

    ! defaults
    mhd_gamma = 5.0_dp / 3.0_dp
    flow = .false.
    incompressible = .false.
    radiative_cooling = .false.
    external_gravity = .false.
    parallel_conduction = .false.
    perpendicular_conduction = .false.
    use_fixed_tc_para = .false.
    use_fixed_tc_perp = .false.
    resistivity = .false.
    use_fixed_resistivity = .false.
    use_eta_dropoff = .false.
    viscosity = .false.
    viscous_heating = .false.
    hall_mhd = .false.
    hall_substitution = .true.
    hall_dropoff = .false.
    elec_inertia = .false.
    inertia_dropoff = .false.

    ncool = 4000
    cooling_curve = "jc_corona"
    fixed_tc_para_value = 0.0_dp
    fixed_tc_perp_value = 0.0_dp
    fixed_resistivity_value = 0.0_dp
    dropoff_edge_dist = 0.0_dp
    dropoff_width = 0.0_dp
    viscosity_value = 0.0_dp
    electron_fraction = 0.5_dp
    dropoff_edge_dist = 0.05_dp
    dropoff_width = 0.1_dp

    read(unit, nml=physicslist, iostat=iostat)

    call settings%physics%set_gamma(mhd_gamma)
    if (incompressible) call settings%physics%set_incompressible()
    if (flow) call settings%physics%flow%enable()
    if (radiative_cooling) call settings%physics%enable_cooling(cooling_curve, ncool)
    if (external_gravity) call settings%physics%enable_gravity()
    if (parallel_conduction) call settings%physics%enable_parallel_conduction( &
      use_fixed_tc_para, fixed_tc_para_value &
    )
    if (perpendicular_conduction) then
      call settings%physics%enable_perpendicular_conduction( &
        use_fixed_tc_perp, fixed_tc_perp_value &
      )
    end if
    if (resistivity) call settings%physics%enable_resistivity( &
      use_fixed_resistivity, fixed_resistivity_value &
    )
    if (viscosity) call settings%physics%enable_viscosity( &
      viscosity_value, viscous_heating &
    )
    if (hall_mhd) call settings%physics%enable_hall( &
      hall_substitution, elec_inertia, electron_fraction &
    )
    settings%physics%dropoff_edge_dist = dropoff_edge_dist
    settings%physics%dropoff_width = dropoff_width
  end subroutine read_physicslist


  subroutine read_equilibriumlist(unit, settings, iostat)
    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    character(len=str_len) :: equilibrium_type, boundary_type
    logical :: use_defaults

    namelist /equilibriumlist/ equilibrium_type, boundary_type, use_defaults

    ! defaults
    equilibrium_type = "adiabatic_homo"
    boundary_type = "wall"
    use_defaults = .true.

    read(unit, nml=equilibriumlist, iostat=iostat)

    call settings%equilibrium%set_equilibrium_type(equilibrium_type)
    call settings%equilibrium%set_boundary_type(boundary_type)
    settings%equilibrium%use_defaults = use_defaults
  end subroutine read_equilibriumlist


  subroutine read_paramlist(unit, iostat)
    use mod_equilibrium_params

    integer, intent(in) :: unit
    integer, intent(out) :: iostat

    namelist /paramlist/  &
      k2, k3, cte_rho0, cte_T0, cte_B01, cte_B02, cte_B03, cte_v02, cte_v03, &
      cte_p0, p1, p2, p3, p4, p5, p6, p7, p8, &
      alpha, beta, delta, theta, tau, lambda, nu, &
      r0, rc, rj, Bth0, Bz0, V, j0, g, eq_bool

    call init_equilibrium_params()
    read(unit, nml=paramlist, iostat=iostat)
  end subroutine read_paramlist


  subroutine read_unitlist(unit, settings, iostat)
    use mod_global_variables, only: cgs_units

    integer, intent(in) :: unit
    type(settings_t), intent(inout) :: settings
    integer, intent(out) :: iostat

    real(dp) :: unit_density, unit_temperature, unit_magneticfield, unit_length
    real(dp) :: mean_molecular_weight

    namelist /unitslist/ &
      cgs_units, unit_density, unit_temperature, unit_magneticfield, unit_length, &
      mean_molecular_weight

    ! defaults
    unit_density = NaN
    unit_temperature = NaN
    unit_magneticfield = NaN
    unit_length = NaN
    mean_molecular_weight = NaN

    read(unit, nml=unitslist, iostat=iostat)

    call check_and_set_supplied_unit_normalisations( &
      unit_density, &
      unit_temperature, &
      unit_magneticfield, &
      unit_length, &
      mean_molecular_weight &
    )
  end subroutine read_unitlist


  !> Checks the unit normalisations that are supplied (if any), sets the
  !! unit normalisations if valid.
  subroutine check_and_set_supplied_unit_normalisations( &
    unit_density, &
    unit_temperature, &
    unit_magneticfield, &
    unit_length, &
    mean_molecular_weight &
  )
    use mod_units, only: set_normalisations
    use mod_check_values, only: is_NaN

    !> supplied unit density
    real(dp), intent(in)  :: unit_density
    !> supplied unit temperature
    real(dp), intent(in)  :: unit_temperature
    !> supplied unit magneticfield
    real(dp), intent(in)  :: unit_magneticfield
    !> supplied unit length
    real(dp), intent(in)  :: unit_length
    !> supplied mean molecular weight
    real(dp), intent(in)  :: mean_molecular_weight

    ! Provide normalisations, if supplied
    if (.not. is_NaN(unit_density) .and. .not. is_NaN(unit_temperature)) then
      call log_message( &
        "unit density and unit temperature cannot both be provided in the parfile!", &
        level="error" &
      )
      return
    end if

    ! set normalisations if density is supplied
    if (.not. is_NaN(unit_density)) then
      if (is_NaN(unit_magneticfield) .or. is_NaN(unit_length)) then
        call log_message( &
          "unit_density found, unit_magneticfield and unit_length also required.", &
          level="error" &
        )
        return
      end if
      if (is_NaN(mean_molecular_weight)) then
        call set_normalisations( &
          new_unit_density=unit_density, &
          new_unit_magneticfield=unit_magneticfield, &
          new_unit_length=unit_length &
        )
      else
        call set_normalisations( &
          new_unit_density=unit_density, &
          new_unit_magneticfield=unit_magneticfield, &
          new_unit_length=unit_length, &
          new_mean_molecular_weight=mean_molecular_weight &
        )
      end if
    end if

    ! set normalisations if temperature is supplied
    if (.not. is_NaN(unit_temperature)) then
      if (is_NaN(unit_magneticfield) .or. is_NaN(unit_length)) then
        call log_message( &
          "unit_temperature found, unit_magneticfield and unit_length also required.", &
          level="error" &
        )
        return
      end if
      if (is_NaN(mean_molecular_weight)) then
        call set_normalisations( &
          new_unit_temperature=unit_temperature, &
          new_unit_magneticfield=unit_magneticfield, &
          new_unit_length=unit_length &
        )
      else
        call set_normalisations( &
          new_unit_temperature=unit_temperature, &
          new_unit_magneticfield=unit_magneticfield, &
          new_unit_length=unit_length, &
          new_mean_molecular_weight=mean_molecular_weight &
        )
      end if
    end if
  end subroutine check_and_set_supplied_unit_normalisations


  !> Called when the eigenfunction subset selection is enabled, this checks if the
  !! global variables are properly set.
  subroutine check_eigenfunction_subset_params(center, radius)
    use mod_check_values, only: is_NaN

    complex(dp), intent(in) :: center
    real(dp), intent(in) :: radius

    if (is_NaN(center)) then
      call log_message("eigenfunction_subset_center must be set!", level="error")
      return
    end if
    if (is_NaN(radius)) then
      call log_message("eigenfunction_subset_radius must be set!", level="error")
      return
    end if
  end subroutine check_eigenfunction_subset_params


  ! LCOV_EXCL_START <this routine is excluded from coverage>
  !> Parses the command line arguments and retrieves the parfile passed.
  !! @warning Throws an error if
  !!
  !! - command line arguments can not be parsed.
  !! - the parfile is not found. @endwarning
  !! @note If no parfile is passed, the code uses a default configuration. @endnote
  subroutine get_parfile(filename_par)
    !> the name of the parfile
    character(len=*), intent(out) :: filename_par

    integer :: num_args, i
    character(len=5*str_len), allocatable :: args(:)
    logical :: file_exists

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
        call log_message("unable to read in command line arguments", level='error')
      end select
    end do

    ! If no parfile supplied, flag an error
    if (filename_par == "") then
      call log_message("no parfile supplied, please provide the file path with the -i flag", level="error")
    end if

    !! Check if supplied file exists
    inquire(file=trim(filename_par), exist=file_exists)
    if (.not. file_exists) then
      call log_message(("parfile not found: " // trim(filename_par)), level='error')
    end if
  end subroutine get_parfile
  ! LCOV_EXCL_STOP

end module mod_input
