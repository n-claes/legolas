module mod_input
  use mod_global_variables
  use mod_equilibrium_params
  implicit none

  private

  integer       :: unit_par = 101

  public :: read_parfile
  public :: get_parfile

contains

  !> Reads in the supplied parfile and sets all global variables accordingly.
  !! @param[in] parfile   The name of the parfile
  subroutine read_parfile(parfile)
    use mod_check_values, only: value_is_zero, value_is_nan
    use mod_units, only: set_normalisations

    character(len=*), intent(in)  :: parfile

    real(dp)    :: mhd_gamma
    real(dp)    :: unit_density, unit_temperature, unit_magneticfield, unit_length
    integer     :: gridpoints

    namelist /physicslist/  mhd_gamma, flow, radiative_cooling, ncool, cooling_curve, &
                            external_gravity, thermal_conduction, use_fixed_tc_para, fixed_tc_para_value, &
                            use_fixed_tc_perp, fixed_tc_perp_value, resistivity, use_fixed_resistivity, fixed_eta_value
    namelist /unitslist/    cgs_units, unit_density, unit_temperature, unit_magneticfield, unit_length
    namelist /gridlist/     geometry, x_start, x_end, gridpoints, mesh_accumulation, ev_1, ev_2, sigma_1, sigma_2
    namelist /equilibriumlist/ equilibrium_type, boundary_type, use_defaults
    namelist /savelist/     run_silent, write_matrices, write_eigenvectors, write_eigenfunctions, &
                            write_equilibrium, show_results, show_matrices, show_eigenfunctions, show_equilibrium
    namelist /filelist/     savename_eigenvalues, savename_efgrid, savename_matrix, savename_eigenvectors, &
                            savename_eigenfunctions, savename_config, savename_equil
    namelist /paramlist/    k2, k3, cte_rho0, cte_T0, cte_B02, cte_B03, cte_v02, cte_v03, cte_p0, &
                            p1, p2, p3, p4, p5, p6, p7, p8, &
                            alpha, beta, delta, theta, tau, lambda, nu, &
                            r0, rc, rj, Bth0, Bz0, V, j0, g

    !! Initialise equilibrium parameters to NaN. These are then read in by the paramlist
    !! or set directly in the submodules.
    call init_equilibrium_params()

    ! if no parfile supplied, return to keep using defaults
    if (parfile == "") then
      return
    end if

    ! initialise local variables
    mhd_gamma = 0.0d0
    gridpoints = 0
    unit_density = NaN
    unit_temperature = NaN
    unit_magneticfield = NaN
    unit_length = NaN

    open(unit_par, file=trim(parfile), status='old')
    !! Start reading namelists, rewind so they can appear out of order
          rewind(unit_par)
          read(unit_par, gridlist, end=1001)

    1001  rewind(unit_par)
          read(unit_par, physicslist, end=1002)

    1002  rewind(unit_par)
          read(unit_par, equilibriumlist, end=1003)

    1003  rewind(unit_par)
          read(unit_par, savelist, end=1004)

    1004  rewind(unit_par)
          read(unit_par, filelist, end=1005)

    1005  rewind(unit_par)
          read(unit_par, paramlist, end=1006)

    1006  rewind(unit_par)
          read(unit_par, unitslist, end=1007)

    1007  close(unit_par)

    !> Set gridpoints and gamma, if supplied
    if (.not. gridpoints == 0) then
      call set_gridpts(gridpoints)
    end if
    if (.not. value_is_zero(mhd_gamma)) then
      call set_gamma(mhd_gamma)
    end if

    !> Provide normalisations, if supplied
    if (.not. value_is_nan(unit_density) .and. .not. value_is_nan(unit_temperature)) then
      error stop "unit density and unit temperature can not both be provided in the par file! (choose one)"
    end if
    if (.not. value_is_nan(unit_density)) then
      if (value_is_nan(unit_magneticfield) .or. value_is_nan(unit_length)) then
        error stop "Unit_density found, but unit_magneticfield and unit_length are also required."
      end if
      call set_normalisations(new_unit_density=unit_density, new_unit_magneticfield=unit_magneticfield, &
                              new_unit_length=unit_length)
    else if (.not. value_is_nan(unit_temperature)) then
      if (value_is_nan(unit_magneticfield) .or. value_is_nan(unit_length)) then
        error stop "Unit_density found, but unit_magneticfield and unit_length are also required."
      end if
      call set_normalisations(new_unit_temperature=unit_temperature, &
                              new_unit_magneticfield=unit_magneticfield, new_unit_length=unit_length)
    end if

  end subroutine read_parfile


  !> Retrieves the parfile passed as command line argument.
  !! @param[out] filename_par   The name of the parfile
  subroutine get_parfile(filename_par)
    character(len=str_len), intent(out) :: filename_par

    integer                             :: num_args, i
    character(len=str_len), allocatable :: args(:)
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
