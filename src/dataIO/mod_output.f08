! =============================================================================
!> This module contains all routines for file opening and file writing.
module mod_output
  use mod_global_variables, only: dp, str_len
  implicit none

  private

  ! IO units -- do not use 0/5/6/7, these are system-reserved
  !> filehandler IO unit for the main data file
  integer, parameter  :: dat_fh = 10
  !> filehandler IO unit for the log file
  integer, parameter  :: log_fh = 20
  !> datfile name
  character(len=str_len) :: datfile_name
  !> logfile name
  character(len=str_len) :: logfile_name


  public :: create_datfile
  public :: datfile_name

contains

  !> Opens a file with a given IO unit and filename.
  !! All files are opened using <tt>access='stream'</tt>,
  !! <tt>status='unknown'</tt> and <tt>action='write'</tt>.
  subroutine open_file(file_unit, filename)
    !> the IO unit of the file to open
    integer, intent(in)           :: file_unit
    !> the filename of the file to open
    character(len=*), intent(in)  :: filename

    open(unit=file_unit, file=filename, access='stream', &
         status='unknown', action='write')
  end subroutine open_file


  !> Builds a filename based on a given base filename and the
  !! output folder defined in the global variables module.
  !! The output folder is prepended to the base filename.
  !! @note    At this point filenames are not yet given extensions.
  subroutine make_filename(base_filename, filename)
    use mod_global_variables, only: output_folder

    !> the base filename to use
    character(len=*), intent(in)  :: base_filename
    !> the filename that is created
    character(len=*), intent(out) :: filename

    filename = trim(trim(output_folder) // "/" // base_filename)
  end subroutine make_filename


  !> Writes the datfile, where eigenfunctions and matrices are optionally
  !! included. First a header is written containing default information
  !! on the configuration, then the actual data.
  !! @note    Eigenfunctions are only written if this is enabled in the
  !!          global variables. @endnote
  !! @note    Matrices are only written if this is enabled in the global variables.
  !!          The matrices are not written entirely to save diskspace: a
  !!          first pass is performed locating the non-zero values, and then the
  !!          values are saved to file in the format
  !!          <tt>(row_idx, column_idx, value)</tt>. @endnote
  !! @note    The extension <tt>".dat"</tt> is appended to the filename. @endnote
  subroutine create_datfile(eigenvalues, matrix_A, matrix_B)
    use mod_global_variables
    use mod_version, only: LEGOLAS_VERSION
    use mod_logging, only: log_message
    use mod_grid, only: grid, grid_gauss
    use mod_equilibrium, only: rho_field, T_field, B_field, v_field, rc_field, &
      kappa_field, eta_field, grav_field, hall_field
    use mod_eigenfunctions, only: ef_grid, ef_names, ef_array
    use mod_check_values, only: value_is_zero
    use mod_equilibrium_params
    use mod_units

    !> the eigenvalues
    complex(dp), intent(in)       :: eigenvalues(:)
    !> the A-matrix
    complex(dp), intent(in)       :: matrix_A(:, :)
    !> the B-matrix
    real(dp), intent(in)          :: matrix_B(:, :)

    character(len=str_len_arr)    :: param_names(33), equil_names(27)
    character(len=2*str_len_arr)  :: unit_names(11)
    integer                       :: i, j, nonzero_A_values, nonzero_B_values

    param_names = [ &
      character(len=str_len_arr) :: "k2", "k3", "cte_rho0", "cte_T0", "cte_B01", &
      "cte_B02", "cte_B03", "cte_v02", "cte_v03", "cte_p0", "p1", "p2", "p3", &
      "p4", "p5", "p6", "p7", "p8", "alpha", "beta", "delta", "theta", "tau", &
      "lambda", "nu", "r0", "rc", "rj", "Bth0", "Bz0", "V", "j0", "g" &
    ]
    equil_names = [ &
      character(len=str_len_arr) :: &
        "rho0", "drho0", &
        "T0", "dT0", &
        "B02", "B03", "dB02", "dB03", "ddB02", "ddB03", "B0", &
        "v02", "v03", "dv02", "dv03", "ddv02", "ddv03", &
        "dLdT", "dLdrho", &
        "kappa_para", "kappa_perp", &
        "eta", "detadT", "detadr", &
        "grav", &
        "Hall", "inertia" &
    ]
    unit_names = [ &
      character(len=2*str_len_arr) :: "unit_length", "unit_time", "unit_density", &
      "unit_velocity", "unit_temperature", "unit_pressure", "unit_magneticfield", &
      "unit_numberdensity", "unit_lambdaT", "unit_conduction", "unit_resistivity" &
    ]

    call make_filename(trim(basename_datfile) // ".dat", datfile_name)
    call open_file(dat_fh, datfile_name)

    ! First we write all header information
    write(dat_fh) "legolas_version", LEGOLAS_VERSION
    write(dat_fh) str_len, str_len_arr, geometry, x_start, x_end, gridpts, &
      gauss_gridpts, matrix_gridpts, ef_gridpts, gamma, equilibrium_type, &
      write_eigenfunctions, write_matrices
    write(dat_fh) size(param_names), len(param_names(1)), param_names
    write(dat_fh) k2, k3, cte_rho0, cte_T0, cte_B01, cte_B02, cte_B03, cte_v02, &
      cte_v03, cte_p0, p1, p2, p3, p4, p5, p6, p7, p8, alpha, beta, delta, &
      theta, tau, lambda, nu, r0, rc, rj, Bth0, Bz0, V, j0, g
    write(dat_fh) size(equil_names), len(equil_names(1)), equil_names
    write(dat_fh) cgs_units
    write(dat_fh) size(unit_names), len(unit_names(1)), unit_names
    write(dat_fh) unit_length, unit_time, unit_density, unit_velocity, &
      unit_temperature, unit_pressure, unit_magneticfield, unit_numberdensity, &
      unit_lambdaT, unit_conduction, unit_resistivity

    ! Next write the data itself
    ! General data: eigenvalues, grids, equilibrium configuration
    write(dat_fh) size(eigenvalues), eigenvalues, grid, grid_gauss
    write(dat_fh) &
      rho_field % rho0, rho_field % d_rho0_dr, &
      T_field % T0, T_field % d_T0_dr, &
      B_field % B02, B_field % B03, &
      B_field % d_B02_dr, B_field % d_B03_dr, &
      eta_field % dd_B02_dr, eta_field % dd_B03_dr, B_field % B0, &
      v_field % v02, v_field % v03, &
      v_field % d_v02_dr, v_field % d_v03_dr, &
      v_field % dd_v02_dr, v_field % dd_v03_dr, &
      rc_field % d_L_dT, rc_field % d_L_drho, &
      kappa_field % kappa_para, kappa_field % kappa_perp, &
      eta_field % eta, eta_field % d_eta_dT, eta_field % d_eta_dr, &
      grav_field % grav, &
      hall_field % hallfactor, hall_field % inertiafactor

    ! Eigenfunction data [optional]
    if (write_eigenfunctions) then
      call log_message("writing eigenfunctions...", level='info')
      write(dat_fh) size(ef_names), ef_names
      write(dat_fh) ef_grid
      do i = 1, nb_eqs
        write(dat_fh) ef_array(i) % eigenfunctions
      end do
    end if

    ! Matrix data [optional]
    if (write_matrices) then
      call log_message("writing matrices...", level='info')
      ! Write non-zero matrix indices and values. Since this varies every run, we first
      ! loop through without writing and count the non-zero values. This number is needed
      ! to correctly read in the values later on (we have to know how many there are).
      nonzero_B_values = 0
      nonzero_A_values = 0
      do j = 1, matrix_gridpts
        do i = 1, matrix_gridpts
          if (.not. value_is_zero(matrix_B(i, j))) then
            nonzero_B_values = nonzero_B_values + 1
          end if
          if (.not. value_is_zero(matrix_A(i, j))) then
            nonzero_A_values = nonzero_A_values + 1
          end if
        end do
      end do
      ! write these numbers to the file
      write(dat_fh) nonzero_B_values, nonzero_A_values
      do j = 1, matrix_gridpts
        do i = 1, matrix_gridpts
          if (.not. value_is_zero(matrix_B(i, j))) then
            write(dat_fh) i, j
            write(dat_fh) matrix_B(i, j)
          end if
        end do
      end do
      ! Write non-zero A matrix indices and values
      do j = 1, matrix_gridpts
        do i = 1, matrix_gridpts
          if (.not. value_is_zero(matrix_A(i, j))) then
            write(dat_fh) i, j
            write(dat_fh) matrix_A(i, j)
          end if
        end do
      end do
    end if

    call log_message("results saved to " // trim(datfile_name), level='info')
    close(dat_fh)

    call create_logfile(eigenvalues)
  end subroutine create_datfile


  !> Creates a logfile. If <tt>basename_logfile</tt> is specified in the datfile,
  !! a logfile is written.
  !! This is a pure textfile containing the real and imaginary parts of the
  !! eigenvalues, written in an exponential format. This is mainly used for testing
  !! purposes but may come in handy to do some quick inspections on the data.
  !! @note    If <tt>basename_logfile</tt> is unspecified in the parfile,
  !!          no logfile is written. @endnote
  !! @note    The extension <tt>".log"</tt> is appended to the filename. @endnote
  subroutine create_logfile(eigenvalues)
    use mod_global_variables, only: basename_logfile
    use mod_logging, only: log_message, exp_fmt

    !> the eigenvalues
    complex(dp), intent(in)   :: eigenvalues(:)
    character(20)             :: real_part, imag_part
    integer   :: i

    if (basename_logfile == "") then
      return
    end if

    call make_filename(trim(basename_logfile) // ".log", logfile_name)
    ! open manually since this is not a binary file
    open(unit=log_fh, file=logfile_name, status='unknown', action='write')
    do i = 1, size(eigenvalues)
      write(real_part, exp_fmt) real(eigenvalues(i))
      write(imag_part, exp_fmt) aimag(eigenvalues(i))
      write(log_fh, *) real_part, ',', imag_part
    end do

    call log_message("eigenvalues logged to " // trim(logfile_name), level='info')
    close(log_fh)
  end subroutine create_logfile

end module mod_output
