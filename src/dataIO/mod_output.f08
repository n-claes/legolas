module mod_output
  use mod_global_variables, only: dp, str_len
  implicit none

  private

  !! IO units -- do not use 0/5/6/7, these are system-reserved
  !> Filehandler IO unit for the main data file
  integer, parameter  :: dat_fh = 10

  !! Format settings
  !> long exponential format
  character(8), parameter    :: form_e = '(e30.20)'
  !> long float format
  character(8), parameter    :: form_f = '(f30.20)'
  !> shorter exponential format
  character(8), parameter    :: form_eout = '(e20.10)'
  !> shorter float format
  character(8), parameter    :: form_fout = '(f20.10)'
  !> integer format
  character(4), parameter    :: form_int  = '(i8)'

  !> name of base output folder
  character(len=7)    :: output_folder = 'output/'
  !> filename extension
  character(len=4)    :: file_extension = '.dat'
  !> datfile name
  character(len=str_len) :: datfile_name

  public :: create_datfile
  public :: startup_info_toconsole
  public :: datfile_name

contains

  !> Opens a file with given unit and filename.
  !! @param[in] file_unit   IO unit for opening the file
  !! @param[in] filename    filename to open
  subroutine open_file(file_unit, filename)
    integer, intent(in)           :: file_unit
    character(len=*), intent(in)  :: filename

    open(unit=file_unit, file=filename, access='stream', &
         status='unknown', action='write')
  end subroutine open_file


  !> Creates a filename, prepending the output directory and appending the
  !! file extension
  !! @param[in] base_filename   the base filename to use
  !! @param[out] filename       the resulting filename, of the form
  !!                            output_folder/base_filename.extension
  subroutine make_filename(base_filename, filename)
    character(len=*), intent(in)  :: base_filename
    character(len=*), intent(out) :: filename

    filename = trim(trim(output_folder) // trim(base_filename) // &
                    trim(file_extension))
  end subroutine make_filename


  subroutine create_datfile(base_filename, eigenvalues, matrix_A, matrix_B)
    use mod_global_variables, only: geometry, x_start, x_end, gridpts, gauss_gridpts, &
                                    matrix_gridpts, ef_gridpts, gamma, equilibrium_type, &
                                    nb_eqs, cgs_units, str_len_arr, run_silent, &
                                    write_matrices, write_eigenfunctions
    use mod_grid, only: grid, grid_gauss
    use mod_equilibrium, only: rho_field, T_field, B_field, v_field, rc_field, kappa_field, eta_field, grav_field
    use mod_eigenfunctions, only: ef_grid, ef_names, ef_array
    use mod_check_values, only: value_is_zero
    use mod_equilibrium_params
    use mod_units

    character(len=*), intent(in)  :: base_filename
    complex(dp), intent(in)       :: eigenvalues(matrix_gridpts)
    complex(dp), intent(in)       :: matrix_A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)          :: matrix_B(matrix_gridpts, matrix_gridpts)

    character(len=16)             :: param_names(32), equil_names(20)
    integer                       :: i, j, nonzero_A_values, nonzero_B_values

    param_names = [character(len=str_len_arr) :: 'k2', 'k3', 'cte_rho0', 'cte_T0', 'cte_B02', 'cte_B03', &
                   'cte_v02', 'cte_v03', 'cte_p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', &
                   'alpha', 'beta', 'delta', 'theta', 'tau', 'lambda', 'nu', 'r0', 'rc', 'rj', &
                   'Bth0', 'Bz0', 'V', 'j0', 'g']
    equil_names = [character(len=str_len_arr) :: 'rho0', 'drho0', 'T0', 'dT0', 'B02', 'B03', 'dB02', 'dB03', &
                   'B0', 'v02', 'v03', 'dv02', 'dv03', 'dLdT', 'dLdrho', 'kappa_para', 'kappa_perp', &
                   'eta', 'detadT', 'grav']

    call make_filename(base_filename, datfile_name)
    call open_file(dat_fh, datfile_name)

    ! First we write all header information
    write(dat_fh) str_len, str_len_arr, geometry, x_start, x_end, gridpts, gauss_gridpts, matrix_gridpts, &
                  ef_gridpts, gamma, equilibrium_type, write_eigenfunctions, write_matrices
    write(dat_fh) size(param_names), param_names
    write(dat_fh) k2, k3, cte_rho0, cte_T0, cte_B02, cte_B03, cte_v02, cte_v03, cte_p0, &
                  p1, p2, p3, p4, p5, p6, p7, p8, alpha, beta, delta, theta, tau, lambda, nu, &
                  r0, rc, rj, Bth0, Bz0, V, j0, g
    write(dat_fh) size(equil_names), equil_names
    write(dat_fh) cgs_units
    write(dat_fh) unit_length, unit_time, unit_density, unit_velocity, unit_temperature, &
                  unit_pressure, unit_magneticfield, unit_numberdensity, unit_luminosity, &
                  unit_conduction, unit_resistivity

    ! Next write the data itself
    ! General data: eigenvalues, grids, equilibrium configuration
    write(dat_fh) eigenvalues, grid, grid_gauss
    write(dat_fh) rho_field % rho0, rho_field % d_rho0_dr, T_field % T0, T_field % d_T0_dr, &
                  B_field % B02, B_field % B03, B_field % d_B02_dr, B_field % d_B03_dr, &
                  B_field % B0, v_field % v02, v_field % v03, v_field % d_v02_dr, v_field % d_v03_dr, &
                  rc_field % d_L_dT, rc_field % d_L_drho, kappa_field % kappa_para, kappa_field % kappa_perp, &
                  eta_field % eta, eta_field % d_eta_dT, grav_field % grav

    ! Eigenfunction data [optional]
    if (write_eigenfunctions) then
      if (.not. run_silent) then
        write(*, *) "Writing eigenfunctions..."
      end if

      write(dat_fh) size(ef_names), ef_names
      write(dat_fh) ef_grid
      do i = 1, nb_eqs
        write(dat_fh) ef_array(i) % eigenfunctions
      end do
    end if

    ! Matrix data [optional]
    if (write_matrices) then
      if (.not. run_silent) then
        write(*, *) "Writing matrices..."
      end if

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

    if (.not. run_silent) then
      write(*, *) "Results saved to ", trim(datfile_name)
    end if

    close(dat_fh)
  end subroutine create_datfile


  !> Converts a Fortran logical ('T' or 'F') to a string ('true', 'false').
  !! @param[in] boolean   Fortran logical to convert
  !! @param[out] boolean_string   'true' if boolean == True, 'false' otherwise
  subroutine logical_tostring(boolean, boolean_string)
    logical, intent(in)             :: boolean
    character(len=20), intent(out)  :: boolean_string

    if (boolean) then
      boolean_string = 'true'
    else
      boolean_string = 'false'
    end if
  end subroutine logical_tostring


  !> Prints basic information of the current configuration to the console.
  subroutine startup_info_toconsole()
    use mod_global_variables, only: geometry, x_start, x_end, gridpts, gauss_gridpts, matrix_gridpts, &
                                    equilibrium_type, boundary_type, gamma, write_matrices, write_eigenfunctions
    use mod_equilibrium_params, only: k2, k3

    character(20)                   :: char

    write(*, *) ""
    write(*, *) "------------------------------"
    write(*, *) "----------- LEGOLAS ----------"
    write(*, *) "------------------------------"
    write(*, *) ""

    write(*, *) "Running with the following configuration:"
    write(*, *) ""

    ! Geometry info
    write(*, *) "-- Geometry settings --"
    write(*, *) "Coordinate system  : ", geometry
    write(char, form_fout) x_start
    write(*, *) "Start              : ", adjustl(char)
    write(char, form_fout) x_end
    write(*, *) "End                : ", adjustl(char)
    write(char, form_int) gridpts
    write(*, *) "Gridpoints         : ", adjustl(char)
    write(char, form_int) gauss_gridpts
    write(*, *) "Gaussian gridpoints: ", adjustl(char)
    write(char, form_int) matrix_gridpts
    write(*, *) "Matrix gridpoints  : ", adjustl(char)
    write(*, *) ""

    ! Equilibrium info
    write(*, *) "-- Equilibrium settings --"
    write(*, *) "Equilibrium type   : ", equilibrium_type
    write(*, *) "Boundary conditions: ", boundary_type
    write(char, form_fout) gamma
    write(*, *) "Gamma              : ", adjustl(char)
    write(char, form_fout) k2
    write(*, *) "Wave number k2     : ", adjustl(char)
    write(char, form_fout) k3
    write(*, *) "Wave number k3     : ", adjustl(char)
    write(*, *) ""

    ! Save info
    write(*, *) "-- DataIO settings --"
    call logical_tostring(write_matrices, char)
    write(*, *) "Write matrices to file       : ", char
    call logical_tostring(write_eigenfunctions, char)
    write(*, *) "Write eigenfunctions to file : ", char

    write(*, *) '----------------------------------------------------'
    write(*, *) ''

  end subroutine startup_info_toconsole

end module mod_output
