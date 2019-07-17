module mod_io
  use mod_global_variables
  use mod_physical_constants
  implicit none


  public

  integer, parameter  :: w_output            = 10
  integer, parameter  :: config              = 30
  integer, parameter  :: mat_a_out           = 40
  integer, parameter  :: mat_b_out           = 50
  integer, parameter  :: eigenf_l_out        = 60
  integer, parameter  :: eigenf_r_out        = 70

  character(8), parameter    :: form_e = '(e30.20)'
  character(8), parameter    :: form_f = '(f30.20)'

  character(8), parameter    :: form_eout = '(e20.10)'
  !character(10), parameter   :: form_eout_long = '(e20.10e3)'
  character(8), parameter    :: form_fout = '(f20.10)'

  character(4), parameter    :: form_int  = '(i8)'

contains

  subroutine open_file(file_no, filename, append)
    integer, intent(in)          :: file_no
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: append

    logical                      :: file_present

    if (append) then
      inquire(file=filename, exist=file_present)

      if (file_present) then
        open(file_no, file=filename, status='old', position='append', &
             action='write')
      else
        open(file_no, file=filename, status='new', action='write')
      end if
    else
      open(file_no, file=filename, status='unknown')
    end if
  end subroutine open_file



  subroutine save_eigenvalues(omega, filenameW, append)
    complex(dp), intent(in)      :: omega(matrix_gridpts)
    character(len=*), intent(in) :: filenameW
    logical, intent(in)          :: append
    integer                      :: i

    character(30)                :: w_char_r, w_char_i
    character(str_len)           :: filenameW_out

    filenameW_out = trim("output/" // trim(filenameW) // ".txt")

    call open_file(w_output, filenameW_out, append)

    do i = 1, matrix_gridpts
      ! Write normalised output
      write(w_char_r, form_e) real(omega(i))
      write(w_char_i, form_e) aimag(omega(i))

      write(w_output, *) w_char_r, ',', w_char_i
    end do

    close(w_output)

  end subroutine save_eigenvalues

  subroutine save_config(filenameCFG)
    character(len=*), intent(in)  :: filenameCFG

    character(20)      :: char
    character(str_len) :: filenameCFG_out

    filenameCFG_out = trim("output/" // trim(filenameCFG) // ".txt")

    call open_file(config, filenameCFG_out, .false.)

    write(config, *) "Equilibrium type   : ", equilibrium_type

    write(config, *) "Flow               : ", flow
    write(config, *) "Radiative cooling  : ", radiative_cooling
    if (radiative_cooling) then
      write(config, *) "Cooling curve      : ", cooling_curve
    end if
    write(config, *) "Thermal conduction : ", thermal_conduction
    write(config, *) "Resistivity        : ", resistivity
    write(config, *) "External gravity   : ", external_gravity
    if (external_gravity) then
      write(config, *) "Gravity strength   : ", gravity_type
    end if

    write(char, form_fout) gamma
    write(config, *) "Gamma              : ", adjustl(char)

    write(config, *) ""
    write(config, *) "Geometry           : ", geometry
    write(char, form_fout) x_start
    write(config, *) "Start              : ", adjustl(char)
    write(char, form_fout) x_end
    write(config, *) "End                : ", adjustl(char)
    write(char, form_int) gridpts
    write(config, *) "Gridpoints         : ", adjustl(char)
    write(char, form_int) matrix_gridpts
    write(config, *) "Matrix gridpoints  : ", adjustl(char)
    write(char, form_int) gauss_gridpts
    write(config, *) "Gaussian gridpoints: ", adjustl(char)
    write(config, *) "Mesh accumulation  : ", mesh_accumulation

    write(config, *) ""

    write(config, *) "CGS units          : ", cgs_units
    write(char, form_eout) unit_length
    write(config, *) "Unit length        : ", adjustl(char)
    write(char, form_eout) unit_time
    write(config, *) "Unit time          : ", adjustl(char)
    write(char, form_eout) unit_density
    write(config, *) "Unit density       : ", adjustl(char)
    write(char, form_eout) unit_velocity
    write(config, *) "Unit velocity      : ", adjustl(char)
    write(char, form_eout) unit_temperature
    write(config, *) "Unit temperature   : ", adjustl(char)
    write(char, form_eout) unit_pressure
    write(config, *) "Unit pressure      : ", adjustl(char)
    write(char, form_eout) unit_magneticfield
    write(config, *) "Unit magnetic field: ", adjustl(char)
    write(char, form_eout) unit_numberdensity
    write(config, *) "Unit numberdensity : ", adjustl(char)
    write(char, form_eout) unit_luminosity
    write(config, *) "Unit luminosity    : ", adjustl(char)
    write(char, form_eout) unit_resistivity
    write(config, *) "Unit resistivity   : ", adjustl(char)

    close(config)
  end subroutine save_config

  subroutine save_matrices(matrix_A, matrix_B, filenameA, filenameB)
    complex(dp), intent(in)   :: matrix_A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)      :: matrix_B(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in)  :: filenameA, filenameB

    character(len=str_len)    :: filenameA_out, filenameB_out

    integer                   :: i, j
    character(30)             :: char_A_r, char_A_i, char_B
    character(8)              :: char_idx1, char_idx2

    filenameA_out = trim("output/" // trim(filenameA) // ".txt")
    filenameB_out = trim("output/" // trim(filenameB) // ".txt")

    call open_file(mat_A_out, filenameA_out, .false.)
    call open_file(mat_B_out, filenameB_out, .false.)

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        write(char_idx1, form_int) i
        write(char_idx2, form_int) j

        write(char_A_r, form_eout) real(matrix_A(i, j))
        write(char_A_i, form_eout) aimag(matrix_A(i, j))
        write(char_B,   form_eout) matrix_B(i, j)

        write(mat_A_out, *) char_idx1, ",", &
                            char_idx2, ",", &
                            char_A_r, ",", char_A_i
        write(mat_B_out, *) char_idx1, ",", &
                            char_idx2, ",", &
                            char_B
      end do
    end do

    close(mat_A_out)
    close(mat_B_out)
  end subroutine save_matrices


  subroutine save_eigenvectors(vl, vr, filenameL, filenameR)
    use mod_eigenfunctions

    complex(dp), intent(in)      :: vl(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)      :: vr(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in) :: filenameL, filenameR

    character(len=str_len)       :: filenameL_out, filenameR_out

    integer                      :: i, j
    character(30)                :: char_l_r, char_l_i, char_r_r, char_r_i
    character(8)                 :: char_idx1, char_idx2

    filenameL_out = trim("output/" // trim(filenameL) // ".txt")
    filenameR_out = trim("output/" // trim(filenameR) // ".txt")

    call open_file(eigenf_l_out, filenameL_out, .false.)
    call open_file(eigenf_r_out, filenameR_out, .false.)

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        write(char_idx1, form_int) i
        write(char_idx2, form_int) j

        write(char_l_r, form_eout) real(vl(i, j))
        write(char_l_i, form_eout) aimag(vl(i, j))
        write(char_r_r, form_eout) real(vr(i, j))
        write(char_r_i, form_eout) aimag(vr(i, j))

        write(eigenf_l_out, *) char_idx1, ",", &
                               char_idx2, ",", &
                               char_l_r, ",", char_l_i
        write(eigenf_r_out, *) char_idx1, ",", &
                               char_idx2, ",", &
                               char_r_r, ",", char_r_i
      end do
    end do

    close(eigenf_l_out)
    close(eigenf_r_out)
  end subroutine save_eigenvectors



  subroutine plot_results()
    if (plot_when_finished) then
      write(*, *) "Plotting results..."
      call execute_command_line('python python/plot_eigenvalues.py')
      call execute_command_line('python python/plot_matrices.py')
    end if
    return
  end subroutine plot_results

end module mod_io
