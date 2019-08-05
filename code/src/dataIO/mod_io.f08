module mod_io
  use mod_global_variables
  use mod_physical_constants
  implicit none


  public

  integer, parameter  :: w_output            = 10
  integer, parameter  :: config              = 30
  integer, parameter  :: mat_a_out           = 40
  integer, parameter  :: mat_b_out           = 50
  integer, parameter  :: eigenvect_l_out     = 60
  integer, parameter  :: eigenvect_r_out     = 70
  integer, parameter  :: eigenvectors_out    = 80
  integer, parameter  :: grid_out            = 90

contains

  subroutine open_file(file_no, filename, append, stream)
    integer, intent(in)          :: file_no
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: append
    logical, intent(in)          :: stream

    logical                      :: file_present

    if (append) then
      inquire(file=filename, exist=file_present)

      if (file_present) then
        if (stream) then
          open(file_no, file=filename, status='old', position='append', &
               access='stream', action='write')
        else
          open(file_no, file=filename, status='old', position='append', &
               action='write')
        end if
      else
        if (stream) then
          open(file_no, file=filename, status='new', access='stream', &
               action='write')
        else
          open(file_no, file=filename, status='new', action='write')
        end if
      end if
    else
      if (stream) then
        open(file_no, file=filename, access='stream', status='unknown')
      else
        open(file_no, file=filename, status='unknown')
      end if
    end if
  end subroutine open_file



  subroutine save_eigenvalues(omega, filenameW, append, stream)
    complex(dp), intent(in)      :: omega(matrix_gridpts)
    character(len=*), intent(in) :: filenameW
    logical, intent(in)          :: append, stream
    integer                      :: i

    character(30)                :: w_char_r, w_char_i
    character(str_len)           :: filenameW_out

    filenameW_out = trim("output/" // trim(filenameW) // ".dat")

    call open_file(w_output, filenameW_out, append, stream)

    if (stream) then
      write(w_output) omega
    else
      do i = 1, matrix_gridpts
        write(w_char_r, form_e) real(omega(i))
        write(w_char_i, form_e) aimag(omega(i))

        write(w_output, *) w_char_r, ',', w_char_i
      end do
    end if

    close(w_output)

  end subroutine save_eigenvalues



  subroutine save_eigenfunctions(var_eigenf)
    use mod_types

    type(eigenf_type), intent(in) :: var_eigenf

    character(str_len)            :: filenameEF_out
    character                     :: idx
    integer                       :: file_out, w_idx

    write(idx, '(i0)') var_eigenf % index

    filenameEF_out = trim("output/eigenfunctions/" // trim(idx) // "_" &
                          // trim(var_eigenf % savename) // ".dat")
    file_out = var_eigenf % write_out

    call open_file(file_out, filenameEF_out, append=.false., stream=.true.)

    do w_idx = 1, matrix_gridpts
      write(file_out) var_eigenf % eigenfunctions(:, w_idx)
    end do
    close(file_out)
  end subroutine save_eigenfunctions


  subroutine save_eigenf_grid(eigenf_grid)
    real(dp), intent(in)  :: eigenf_grid(eigenf_gridpts)

    call open_file(grid_out, 'output/grid.dat', append=.false., stream=.true.)

    write(grid_out) eigenf_grid
    close(grid_out)
  end subroutine save_eigenf_grid



  subroutine save_config(filenameCFG)
    character(len=*), intent(in)  :: filenameCFG

    character(20)      :: char
    character(str_len) :: filenameCFG_out

    filenameCFG_out = trim("output/" // trim(filenameCFG) // ".txt")

    call open_file(config, filenameCFG_out, .false., .false.)

    write(config, *) "Equilibrium type   : ", equilibrium_type

    write(config, *) "Flow               : ", flow
    write(config, *) "Radiative cooling  : ", radiative_cooling
    if (radiative_cooling) then
      write(config, *) "Cooling curve      : ", cooling_curve
    end if
    write(config, *) "Thermal conduction : ", thermal_conduction
    write(config, *) "Resistivity        : ", resistivity
    if (resistivity) then
      write(config, *) "Fixed Resistivity  : ", use_fixed_resistivity
      if (use_fixed_resistivity) then
        write(char, form_fout) fixed_eta_value
        write(config, *) "Fixed eta value    : ", adjustl(char)
      end if
    end if
    write(config, *) "External gravity   : ", external_gravity
    if (external_gravity) then
      write(config, *) "Custom gravity     : ", use_custom_gravity
      if (use_custom_gravity) then
        write(char, form_fout) custom_g_value
        write(config, *) "g                  : ", adjustl(char)
      else
        write(config, *) "Gravity strength   : ", gravity_type
      end if
    end if

    write(char, form_fout) gamma
    write(config, *) "Gamma              : ", adjustl(char)
    write(char, form_fout) k2
    write(config, *) "k2                 : ", adjustl(char)
    write(char, form_fout) k3
    write(config, *) "k3                 : ", adjustl(char)


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
    write(char, form_int) eigenf_gridpts
    write(config, *) "Eigenfunction gridpoints: ", adjustl(char)
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

    write(config, *) ""

    write(config, *) "Write matrices      : ", write_AB
    write(config, *) "Write eigenvectors  : ", write_eigenvectors
    write(config, *) "Write eigenfunctions: ", write_eigenfunctions
    write(config, *) "Plot when finished  : ", plot_when_finished
    write(config, *) "Plot matrices       : ", plot_AB
    write(config, *) "Plot eigenfunctions : ", plot_eigenfunctions

    close(config)
  end subroutine save_config

  subroutine save_matrices(matrix_A, matrix_B, filenameA, filenameB)
    complex(dp), intent(in)   :: matrix_A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)      :: matrix_B(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in)  :: filenameA, filenameB

    character(len=str_len)    :: filenameA_out, filenameB_out

    integer                   :: i

    filenameA_out = trim("output/" // trim(filenameA) // ".dat")
    filenameB_out = trim("output/" // trim(filenameB) // ".dat")

    call open_file(mat_A_out, filenameA_out, append=.false., stream=.true.)
    call open_file(mat_B_out, filenameB_out, append=.false., stream=.true.)

    do i = 1, matrix_gridpts
      write(mat_B_out) matrix_B(i, :)
      write(mat_A_out) matrix_A(i, :)
    end do

    close(mat_A_out)
    close(mat_B_out)
  end subroutine save_matrices


  subroutine save_eigenvectors(vl, vr, filenameL, filenameR)
    complex(dp), intent(in)      :: vl(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)      :: vr(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in) :: filenameL, filenameR

    character(len=str_len)       :: filenameL_out, filenameR_out

    integer                      :: i

    filenameL_out = trim("output/" // trim(filenameL) // ".dat")
    filenameR_out = trim("output/" // trim(filenameR) // ".dat")

    call open_file(eigenvect_l_out, filenameL_out, append=.false., stream=.true.)
    call open_file(eigenvect_r_out, filenameR_out, append=.false., stream=.true.)

    do i = 1, matrix_gridpts
      write(eigenvect_l_out) vl(i, :)
      write(eigenvect_r_out) vr(i, :)
    end do

    close(eigenvect_l_out)
    close(eigenvect_r_out)
  end subroutine save_eigenvectors

end module mod_io
