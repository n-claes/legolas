module mod_info
  use mod_global_variables
  implicit none

  private

  character(20)   :: char

  public :: show_startup_info
  public :: show_final_info


contains

  subroutine show_startup_info()


    write(*, *) "------------------------------"
    write(*, *) "----------- LEGOLAS ----------"
    write(*, *) "------------------------------"
    write(*, *) ""

    write(*, *) "Running with the following configuration:"
    write(*, *)

    ! Geometry info
    write(*, *) "-- Geometry settings --"
    write(*, *) "Coordinate system: ", geometry
    write(char, form_fout) x_start
    write(*, *) "Start              : ", adjustl(char)
    write(char, form_fout) x_end
    write(*, *) "End                : ", adjustl(char)
    write(char, form_int) gridpts
    write(*, *) "Gridpoints         : ", adjustl(char)
    write(char, form_int) matrix_gridpts
    write(*, *) "Matrix gridpoints  : ", adjustl(char)
    write(char, form_int) gauss_gridpts
    write(*, *) "Gaussian gridpoints: ", adjustl(char)
    write(*, *) ""

    ! Equilibrium info
    write(*, *) "-- Equilibrium settings --"
    write(*, *) "Equilibrium type   : ", equilibrium_type
    write(*, *) "Boundary conditions: ", boundary_type
    write(char, form_fout) k2
    write(*, *) "Wave number k2     : ", adjustl(char)
    write(char, form_fout) k3
    write(*, *) "Wave number k3     : ", adjustl(char)
    write(*, *) ""

    ! Physics info
    write(*, *) "-- Physics settings --"
    write(char, form_fout) gamma
    write(*, *) "Gamma              : ", adjustl(char)
    call get_bool(flow, char)
    write(*, *) "Flow               : ", char
    call get_bool(radiative_cooling, char)
    write(*, *) "Radiative coling   : ", char
    if (radiative_cooling) then
      write(*, *) "  Cooling curve  : ", cooling_curve
    end if
    call get_bool(thermal_conduction, char)
    write(*, *) "Thermal conduction : ", char
    call get_bool(resistivity, char)
    write(*, *) "Resistivity        : ", char
    if (resistivity) then
      call get_bool(use_fixed_resistivity, char)
      write(*, *) "  Fixed resistivity: ", char
      if (use_fixed_resistivity) then
        write(char, form_fout) fixed_eta_value
        write(*, *) "  Value for eta    : ", adjustl(char)
      end if
    end if
    call get_bool(external_gravity, char)
    write(*, *) "External gravity   : ", char
    call get_bool(cgs_units, char)
    write(*, *) "Using CGS units    : ", char
    write(*, *) ""

    ! Save info
    write(*, *) "-- DataIO settings --"
    call get_bool(write_matrices, char)
    write(*, *) "Write matrices to file       : ", char
    call get_bool(write_eigenvectors, char)
    write(*, *) "Write eigenvectors to file   : ", char
    call get_bool(write_eigenfunctions, char)
    write(*, *) "Write eigenfunctions to file : ", char
    call get_bool(show_results, char)
    write(*, *) "Showing results              : ", char
    call get_bool(show_matrices, char)
    write(*, *) "Showing matrices             : ", char
    call get_bool(show_eigenfunctions, char)
    write(*, *) "Showing eigenfunctions       : ", char

    call print_line()

  end subroutine show_startup_info


  subroutine show_final_info()
    return
    !use mod_timer
    !real(dp)     :: total_runtime, total_matrixtime, total_writetime
    !character(6) :: form_fout_short = '(f8.1)'
    !
    !write(*,*) "Finished Legolas"
    !
    ! call get_total_runtime(total_runtime)
    !
    ! if (total_runtime < 10.0d0) then
    !   write(*, *) ""
    !   write(*, *) "Finished ESoNaS."
    !   return
    ! end if
    ! total_matrixtime = time_invert + time_multiply + time_QR
    ! total_writetime  = time_write_omegas + time_write_matrices + &
    !                    time_write_eigenvectors + time_write_eigenfunctions
    !
    ! call print_line()
    ! write(*, *) "Division of total runtime:"
    ! write(*, *) ""
    !
    ! write(char, form_fout_short) time_initialisation
    ! write(*, *) "Initialisation             : ", trim(adjustl(char)), " s"
    !
    ! write(char, form_fout_short) time_matrices
    ! write(*, *) "Matrix Creation            : ", trim(adjustl(char)), " s"
    !
    ! write(char, form_fout_short) total_matrixtime
    ! write(*, *) "Solvers                    : ", trim(adjustl(char)), " s"
    ! write(*, *) "    of which"
    ! write(char, form_fout_short) time_invert
    ! write(*, *) "  - Matrix B inversion     : ", trim(adjustl(char)), " s"
    ! write(char, form_fout_short) time_multiply
    ! write(*, *) "  - Matrix multiplication  : ", trim(adjustl(char)), " s"
    ! write(char, form_fout_short) time_QR
    ! write(*, *) "  - QR algorithm           : ", trim(adjustl(char)), " s"
    !
    ! write(char, form_fout_short) total_writetime
    ! write(*, *) "Data IO                    : ", trim(adjustl(char)), " s"
    ! write(*, *) "    of which"
    ! write(char, form_fout_short) time_write_omegas
    ! write(*, *) "  - Writing omegas         : ", trim(adjustl(char)), " s"
    ! write(char, form_fout_short) time_write_matrices
    ! write(*, *) "  - Writing matrices       : ", trim(adjustl(char)), " s"
    ! write(char, form_fout_short) time_write_eigenvectors
    ! write(*, *) "  - Writing eigenvectors   : ", trim(adjustl(char)), " s"
    ! write(char, form_fout_short) time_write_eigenfunctions
    ! write(*, *) "  - Writing eigenfunctions : ", trim(adjustl(char)), " s"
    !
    ! write(*, *) "----------------------------------------------------"
    ! write(char, form_fout_short) total_runtime
    ! write(*, *) "Finished ESoNaS in ", trim(adjustl(char)), " seconds"
    ! write(*, *) ""
  end subroutine show_final_info



  subroutine get_bool(bool, bool_string)
    logical, intent(in)            :: bool
    character(len=20), intent(out) :: bool_string

    if (bool) then
      bool_string = "true"
    else
      bool_string = "false"
    end if
  end subroutine get_bool


  subroutine print_line()
    write(*, *) ""
    write(*, *) "----------------------------------------------------"
    write(*, *) ""
  end subroutine print_line




end module mod_info
