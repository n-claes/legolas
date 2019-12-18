module mod_output
  use mod_global_variables
  implicit none

  public

  !! IO units -- do not use 0/5/6/7, these are system-reserved
  !> IO unit for eigenvalues array omega
  integer, parameter  :: omega_unit = 10
  !> IO unit for configuration file
  integer, parameter  :: config_unit = 20
  !> IO unit for matrix A
  integer, parameter  :: mat_a_unit = 30
  !> IO unit for matrix B
  integer, parameter  :: mat_b_unit = 40
  !> IO (base) unit for eigenfunctions (incremented by eigenfunction index 1-8)
  integer, parameter  :: base_ef_unit = 50
  !> IO unit for left eigenvectors
  integer, parameter  :: eigenvector_l_unit = 70
  !> IO unit for right eigenvectors
  integer, parameter  :: eigenvector_r_unit = 80
  !> IO unit for eigenfunction grid (2*gridpts - 1 = ef_gridpts)
  integer, parameter  :: ef_grid_unit = 90

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

contains

  subroutine open_file(file_unit, filename)
    integer, intent(in)           :: file_unit
    character(len=*), intent(in)  :: filename

    open(unit=file_unit, file=filename, access='stream', &
         status='unknown', action='write')
  end subroutine open_file


  subroutine make_filename(base_filename, filename)
    character(len=*), intent(in)  :: base_filename
    character(len=*), intent(out) :: filename

    filename = trim(trim(output_folder) // trim(base_filename) // &
                    trim(file_extension))
  end subroutine make_filename


  subroutine omegas_tofile(omega, base_filename)
    complex(dp), intent(in)       :: omega(matrix_gridpts)
    character(len=*), intent(in)  :: base_filename

    character(str_len)            :: filename

    call make_filename(base_filename, filename)

    call open_file(omega_unit, filename)
    write(unit=omega_unit) omega

    close(unit=omega_unit)
  end subroutine omegas_tofile


  subroutine matrices_tofile(matrix_A, matrix_B, base_filenameA, base_filenameB)
    complex(dp), intent(in)       :: matrix_A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(in)          :: matrix_B(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in)  :: base_filenameA, base_filenameB

    character(str_len)            :: filenameA, filenameB
    integer                       :: i

    call make_filename(base_filenameA, filenameA)
    call make_filename(base_filenameB, filenameB)

    call open_file(mat_a_unit, filenameA)
    call open_file(mat_b_unit, filenameB)

    do i = 1, matrix_gridpts
      write(mat_a_unit) matrix_A(i, :)
      write(mat_b_unit) matrix_B(i, :)
    end do

    close(mat_a_unit)
    close(mat_b_unit)
  end subroutine matrices_tofile


  subroutine eigenvectors_tofile(ev_l, ev_r, base_filenameL, base_filenameR)
    complex(dp), intent(in)       :: ev_l(matrix_gridpts, matrix_gridpts)
    complex(dp), intent(in)       :: ev_r(matrix_gridpts, matrix_gridpts)
    character(len=*), intent(in)  :: base_filenameL, base_filenameR

    character(str_len)            :: filenameL, filenameR
    integer                       :: i

    call make_filename(base_filenameL, filenameL)
    call make_filename(base_filenameR, filenameR)

    call open_file(eigenvector_l_unit, filenameL)
    call open_file(eigenvector_r_unit, filenameR)

    do i = 1, matrix_gridpts
      write(eigenvector_l_unit) ev_l(i, :)
      write(eigenvector_r_unit) ev_r(i, :)
    end do

    close(eigenvector_l_unit)
    close(eigenvector_r_unit)
  end subroutine eigenvectors_tofile


  subroutine ef_grid_tofile(ef_grid, base_filename)
    real(dp), intent(in)            :: ef_grid(ef_gridpts)
    character(len=*), intent(in)    :: base_filename

    character(str_len)              :: filename

    call make_filename(base_filename, filename)

    call open_file(ef_grid_unit, filename)

    write(ef_grid_unit) ef_grid
    close(ef_grid_unit)
  end subroutine ef_grid_tofile


  subroutine eigenfunction_tofile(single_ef)
    use mod_types

    type(ef_type), intent(in)       :: single_ef

    character                       :: idx
    character(str_len)              :: base_filename, filename
    integer                         :: single_ef_unit, w_idx

    ! current eigenfunction index (1 = rho, 2 = v1 etc.)
    write(idx, '(i0)') single_ef % index
    ! get output unit for this eigenfunction, increment base unit with index
    single_ef_unit = base_ef_unit + single_ef % index

    ! extend base filename with output folder and index
    base_filename = 'eigenfunctions/' // trim(idx) // '_' &
                    // trim(var_eigenf % savename)
    call make_filename(base_filename_ext, filename)

    call open_file(single_ef_unit, filename)

    do w_idx = 1, matrix_gridpts
      write(single_ef_unit) single_ef % eigenfunctions(:, w_idx)
    end do
    close(single_ef_unit)
  end subroutine eigenfunction_tofile


end module mod_output
