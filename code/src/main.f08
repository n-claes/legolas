program esonas
  implicit none

  character(:), allocatable     :: geometry
  integer                       :: gridpts, matrix_gridpts
  complex, allocatable          :: matrix_A(:, :)
  real, allocatable             :: matrix_B(:, :)


  geometry   = "cylindrical"
  gridpts = 1
  matrix_gridpts = 16 * gridpts

  allocate(matrix_A(matrix_gridpts, matrix_gridpts))
  allocate(matrix_B(matrix_gridpts, matrix_gridpts))

  call setup_B(matrix_gridpts, matrix_B)

  ! First find a way to set global parameters?
  ! similar to use mod_global_parameters in amrvac



contains

  subroutine setup_B(matrix_gridpts, matrix_B)
    use setup_matrix_b

    integer, intent(in) ::  matrix_gridpts
    real, intent(inout) ::  matrix_B(matrix_gridpts, matrix_gridpts)

    call construct_B(matrix_gridpts, matrix_B)


    return

  end subroutine setup_B


end program esonas
