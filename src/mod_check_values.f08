module mod_check_values
  use mod_global_variables
  implicit none

  interface check_small_values
    module procedure small_values_real_array
    module procedure small_values_complex_array
    module procedure small_values_real_matrix
    module procedure small_values_complex_matrix
  end interface check_small_values

contains

  subroutine small_values_real_array(array)
    real(dp), intent(inout) :: array(:)
    integer                 :: i

    do i = 1, size(array)
      if (abs(array(i)) < dp_LIMIT) then
        array(i) = 0.0d0
      end if
    end do
  end subroutine small_values_real_array

  subroutine small_values_complex_array(array)
    complex(dp), intent(inout) :: array(:)
    integer                    :: i
    real(dp)                   :: a_real, a_imag

    do i = 1, size(array)
      a_real = real(array(i))
      a_imag = aimag(array(i))

      if (abs(a_real) < dp_LIMIT) then
        a_real = 0.0d0
      end if
      if (abs(a_imag) < dp_LIMIT) then
        a_imag = 0.0d0
      end if

      array(i) = cmplx(a_real, a_imag, kind=dp)
    end do
  end subroutine small_values_complex_array

  subroutine small_values_real_matrix(matrix)
    real(dp), intent(inout)    :: matrix(:, :)
    integer                    :: i, j

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        if (abs(matrix(i, j)) < dp_LIMIT) then
          matrix(i, j) = 0.0d0
        end if
      end do
    end do
  end subroutine small_values_real_matrix

  subroutine small_values_complex_matrix(matrix)
    complex(dp), intent(inout) :: matrix(:, :)
    integer                    :: i, j
    real(dp)                   :: a_real, a_imag

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        a_real = real(matrix(i, j))
        a_imag = aimag(matrix(i, j))

        if (abs(a_real) < dp_LIMIT) then
          a_real = 0.0d0
        end if
        if (abs(a_imag) < dp_LIMIT) then
          a_imag = 0.0d0
        end if

        matrix(i, j) = cmplx(a_real, a_imag, kind=dp)
      end do
    end do
  end subroutine small_values_complex_matrix

  subroutine check_negative_array(array, variable_name)
    real(dp), intent(in)          :: array(:)
    character(len=*), intent(in)  :: variable_name
    integer                       :: i

    do i = 1, size(array)
      if (array(i) < 0.0d0) then
        write(*, *) "WARNING: ", trim(variable_name), " is negative somewhere!"
        stop
      end if
    end do
  end subroutine check_negative_array


end module mod_check_values
