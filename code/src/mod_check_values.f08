module mod_check_values
  use mod_global_variables
  implicit none

contains

  subroutine check_small_values(array)
    complex(dp), intent(inout)  :: array(matrix_gridpts)

    real(dp)                    :: w_real, w_imag
    integer                     :: i

    do i = 1, matrix_gridpts
      w_real = real(array(i))
      w_imag = aimag(array(i))

      if (abs(w_real) < dp_LIMIT) then
        w_real = 0.0d0
      end if

      if (abs(w_imag) < dp_LIMIT) then
        w_imag = 0.0d0
      end if

      array(i) = cmplx(w_real, w_imag, kind=dp)
    end do
  end subroutine check_small_values


  subroutine check_small_values_matrix(matrix)
    complex(dp), intent(inout)  :: matrix(matrix_gridpts, matrix_gridpts)
    real(dp)                    :: v_real, v_imag

    integer                     :: i, j

    do i = 1, matrix_gridpts
      do j = 1, matrix_gridpts
        v_real = real(matrix(i, j))
        v_imag = aimag(matrix(i, j))

        if (abs(v_real) < dp_LIMIT) then
          v_real = 0.0d0
        end if

        if (abs(v_imag) < dp_LIMIT) then
          v_imag = 0.0d0
        end if

        matrix(i, j) = cmplx(v_real, v_imag, kind=dp)
      end do
    end do
  end subroutine check_small_values_matrix


  subroutine remove_large_eigenvalues(omega)
    complex(dp), intent(inout)  :: omega(matrix_gridpts)

    !! \TODO: when changing to ZBIG in matrix A, remove these
    !!        eigenvalues after?
    return
  end subroutine remove_large_eigenvalues

end module mod_check_values
