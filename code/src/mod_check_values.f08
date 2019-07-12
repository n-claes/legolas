module mod_check_values
  use mod_global_variables
  implicit none

contains

  subroutine check_small_omega(omega)
    complex(dp), intent(inout)  :: omega(matrix_gridpts)

    real(dp)                    :: w_real, w_imag
    integer                     :: i

    do i = 1, matrix_gridpts
      w_real = real(omega(i))
      w_imag = aimag(omega(i))

      if (abs(w_real) < dp_LIMIT) then
        w_real = 0.0d0
      end if

      if (abs(w_imag) < dp_LIMIT) then
        w_imag = 0.0d0
      end if

      omega(i) = cmplx(w_real, w_imag, kind=dp)
    end do
  end subroutine check_small_omega


  subroutine remove_large_eigenvalues(omega)
    complex(dp), intent(inout)  :: omega(matrix_gridpts)

    !! \TODO: when changing to ZBIG in matrix A, remove these
    !!        eigenvalues after?
    return
  end subroutine remove_large_eigenvalues

end module mod_check_values
