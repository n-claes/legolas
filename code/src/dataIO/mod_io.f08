module mod_io
  use mod_global_variables
  implicit none


  public

  integer, parameter  :: w_output = 10

contains

  subroutine save_eigenvalues(omega)
    complex(dp), intent(in) :: omega(matrix_gridpts)
    integer                 :: i

    open (w_output, file='eigenvalues.txt', status='unknown')

    do i = 1, matrix_gridpts
      write(w_output, *) omega(i)
    end do

  end subroutine save_eigenvalues

end module mod_io
