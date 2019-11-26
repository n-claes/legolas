module mod_timer
  use mod_global_variables, only: dp
  implicit none

  real(dp), protected :: time_initialisation = 0.0d0
  real(dp), protected :: time_matrices = 0.0d0
  real(dp), protected :: time_invert = 0.0d0
  real(dp), protected :: time_multiply = 0.0d0
  real(dp), protected :: time_QR = 0.0d0
  real(dp), protected :: time_write_omegas = 0.0d0
  real(dp), protected :: time_write_matrices = 0.0d0
  real(dp), protected :: time_write_eigenvectors = 0.0d0
  real(dp), protected :: time_write_eigenfunctions = 0.0d0

  real(dp)            :: start_time
  real(dp)            :: end_time

contains

  subroutine set_time_start()
    call cpu_time(start_time)
  end subroutine set_time_start

  subroutine set_time_initialisation()
    call cpu_time(end_time)
    time_initialisation = end_time - start_time
  end subroutine set_time_initialisation

  subroutine set_time_matrices()
    call cpu_time(end_time)
    time_matrices = end_time - start_time
  end subroutine set_time_matrices

  subroutine set_time_invert()
    call cpu_time(end_time)
    time_invert = end_time - start_time
  end subroutine set_time_invert

  subroutine set_time_multiply()
    call cpu_time(end_time)
    time_multiply = end_time - start_time
  end subroutine set_time_multiply

  subroutine set_time_QR()
    call cpu_time(end_time)
    time_QR = end_time - start_time
  end subroutine set_time_QR

  subroutine set_time_write_omegas()
    call cpu_time(end_time)
    time_write_omegas = end_time - start_time
  end subroutine set_time_write_omegas

  subroutine set_time_write_matrices()
    call cpu_time(end_time)
    time_write_matrices = end_time - start_time
  end subroutine set_time_write_matrices

  subroutine set_time_write_eigenvectors()
    call cpu_time(end_time)
    time_write_eigenvectors = end_time - start_time
  end subroutine set_time_write_eigenvectors

  subroutine set_time_write_eigenfunctions()
    call cpu_time(end_time)
    time_write_eigenfunctions = end_time - start_time
  end subroutine set_time_write_eigenfunctions

  subroutine get_total_runtime(total_time)
    real(dp), intent(out)   :: total_time

    total_time = time_initialisation + time_matrices + time_invert + &
                 time_multiply + time_QR + time_write_omegas + &
                 time_write_matrices + time_write_eigenvectors + &
                 time_write_eigenfunctions
  end subroutine get_total_runtime


end module mod_timer
