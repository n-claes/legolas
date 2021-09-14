!======================================
!> Module to provide timing facilities.
!!
!!### Usage
!! The simplest way is using the module variable.
!!
!! Using wall time
!!```fortran
!!call tic()
!!call slow_operation()
!!call toc("slow operation is done")
!!```
!! or CPU time
!!```fortran
!!call cputic()
!!call slow_operation()
!!call cputoc("slow operation is done")
!!```
!!
!! Alternatively, if thread safety is required,
!! you can keep track of the start time yourself.
!!```fortran
!!integer :: start_time
!!call tic(start_time)
!!call slow_operation()
!!call toc("slow operation is done", start_time)
!!```
!! or
!!```fortran
!!real :: start_time
!!call cputic(start_time)
!!call slow_operation()
!!call cputoc("slow operation is done", start_time)
!!```
module mod_timing
  use mod_global_variables, only: dp
  use mod_logging, only: log_message, str
  implicit none

  private
  public tic, toc, cputic, cputoc

  !> Wall clock time at last call of `tic()`.
  !! Intially -1000 to more easily detect incorrect use.
  integer, save :: wall_clock_tic = -1000
  !> CPU time at last call of `cputic()`.
  !! Intially -1000 to more easily detect incorrect use.
  real,    save :: cpu_clock_tic = -1000.0

contains

  !> Subroutine to start a wall clock timer.
  !!
  !! `tic()` starts the wall clock timer. The subroutines records the wall clock
  !! time at execution of this command. Use `toc(message)` to log the elapsed time
  !! since the last call of `tic()`.
  !!
  !! @warning As `tic()` and `toc(message)` rely on a module variable, they are
  !! not thread safe. @endwarning
  !!
  !! `tic(start_time)` writes the value of the wall clock time at execution
  !! of this command to **start_time**. Calling this command, shall not
  !! change the module variable associated with calling `tic()`.
  subroutine tic(start_time)
    !> Optional output for the start time. If not present the time is
    !! written to a module variable.
    integer, optional, intent(out) :: start_time

    if (present(start_time)) then
        call system_clock(start_time)
    else
        call system_clock(wall_clock_tic)
    end if
  end subroutine tic

  !> Subroutine to end a wall clock timer.
  !!
  !! `toc(message)` logs the elapsed wall clock time in seconds since the
  !! most recent call to `tic()` along with **message** as a debug message.
  !!
  !! @warning As `tic()` and `toc(message)` rely on a module variable, they are
  !! not thread safe. @endwarning
  !!
  !! @warning If `toc(message)` is called without first calling `tic()` the result
  !! will be meaningless. @endwarning
  !!
  !! `toc(message, start_time)` logs the elapsed wall clock time in seconds since
  !! the call of `tic(start_time)`, along with **message**.
  subroutine toc(message, start_time)
    !> Message to log along elapsed time.
    character(*),      intent(in) :: message
    !> Optional starting time. If not present the time recorded in
    !! the module variable is used.
    integer, optional, intent(in) :: start_time

    integer  :: selected_start_time, rate, end_time
    real(dp) :: elapsed_time

    if (present(start_time)) then
        selected_start_time = start_time
    else
        selected_start_time = wall_clock_tic
    end if

    call system_clock(end_time, rate)

    elapsed_time = real(end_time - selected_start_time, kind=dp) / rate
    call log_message(message // " (" // str(elapsed_time, fmt="f20.3") // " s)", &
                   & level="debug")
  end subroutine toc

  !> Subroutine to start a CPU timer.
  !!
  !! `cputic()` starts the CPU timer. The subroutines records the processor time
  !! at execution of this command. Use `cputoc(message)` to log the elapsed time
  !! since the last call of `cputic()`.
  !!
  !! @warning As `cputic()` and `cputoc(message)` rely on a module variable, they
  !! are not thread safe. @endwarning
  !!
  !! `cputic(start_time)` writes the value of the processor time at execution
  !! of this command to **start_time**. Calling this command, shall not change
  !! the module variable associated with calling `cputic()`.
  subroutine cputic(start_time)
    !> Optional output for the start time. If not present the time is
    !! written to a module variable.
    real, optional, intent(out) :: start_time

    if (present(start_time)) then
        call cpu_time(start_time)
    else
        call cpu_time(cpu_clock_tic)
    end if
  end subroutine cputic

  !> Subroutine to end a CPU timer.
  !!
  !! `cputoc(message)` logs the elapsed CPU time in seconds since the most recent
  !! call to `cputic()` along with **message** as a debug message.
  !!
  !! @warning As `cputic()` and `cputoc(message)` rely on a module variable, they
  !! are not thread safe. @endwarning
  !!
  !! @warning If `cputoc(message)` is called without first calling `cputic()` the
  !! result will be meaningless. @endwarning
  !!
  !! `cputoc(message, start_time)` logs the elapsed CPU time in seconds since the
  !! call of `cputic(start_time)`, along with **message**.
  subroutine cputoc(message, start_time)
    !> Message to log along elapsed time.
    character(*),   intent(in) :: message
    !> Optional starting time. If not present the time recorded in
    !! the module variable is used.
    real, optional, intent(in) :: start_time

    real     :: selected_start_time, end_time
    real(dp) :: elapsed_time

    if (present(start_time)) then
        selected_start_time = start_time
    else
        selected_start_time = cpu_clock_tic
    end if

    call cpu_time(end_time)

    elapsed_time = end_time - selected_start_time
    call log_message(message // " (" // str(elapsed_time, fmt="f20.3") // " s, CPU time)", &
                   & level="debug")
  end subroutine cputoc

end module mod_timing
