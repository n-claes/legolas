module mod_logging
  use mod_global_variables, only: logging_level
  implicit none

  !! Format settings
  !> exponential format
  character(8), parameter    :: exp_fmt = '(e12.8)'
  !> shorter float format
  character(8), parameter    :: dp_fmt = '(f12.8)'
  !> integer format
  character(4), parameter    :: int_fmt  = '(i8)'
  !> character used as variable to log non-strings
  character(12) :: char_log

  ! logging level:
  ! 0: nothing logged, only errors are thrown
  ! 1: errors + warnings are logged
  ! 2: info + errors + warnings are logged
  ! 3: everything is logged
  ! TODO: maybe create a .log file which saves everything

  private

  public :: log_message
  public :: print_console_info
  public :: print_whitespace
  public :: char_log, exp_fmt, dp_fmt, int_fmt

contains

  subroutine log_message(msg, level)
    character(len=*), intent(in)  :: msg, level

    select case(level)
    case('error')
      write(*, *) "[   ERROR   ] ", msg
      error stop
    case('warning')
      if (logging_level >= 1) then
        write(*, *) "[  WARNING  ] ", msg
      end if
    case('info')
      if (logging_level >= 2) then
        write(*, *) "[   INFO    ] ", msg
      end if
    case('debug')
      if (logging_level >=3) then
        write(*, *) "[   DEBUG   ] ", msg
      end if
    case default
      write(*, *) "level argument should be 'error', 'warning', 'info' or 'debug'."
    end select
  end subroutine log_message

  subroutine print_console_info()
    use mod_global_variables
    use mod_equilibrium_params, only: k2, k3

    if (logging_level <= 1) then
      return
    end if

    call print_whitespace(1)
    write(*, *) " _        _______  _______  _______  _        _______  _______ "
    write(*, *) "( \      (  ____ \(  ____ \(  ___  )( \      (  ___  )(  ____ \"
    write(*, *) "| (      | (    \/| (    \/| (   ) || (      | (   ) || (    \/"
    write(*, *) "| |      | (__    | |      | |   | || |      | (___) || (_____ "
    write(*, *) "| |      |  __)   | | ____ | |   | || |      |  ___  |(_____  )"
    write(*, *) "| |      | (      | | \_  )| |   | || |      | (   ) |      ) |"
    write(*, *) "| (____/\| (____/\| (___) || (___) || (____/\| )   ( |/\____) |"
    write(*, *) "(_______/(_______/(_______)(_______)(_______/|/     \|\_______)"
    call print_whitespace(2)

    write(*, *) "Running with the following configuration:"
    call print_whitespace(1)

    ! Geometry info
    write(*, *) "-- Geometry settings --"
    write(*, *) "Geometry           : ", geometry
    write(char_log, dp_fmt) x_start
    write(*, *) "Grid start         : ", adjustl(char_log)
    write(char_log, dp_fmt) x_end
    write(*, *) "Grid end           : ", adjustl(char_log)
    write(char_log, int_fmt) gridpts
    write(*, *) "Gridpoints         : ", adjustl(char_log)
    write(char_log, int_fmt) gauss_gridpts
    write(*, *) "Gaussian gridpoints: ", adjustl(char_log)
    write(char_log, int_fmt) matrix_gridpts
    write(*, *) "Matrix gridpoints  : ", adjustl(char_log)
    call print_whitespace(1)

    ! Equilibrium info
    write(*, *) "-- Equilibrium settings --"
    write(*, *) "Equilibrium type   : ", equilibrium_type
    write(*, *) "Boundary conditions: ", boundary_type
    write(char_log, dp_fmt) gamma
    write(*, *) "Gamma              : ", adjustl(char_log)
    write(char_log, dp_fmt) k2
    write(*, *) "Wave number k2     : ", adjustl(char_log)
    write(char_log, dp_fmt) k3
    write(*, *) "Wave number k3     : ", adjustl(char_log)
    call print_whitespace(1)

    ! Save info
    write(*, *) "-- DataIO settings --"
    call logical_tostring(write_matrices, char_log)
    write(*, *) "Write matrices to file       : ", char_log
    call logical_tostring(write_eigenfunctions, char_log)
    write(*, *) "Write eigenfunctions to file : ", char_log
    call print_whitespace(2)

  end subroutine print_console_info

  !> Converts a Fortran logical ('T' or 'F') to a string ('true', 'false').
  !! @param[in] boolean   Fortran logical to convert
  !! @param[out] boolean_string   'true' if boolean == True, 'false' otherwise
  subroutine logical_tostring(boolean, boolean_string)
    logical, intent(in)             :: boolean
    character(len=12), intent(out)  :: boolean_string

    if (boolean) then
      boolean_string = 'True'
    else
      boolean_string = 'False'
    end if
  end subroutine logical_tostring

  subroutine print_whitespace(lines)
    integer, intent(in) :: lines
    integer :: i

    if (logging_level >= 1) then
      do i = 1, lines
        write(*, *) ""
      end do
    end if
  end subroutine print_whitespace

end module mod_logging
