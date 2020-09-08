! =============================================================================
!> Main handler for console print statements. The level of information
!! printed to the console depends on the corresponding global variable
!! called <tt>logging_level</tt> defined in the parfile.
!! @note Values for <tt>logging_level</tt> can be set to
!!
!! - If <tt>logging_level = 0</tt>: only critical errors are printed, everything else is suppressed.
!! - If <tt>logging_level = 1</tt>: only errors and warnings are printed.
!! - If <tt>logging_level = 2</tt>: errors, warnings and info messages are printed. This
!!                                  is the default value
!! - If <tt>logging_level = 3+</tt>: prints all of the above, including debug messages. @endnote
module mod_logging
  use mod_global_variables, only: logging_level
  implicit none

  !> exponential format
  character(8), parameter    :: exp_fmt = '(e20.8)'
  !> shorter float format
  character(8), parameter    :: dp_fmt = '(f20.8)'
  !> integer format
  character(4), parameter    :: int_fmt  = '(i8)'
  !> character used as variable to log non-strings
  character(20) :: char_log

  private

  public :: log_message
  public :: print_logo
  public :: print_console_info
  public :: print_whitespace
  public :: char_log, exp_fmt, dp_fmt, int_fmt

contains


  !> Logs messages to the console. Every message will be prepended by
  !! [  LEVEL  ] to indicate its type.
  !! @warning An error is thrown if a wrong level is passed. @endwarning
  !! @note The argument <tt>level</tt> can be 'error', 'warning', 'info' or 'debug'.
  !!       The 'error' level corresponds to throwing a critical error and stops code execution.
  subroutine log_message(msg, level)
    use mod_exceptions, only: raise_exception

    !> the message to print to the console
    character(len=*), intent(in)  :: msg
    !> the level (severity) of the message
    character(len=*), intent(in)  :: level

    select case(level)
    case('error')
      call raise_exception(msg)
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
      write(*, *) "[   ERROR   ] level argument should be 'error', 'warning', 'info' or 'debug'."
      error stop
    end select
  end subroutine log_message


  !> Prints the Legolas logo to the console.
  !! The logo is wrapped in 1 whitespace at the top and
  !! two at the bottom. Only for logging level 'warning' (1) and above
  subroutine print_logo()
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
  end subroutine print_logo


  !> Prints various console messages showing geometry, grid parameters,
  !! equilibrium parameters etc. Only for logging level "info" or above.
  subroutine print_console_info()
    use mod_global_variables
    use mod_equilibrium_params, only: k2, k3

    if (logging_level <= 1) then
      return
    end if

    call log_message("--------------- Configuration ---------------", level='info')

    call log_message("                Grid settings", level='info')
    call log_message("geometry             : " // adjustl(trim(geometry)), level='info')
    write(char_log, dp_fmt) x_start
    call log_message("grid start           : " // adjustl(char_log), level='info')
    write(char_log, dp_fmt) x_end
    call log_message("grid end             : " // adjustl(char_log), level='info')
    write(char_log, int_fmt) gridpts
    call log_message("gridpoints (base)    : " // adjustl(char_log), level='info')
    write(char_log, int_fmt) gauss_gridpts
    call log_message("gridpoints (Gauss)   : " // adjustl(char_log), level='info')
    write(char_log, int_fmt) matrix_gridpts
    call log_message("gridpoints (matrix)  : " // adjustl(char_log), level='info')

    call log_message("             Equilibrium settings", level='info')
    call log_message("selected equilibrium : " // adjustl(trim(equilibrium_type)), level='info')
    call log_message("boundary conditions  : " // adjustl(trim(boundary_type)), level='info')
    write(char_log, dp_fmt) k2
    call log_message("wave number k2       : " // adjustl(char_log), level='info')
    write(char_log, dp_fmt) k3
    call log_message("wave number k3       : " // adjustl(char_log), level='info')

    call log_message("               DataIO settings", level='info')
    call log_message("datfile name         : " // adjustl(trim(basename_datfile)), level='info')
    call log_message("output folder        : " // adjustl(trim(output_folder)), level='info')
    call logical_tostring(write_matrices, char_log)
    call log_message("write matrices       : " // adjustl(char_log), level='info')
    call logical_tostring(write_eigenfunctions, char_log)
    call log_message("write eigenfunctions : " // adjustl(char_log), level='info')
    call log_message("---------------------------------------------", level='info')
    call print_whitespace(1)
  end subroutine print_console_info


  !> Converts a given Fortran logical to a string "true" or "false".
  subroutine logical_tostring(boolean, boolean_string)
    !> logical to convert
    logical, intent(in)             :: boolean
    !> <tt>True</tt> if boolean == True, <tt>False</tt> otherwise
    character(len=20), intent(out)  :: boolean_string

    if (boolean) then
      boolean_string = 'True'
    else
      boolean_string = 'False'
    end if
  end subroutine logical_tostring


  !> Prints an empty line to the console.
  !! Only if logging level is 'warning' or above.
  subroutine print_whitespace(lines)
    !> amount of empty lines to print
    integer, intent(in) :: lines
    integer :: i

    if (logging_level >= 1) then
      do i = 1, lines
        write(*, *) ""
      end do
    end if
  end subroutine print_whitespace

end module mod_logging
