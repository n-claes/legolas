! =============================================================================
!> Main handler for console print statements. The level of information
!! printed to the console depends on the corresponding global variable
!! called <tt>logging_level</tt> defined in the parfile.
!! @note Values for <tt>logging_level</tt> can be set to
!!
!! - If <tt>logging_level = 0</tt>: only critical errors are printed, everything else
!!                                  is suppressed.
!! - If <tt>logging_level = 1</tt>: only errors and warnings are printed.
!! - If <tt>logging_level = 2</tt>: errors, warnings and info messages are printed.
!!                                  This is the default value.
!! - If <tt>logging_level = 3+</tt>: prints all of the above,
!!                                   including debug messages. @endnote
module mod_logging
  use mod_global_variables, only: logging_level, str_len
  use mod_painting, only: paint_string
  implicit none

  !> exponential format
  character(8), parameter    :: exp_fmt = '(e20.8)'
  !> shorter float format
  character(8), parameter    :: dp_fmt = '(f20.8)'
  !> integer format
  character(4), parameter    :: int_fmt  = '(i8)'
  !> a convenient "tostring" interface, used for easy console writing
  interface str
    module procedure logical_tostr
    module procedure int_tostr
    module procedure float_tostr
    module procedure complex_tostr
    module procedure character_arr_tostr
  end interface str

  !> logical used to (locally) force a prefix override
  logical :: override_prefix_to_false = .false.

  private

  public :: log_message
  public :: str
  public :: exp_fmt, dp_fmt, int_fmt
  public :: override_prefix_to_false

contains


  !> Logs messages to the console. Every message will be prepended by
  !! [  LEVEL  ] to indicate its type. If this is not desired, set
  !! <tt>use_prefix = .false.</tt>.
  !! @warning An error is thrown if a wrong level is passed. @endwarning
  !! @note The argument <tt>level</tt> can be 'error', 'warning', 'info' or 'debug'.
  !!       The 'error' level corresponds to throwing a critical error and
  !!       stops code execution.
  !!       Error messages are printed in red, warnings in yellow, info messages have
  !!       default colouring and debug messages are in green.
  subroutine log_message(msg, level, use_prefix) ! LCOV_EXCL_START
    use mod_exceptions, only: raise_exception

    !> the message to print to the console
    character(len=*), intent(in)  :: msg
    !> the level (severity) of the message, default is <tt>"info"</tt>
    character(len=*), intent(in), optional  :: level
    !> prefixes message type to string, default is <tt>.true.</tt>
    logical, intent(in), optional :: use_prefix

    ! need a bit more room here, we trim anyway when printing
    character(len=2*str_len) :: msg_painted
    character(:), allocatable :: loglevel
    logical :: add_prefix

    add_prefix = .true.
    if (present(use_prefix)) then
      add_prefix = use_prefix
    end if
    ! override prefix if desired
    if (override_prefix_to_false) then
      add_prefix = .false.
    end if
    if (present(level)) then
      loglevel = level
    else
      loglevel = "info"
    end if

    select case(loglevel)
    case("error")
      call raise_exception(msg)
    case("warning")
      if (logging_level >= 1) then
        if (add_prefix) then
          call paint_string(" WARNING | " // msg, "yellow", msg_painted)
        else
          call paint_string("           " // msg, "yellow", msg_painted)
        end if
        write(*, *) trim(msg_painted)
      end if
    case("info")
      if (logging_level >= 2) then
        if (add_prefix) then
          write(*, *) " INFO    | " // trim(msg)
        else
          write(*, *) "           " // trim(msg)
        end if
      end if
    case("debug")
      if (logging_level >=3) then
        if (add_prefix) then
          call paint_string(" DEBUG   | " // msg, "green", msg_painted)
        else
          call paint_string("         | " // msg, "green", msg_painted)
        end if
        write(*, *) trim(msg_painted)
      end if
    case default
      call raise_exception( &
        "argument 'level' should be 'error', 'warning', 'info' or 'debug'" &
      )
      error stop
    end select
  end subroutine log_message ! LCOV_EXCL_STOP


  !> Converts a given logical to a string "True" or "False".
  function logical_tostr(boolean) result(result_str)
    !> logical to convert
    logical, intent(in) :: boolean
    !> return string, made allocatable so it has same length as input
    character(:), allocatable :: result_str

    if (boolean) then
      result_str = "True"
    else
      result_str = "False"
    end if
  end function logical_tostr


  !> Converts a given integer to a string, the default format is "i8".
  function int_tostr(value, fmt) result(result_str)
    !> integer to convert
    integer, intent(in) :: value
    !> optional format used for writing, default "i8"
    character(len=*), intent(in), optional  :: fmt
    !> return string, made allocatable so it has same length as input
    character(:), allocatable :: result_str
    character(len=20) :: format, char_log

    if (present(fmt)) then
      format = "(" // trim(fmt) // ")"
    else
      format = int_fmt
    end if
    write(char_log, format) value
    result_str = trim(adjustl(char_log))
  end function int_tostr


  !> Converts a given float to a string, the default format is "f20.8".
  function float_tostr(value, fmt) result(result_str)
    use mod_global_variables, only: dp

    !> float to convert
    real(dp), intent(in)  :: value
    !> optional format use for writing, default "f20.8"
    character(len=*), intent(in), optional  :: fmt
    !> return string, made allocatable so it has same length as input
    character(:), allocatable :: result_str
    character(len=20) :: format, char_log

    if (present(fmt)) then
      format = "(" // trim(fmt) // ")"
    else
      format = dp_fmt
    end if
    write(char_log, format) value
    result_str = trim(adjustl(char_log))
  end function float_tostr


  !> Converts a given complex number to a string, the default format is "f20.8".
  !! This will be printed in the form xxxx + xxxxi.
  function complex_tostr(value, fmt) result(result_str)
    use mod_global_variables, only: dp

    !> complex to convert
    complex(dp), intent(in) :: value
    !> optional format use for writing, default "f20.8"
    character(len=*), intent(in), optional  :: fmt
    !> return string, made allocatable so it has same length as input
    character(:), allocatable :: result_str
    character(len=20) :: format, char_log, char_log2

    if (present(fmt)) then
      format = "(" // trim(fmt) // ")"
    else
      format = "(f18.8)"
    end if
    write(char_log, format) real(value)
    write(char_log2, '(SP,' // format // ',"i")') aimag(value)
    result_str = trim(adjustl(char_log)) // trim(adjustl(char_log2))
  end function complex_tostr


  !> Converts an array of characters to a string.
  function character_arr_tostr(array) result(result_str)
    !> the array to convert
    character(len=*), intent(in)  :: array(:)
    !> returned result, trimmed
    character(:), allocatable :: result_str
    integer :: i

    result_str = "["
    do i = 1, size(array)
      result_str = result_str // trim(array(i)) // ", "
    end do
    result_str = result_str(:len(result_str) - 2) // "]"
  end function character_arr_tostr

end module mod_logging
