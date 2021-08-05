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
  !> character used as variable to log non-strings
  character(20) :: char_log, char_log2
  !> a convenient "tostring" interface, used for easy console writing
  interface str
    module procedure int_tostr
    module procedure float_tostr
    module procedure complex_tostr
  end interface str

  private

  public :: log_message
  public :: print_logo
  public :: print_console_info
  public :: print_whitespace
  public  :: str
  public :: char_log, char_log2, exp_fmt, dp_fmt, int_fmt

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
  subroutine log_message(msg, level, use_prefix)
    use mod_exceptions, only: raise_exception

    !> the message to print to the console
    character(len=*), intent(in)  :: msg
    !> the level (severity) of the message
    character(len=*), intent(in)  :: level
    !> prefixes message type to string, default is <tt>.true.</tt>
    logical, intent(in), optional :: use_prefix

    ! need a bit more room here, we trim anyway when printing
    character(len=2*str_len) :: msg_painted
    logical                  :: add_prefix

    add_prefix = .true.
    if (present(use_prefix)) then
      add_prefix = use_prefix
    end if

    select case(level)
    case("error")
      call raise_exception(msg)
    case("warning")
      if (logging_level >= 1) then ! LCOV_EXCL_START <we don't print info at testing>
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
      error stop ! LCOV_EXCL_STOP
    end select
  end subroutine log_message


  ! LCOV_EXCL_START <logo is never printed during testing>
  !> Prints the Legolas logo to the console.
  !! The logo is wrapped in 1 whitespace at the top and
  !! two at the bottom. Only for logging level 'warning' (1) and above
  subroutine print_logo()
    use mod_version, only: LEGOLAS_VERSION

    !> array containing the different logo lines
    character(len=str_len) :: logo(11)
    !> whitespace prepended to logo
    character(len=3)       :: spaces_logo = ""
    !> whitespace prepended to versioning
    character(len=57)      :: spaces_v = ""
    integer :: i

    if (logging_level <= 1) then
      return
    end if

    logo(1)  = " __       ________  ________   _______   __          ___     __________ "
    logo(2)  = "|  |     |   ____ \|   ____ \ /   _   \ |  |        /   \   |   ______ \"
    logo(3)  = "|  |     |  |    \/|  |    \/|   / \   ||  |       /  _  \  |  |      \/"
    logo(4)  = "|  |     |  |__    |  |      |  |   |  ||  |      /  / \  \ |  \_______ "
    logo(5)  = "|  |     |   __/   |  | ____ |  |   |  ||  |     /  /   \  \\_______   \"
    logo(6)  = "|  |     |  |      |  | \_  ||  |   |  ||  |    /  /     \  \       |  |"
    logo(7)  = "|  |_____|  |____/\|  |___| ||   \_/   ||  |___/  /  /\___\  \      |  |"
    logo(8)  = "|_______/|________/|________| \_______/ |________/   \_______/      |  |"
    logo(9)  = "                                                        /\__________/  |"
    logo(10) = "Large Eigensystem Generator for One-dimensional pLASmas \______________/"
    logo(11) = spaces_v // "v. " // trim(adjustl(LEGOLAS_VERSION))

    call print_whitespace(1)
    do i = 1, size(logo)
      call paint_string(spaces_logo // trim(logo(i)), "cyan", logo(i))
      write(*, *) logo(i)
    end do
    call print_whitespace(2)
  end subroutine print_logo
  ! LCOV_EXCL_STOP


  ! LCOV_EXCL_START <we don't print info to console when testing>
  !> Prints various console messages showing geometry, grid parameters,
  !! equilibrium parameters etc. Only for logging level "info" or above.
  subroutine print_console_info()
    use mod_global_variables
    use mod_equilibrium_params, only: k2, k3

    if (logging_level <= 1) then
      return
    end if

    call log_message("---------------------------------------------", level="info")
    call log_message( &
      "              << Grid settings >>", level="info", use_prefix=.false. &
    )
    call log_message( &
      "geometry             : " // trim(adjustl(geometry)), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, dp_fmt) x_start
    call log_message( &
      "grid start           : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, dp_fmt) x_end
    call log_message( &
      "grid end             : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, int_fmt) gridpts
    call log_message( &
      "gridpoints (base)    : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, int_fmt) gauss_gridpts
    call log_message( &
      "gridpoints (Gauss)   : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, int_fmt) matrix_gridpts
    call log_message( &
      "gridpoints (matrix)  : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )

    call log_message( &
      "          << Equilibrium settings >>", level="info", use_prefix=.false. &
    )
    call log_message( &
      "selected equilibrium : " // trim(adjustl(equilibrium_type)), &
      level="info", &
      use_prefix=.false. &
    )
    call log_message( &
      "boundary conditions  : " // trim(adjustl(boundary_type)), &
      level="info", &
      use_prefix=.false. &
    )
    write(char_log, dp_fmt) k2
    call log_message( &
      "wave number k2       : " // adjustl(char_log), level="info", use_prefix=.false. &
    )
    write(char_log, dp_fmt) k3
    call log_message( &
      "wave number k3       : " // adjustl(char_log), level="info", use_prefix=.false. &
    )

    call log_message( &
      "            << Physics settings >>", level="info", use_prefix=.false. &
    )
    if (flow) then
      call logical_tostring(flow, char_log)
      call log_message( &
        "flow                 : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (external_gravity) then
      call logical_tostring(external_gravity, char_log)
      call log_message( &
        "external gravity     : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (radiative_cooling) then
      call logical_tostring(radiative_cooling, char_log)
      call log_message( &
        "radiative cooling    : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (thermal_conduction) then
      call logical_tostring(thermal_conduction, char_log)
      call log_message( &
        "thermal conduction   : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (resistivity) then
      call logical_tostring(resistivity, char_log)
      call log_message( &
        "resistivity          : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (viscosity) then
      call logical_tostring(viscosity, char_log)
      call log_message( &
        "viscosity            : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if
    if (viscous_heating) then
      call logical_tostring(viscous_heating, char_log)
      call log_message( &
        "viscous heating      : " // adjustl(char_log), &
        level="info", &
        use_prefix=.false. &
      )
    end if

    call log_message( &
      "            << Solver settings >>", level="info", use_prefix=.false. &
    )
    call log_message( &
      "solver               : " // trim(adjustl(solver)), &
      level="info", &
      use_prefix=.false. &
    )
    if (solver == "arnoldi") then
      call log_message( &
        "ARPACK mode          : " // trim(adjustl(arpack_mode)), &
        level="info", &
        use_prefix=.false. &
      )
      if (arpack_mode == "shift-invert") then
        write(char_log, dp_fmt) real(sigma)
        write(char_log2, '(SP,f18.8,"i")') aimag(sigma)
      call log_message( &
        "sigma value          : " // trim(adjustl(char_log)) &
          // trim(adjustl(char_log2)), &
        level="info", &
        use_prefix=.false. &
      )
      end if
      write(char_log, int_fmt) number_of_eigenvalues
      call log_message( &
        "number of eigenvalues: " // trim(adjustl(char_log)), &
        level="info", &
        use_prefix=.false. &
      )
      call log_message( &
        "which eigenvalues    : " // trim(adjustl(which_eigenvalues)), &
        level="info", &
        use_prefix=.false. &
      )
    end if

    call log_message( &
      "            << DataIO settings >>", level="info", use_prefix=.false. &
    )
    call log_message(&
      "datfile name         : " // trim(adjustl(basename_datfile)), &
      level="info", &
      use_prefix=.false. &
    )
    call log_message( &
      "output folder        : " // trim(adjustl(output_folder)), &
      level="info", &
      use_prefix=.false. &
    )
    call logical_tostring(write_matrices, char_log)
    call log_message( &
      "write matrices       : " // adjustl(char_log), level="info", use_prefix=.false. &
    )
    call logical_tostring(write_eigenfunctions, char_log)
    call log_message( &
      "write eigenfunctions : " // adjustl(char_log), level="info", use_prefix=.false. &
    )
    call logical_tostring(write_derived_eigenfunctions, char_log)
    call log_message( &
      "write derived eigenfunctions : " // adjustl(char_log), &
      level="info", &
      use_prefix=.false. &
    )
    call log_message( &
      "---------------------------------------------", &
      level="info", &
      use_prefix=.false. &
    )
    call print_whitespace(1)
  end subroutine print_console_info
  ! LCOV_EXCL_STOP


  ! LCOV_EXCL_START
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
  ! LCOV_EXCL_STOP


  !> Converts a given integer to a string, the default format is "i8".
  function int_tostr(value, fmt) result(result_str)
    !> integer to convert
    integer, intent(in) :: value
    !> optional format used for writing, default "i8"
    character(len=*), intent(in), optional  :: fmt
    !> return string, made allocatable so it has same length as input
    character(:), allocatable :: result_str
    character(len=20) :: format

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
    character(len=20) :: format

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
    character(len=20) :: format

    if (present(fmt)) then
      format = "(" // trim(fmt) // ")"
    else
      format = "(f18.8)"
    end if
    write(char_log, format) real(value)
    write(char_log2, '(SP,' // format // ',"i")') aimag(value)
    result_str = trim(adjustl(char_log)) // trim(adjustl(char_log2))
  end function complex_tostr


  ! LCOV_EXCL_START <not used during testing>
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
  ! LCOV_EXCL_STOP

end module mod_logging
