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
  end interface str

  !> logical used to (locally) force a prefix override
  logical :: override_prefix_to_false = .false.

  private

  public :: log_message
  public :: print_logo
  public :: print_console_info
  public :: print_whitespace
  public :: str
  public :: exp_fmt, dp_fmt, int_fmt

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

    call log_message("---------------------------------------------")

    ! we temporarily force the use of logging prefix to false
    override_prefix_to_false = .true.

    call log_message("              << Grid settings >>")
    call log_message("geometry           : " // trim(adjustl(geometry)))
    call log_message("grid start         : " // str(x_start))
    call log_message("grid end           : " // str(x_end))
    call log_message("gridpoints (base)  : " // str(gridpts))
    call log_message("gridpoints (Gauss) : " // str(gauss_gridpts))
    call log_message("gridpoints (matrix): " // str(matrix_gridpts))

    call log_message("          << Equilibrium settings >>")
    call log_message("equilibrium    : " // trim(adjustl(equilibrium_type)))
    call log_message("wave number k2 : " // str(k2))
    call log_message("wave number k3 : " // str(k3))
    call log_message("default params : " // str(use_defaults))

    call log_message("            << Physics settings >>")
    call log_message("flow               : " // str(flow))
    call log_message("external gravity   : " // str(external_gravity))

    call log_message("radiative cooling  : " // str(radiative_cooling))
    if (radiative_cooling) then
      call log_message("    cooling curve : " // trim(adjustl(cooling_curve)))
    end if

    call log_message("thermal conduction : " // str(thermal_conduction))
    if (thermal_conduction) then
      if (use_fixed_tc_para) then
        call log_message("    fixed parallel value : " // str(fixed_tc_para_value))
      end if
      if (use_fixed_tc_perp) then
        call log_message("    fixed perpendicular value : " // str(fixed_tc_perp_value))
      end if
    end if
    call log_message("resistivity        : " // str(resistivity))
    if (use_fixed_resistivity) then
      call log_message("    fixed eta value : " // str(fixed_eta_value))
    end if

    call log_message("viscosity          : " // str(viscosity))
    if (viscosity) then
      call log_message("    viscosity value : " // str(viscosity_value))
      call log_message("    viscous heating : " // str(viscous_heating))
    end if

    call log_message("hall mhd           : " // str(hall_mhd))
    if (hall_mhd) then
      call log_message("    electron fraction : " // str(electron_fraction))
      call log_message("    electron inertia  : " // str(elec_inertia))
    end if

    call log_message("            << Solver settings >>")
    call log_message("solver      : " // trim(adjustl(solver)))
    if (solver == "arnoldi") then
      call log_message("ARPACK mode : " // trim(adjustl(arpack_mode)))
      if (arpack_mode == "shift-invert") then
        call log_message("sigma value : " // str(sigma))
      end if
      call log_message("number of eigenvalues : " // str(number_of_eigenvalues))
      call log_message("which eigenvalues     : " // trim(adjustl(which_eigenvalues)))
    end if

    call log_message("            << DataIO settings >>")
    call log_message("datfile name         : " // trim(adjustl(basename_datfile)))
    call log_message("output folder        : " // trim(adjustl(output_folder)))
    call log_message("write matrices       : " // str(write_matrices))
    call log_message("write eigenfunctions : " // str(write_eigenfunctions))
    call log_message( &
      "write derived eigenfunctions : " // str(write_derived_eigenfunctions) &
    )
    call log_message("write eigenfunction subset : " // str(write_eigenfunction_subset))
    if (write_eigenfunction_subset) then
      call log_message("    subset center : " // str(eigenfunction_subset_center))
      call log_message("    subset radius : " // str(eigenfunction_subset_radius))
    end if
    call log_message("---------------------------------------------")
    call print_whitespace(1)
    ! reset override
    override_prefix_to_false = .false.
  end subroutine print_console_info
  ! LCOV_EXCL_STOP


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
