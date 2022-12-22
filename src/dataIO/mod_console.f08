module mod_console
  use mod_global_variables
  use mod_logging, only: log_message, override_prefix_to_false, str, exp_fmt
  implicit none

  private

  public :: print_logo
  public :: print_console_info
  public :: print_whitespace

contains

  ! LCOV_EXCL_START (we don't test printing statements)

  !> Prints the Legolas logo to the console.
  !! The logo is wrapped in 1 whitespace at the top and
  !! two at the bottom. Only for logging level 'warning' (1) and above
  subroutine print_logo()
    use mod_version, only: LEGOLAS_VERSION
    use mod_painting, only: paint_string

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

    logo(1) = " __       ________  ________   _______   __          ___     __________ "
    logo(2) = "|  |     |   ____ \|   ____ \ /   _   \ |  |        /   \   |   ______ \"
    logo(3) = "|  |     |  |    \/|  |    \/|   / \   ||  |       /  _  \  |  |      \/"
    logo(4) = "|  |     |  |__    |  |      |  |   |  ||  |      /  / \  \ |  \_______ "
    logo(5) = "|  |     |   __/   |  | ____ |  |   |  ||  |     /  /   \  \\_______   \"
    logo(6) = "|  |     |  |      |  | \_  ||  |   |  ||  |    /  /     \  \       |  |"
    logo(7) = "|  |_____|  |____/\|  |___| ||   \_/   ||  |___/  /  /\___\  \      |  |"
    logo(8) = "|_______/|________/|________| \_______/ |________/   \_______/      |  |"
    logo(9) = "                                                        /\__________/  |"
    logo(10) ="Large Eigensystem Generator for One-dimensional pLASmas \______________/"
    logo(11) = spaces_v // "v. " // trim(adjustl(LEGOLAS_VERSION))

    call print_whitespace(1)
    do i = 1, size(logo)
      call paint_string(spaces_logo // trim(logo(i)), "cyan", logo(i))
      write(*, *) logo(i)
    end do
    call print_whitespace(2)
  end subroutine print_logo


  ! LCOV_EXCL_START <we don't print info to console when testing>
  !> Prints various console messages showing geometry, grid parameters,
  !! equilibrium parameters etc. Only for logging level "info" or above.
  subroutine print_console_info(settings)
    use mod_global_variables
    use mod_equilibrium_params, only: k2, k3
    use mod_settings, only: settings_t

    type(settings_t), intent(in) :: settings

    if (logging_level <= 1) then
      return
    end if

    call log_message("---------------------------------------------")

    ! we temporarily force the use of logging prefix to false
    override_prefix_to_false = .true.

    call log_message("              << Grid settings >>")
    call log_message("geometry           : " // settings%grid%get_geometry())
    call log_message("grid start         : " // str(settings%grid%get_grid_start()))
    call log_message("grid end           : " // str(settings%grid%get_grid_end()))
    call log_message("gridpoints (base)  : " // str(settings%grid%get_gridpts()))
    call log_message("gridpoints (Gauss) : " // str(settings%grid%get_gauss_gridpts()))
    call log_message("gridpoints (matrix): " // str(settings%dims%get_dim_matrix()))

    call log_message("          << Equilibrium settings >>")
    call log_message( &
      "equilibrium    : " // settings%equilibrium%get_equilibrium_type() &
    )
    call log_message("wave number k2 : " // str(k2))
    call log_message("wave number k3 : " // str(k3))
    call log_message("default params : " // str(settings%equilibrium%use_defaults))

    call log_message("            << Physics settings >>")
    call log_message("flow               : " // str(settings%physics%flow%is_enabled()))
    call log_message("external gravity   : " // str(external_gravity))

    call log_message( &
      "radiative cooling  : " // str(settings%physics%cooling%is_enabled()) &
    )
    if (settings%physics%cooling%is_enabled()) then
      call log_message( &
        "    cooling curve : " &
        // trim(adjustl(settings%physics%cooling%get_cooling_curve())) &
      )
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
      call log_message("    by substitution   : " // str(hall_substitution))
      call log_message("    electron fraction : " // str(electron_fraction))
      call log_message("    electron inertia  : " // str(elec_inertia))
    end if

    call log_message("            << Solver settings >>")
    call log_message("solver      : " // settings%solvers%get_solver())
    if (settings%solvers%get_solver() == "arnoldi") then
      call log_message("ARPACK mode : " // settings%solvers%get_arpack_mode())
      if (settings%solvers%get_arpack_mode() == "shift-invert") then
        call log_message("sigma : " // str(settings%solvers%sigma))
      end if
      call log_message( &
        "number of eigenvalues : " // str(settings%solvers%number_of_eigenvalues) &
      )
      call log_message( &
        "which eigenvalues     : " // settings%solvers%which_eigenvalues &
      )
      call log_message("maxiter               : " // str(settings%solvers%maxiter))
      call log_message( &
        "tolerance             : " // str(settings%solvers%tolerance, exp_fmt) &
      )
    end if
    if (settings%solvers%get_solver() == "inverse-iteration") then
      call log_message("sigma     : " // str(settings%solvers%sigma))
      call log_message("maxiter   : " // str(settings%solvers%maxiter))
      call log_message("tolerance : " // str(settings%solvers%tolerance, exp_fmt))
    end if

    call log_message("            << DataIO settings >>")
    call log_message("datfile name         : " // settings%io%get_basename_datfile())
    call log_message("output folder        : " // settings%io%get_output_folder())
    call log_message("write matrices       : " // str(settings%io%write_matrices))
    call log_message("write eigenvectors   : " // str(settings%io%write_eigenvectors))
    call log_message("write residuals      : " // str(settings%io%write_residuals))
    call log_message("write eigenfunctions : " // str(settings%io%write_eigenfunctions))
    call log_message( &
      "write derived eigenfunctions : " &
      // str(settings%io%write_derived_eigenfunctions) &
    )
    call log_message( &
      "write eigenfunction subset : " // str(settings%io%write_ef_subset) &
    )
    if (settings%io%write_ef_subset) then
      call log_message("    subset center : " // str(settings%io%ef_subset_center))
      call log_message("    subset radius : " // str(settings%io%ef_subset_radius))
    end if
    call log_message("---------------------------------------------")
    call print_whitespace(1)
    ! reset override
    override_prefix_to_false = .false.
  end subroutine print_console_info


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

end module mod_console
