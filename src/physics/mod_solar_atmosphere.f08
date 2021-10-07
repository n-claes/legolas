! =============================================================================
!> Module to set a realistic solar atmosphere, using tabulated density and
!! temperature profiles (see <tt>mod_atmosphere_curves</tt>), in Cartesian
!! geometries only.
module mod_solar_atmosphere
  use mod_global_variables, only: dp, gauss_gridpts
  use mod_logging, only: log_message, str
  implicit none

  !> interpolated heights from atmosphere tables
  real(dp), allocatable :: h_interp(:)
  !> interpolated temperatures from atmosphere tables
  real(dp), allocatable :: T_interp(:)
  !> interpolated numberdensity from atmosphere tables
  real(dp), allocatable :: nh_interp(:)
  !> derivative of interpolated temperature with respect to height
  real(dp), allocatable :: dT_interp(:)
  !> amount of points used for interpolation, defaults to <tt>ncool</tt>
  integer   :: nbpoints
  !> integrated density profile
  real(dp), allocatable :: rho_values(:)
  !> derivative of integrated density profile
  real(dp), allocatable :: drho_values(:)

  !> chosen profile for B02(x)
  procedure(oned_profile), pointer :: b02_prof
  !> chosen profile for dB02(x)
  procedure(oned_profile), pointer :: db02_prof
  !> chosen profile for B03(x)
  procedure(oned_profile), pointer :: b03_prof
  !> chosen profile for dB03(x)
  procedure(oned_profile), pointer :: db03_prof
  !> chosen profile for g(x)
  procedure(oned_profile), pointer :: gravity_prof

  abstract interface
    function oned_profile(x)
      use mod_global_variables, only: dp
      real(dp), intent(in) :: x(:)
      real(dp)  :: oned_profile(size(x))
    end function oned_profile
  end interface

  private

  public :: set_solar_atmosphere

contains


  !> Sets the density, temperature, gravity and magnetic field attributes of
  !! the respective fields to a realistic solar atmosphere profile.
  !! This routine first interpolates the temperature and numberdensity table at
  !! <tt>n_interp</tt> resolution, then solves the following ODE for the density:
  !! $$ \rho'(x) = -\frac{T'(x) + g(x)}{T(x)}\rho(x) -
  !!               \frac{B_{02}(x)B'_{02}(x) + B_{03}(x)B'_{03}(x)}{T(x)} $$
  !! using a fifth order Runge-Kutta method.
  !! If the optional argument <tt>save_to</tt> is provided then the density profiles
  !! are saved to that file, which can be loaded back in on subsequent runs through
  !! the optional argument <tt>load_from</tt>. The integration is done over the entire
  !! table, the curve is sampled on the Gaussian grid, meaning that grid variations
  !! can all use the same result.
  !! @warning   Throws an error if the geometry is not Cartesian. @endwarning
  subroutine set_solar_atmosphere( &
    f_b02, f_db02, f_b03, f_db03, f_g, n_interp, load_from, save_to &
  )
    use mod_global_variables, only: ncool, geometry
    use mod_grid, only: grid_gauss
    use mod_interpolation, only: lookup_table_value, get_numerical_derivative
    use mod_integration, only: integrate_ode_rk
    use mod_equilibrium, only: rho_field, T_field, B_field, grav_field

    !> function reference for calculation of B02
    procedure(oned_profile), optional :: f_b02
    !> function reference for calculation of B02'
    procedure(oned_profile), optional :: f_db02
    !> function reference for calculation of B03
    procedure(oned_profile), optional :: f_b03
    !> function reference for calculation of B03'
    procedure(oned_profile), optional :: f_db03
    !> function reference for calculation of gravitational profile
    procedure(oned_profile), optional :: f_g
    !> points used for interpolation, defaults to <tt>ncool</tt> if not present
    integer, intent(in), optional :: n_interp
    !> if present loads the (previously) integrated density profile from this
    character(len=*), intent(in), optional  :: load_from
    !> if present, saves the integrated density profile to this
    character(len=*), intent(in), optional  :: save_to

    integer   :: i
    real(dp)  :: x, rhoinit
    real(dp), allocatable  :: axvalues(:), bxvalues(:)
    logical, save :: loaded

    ! check for presence of custom profiles, if not, use default ones
    if (present(f_b02)) then
      if (.not. present(f_db02)) then
        call log_message("solar atmosphere: B02 defined but no dB02", level="error")
        return
      end if
      b02_prof => f_b02
      db02_prof => f_db02
    else
      b02_prof => default_b02_profile
      db02_prof => default_db02_profile
    end if
    if (present(f_b03)) then
      if (.not. present(f_db03)) then
        call log_message("solar atmosphere: B03 defined but no dB03", level="error")
        return
      end if
      b03_prof => f_b03
      db03_prof => f_db03
    else
      b03_prof => default_b03_profile
      db03_prof => default_db03_profile
    end if
    if (present(f_g)) then
      gravity_prof => f_g
    else
      gravity_prof => default_gravity_profile
    end if

    nbpoints = ncool
    if (present(n_interp)) then
      nbpoints = n_interp
    end if

    if (geometry /= "Cartesian") then
      call log_message( &
        "solar atmosphere can only be set in Cartesian geometries!", level="error" &
      )
      return
    end if

    call log_message("setting solar atmosphere...", level="info")
    ! load pre-existing profiles from file
    if (present(load_from)) then
      ! this allocates and sets rho_values and drho_values
      call load_profile_from_file(load_from, loading_ok=loaded)
      if (.not. loaded) then
        call log_message("solar atmosphere: profile loading failed!", level="error")
        return
      end if
    else
      allocate(rho_values(nbpoints), drho_values(nbpoints))
      loaded = .false.
    end if

    ! interpolate atmospheric tables
    allocate(h_interp(nbpoints), T_interp(nbpoints), dT_interp(nbpoints))
    allocate(nh_interp(nbpoints))
    ! curves are normalised on return
    call create_atmosphere_curves()

    ! only do integration if no profile was loaded
    if (.not. loaded) then
      ! fill ODE functions, use high resolution arrays
      allocate(axvalues(nbpoints), bxvalues(nbpoints))
      axvalues = -(dT_interp + gravity_prof(h_interp)) / (T_interp)
      bxvalues = ( &
        -( &
          b02_prof(h_interp) * db02_prof(h_interp) &
          + b03_prof(h_interp) * db03_prof(h_interp) &
        ) &
        / T_interp &
      )
      ! set initial density value (numberdensity = density in normalised units)
      rhoinit = nh_interp(1)
      ! solve differential equation
      call log_message( &
        "solving equilibrium ODE for density (" // str(nbpoints) // " points)", &
        level="info", &
        use_prefix=.false. &
      )
      ! do integration
      call integrate_ode_rk( &
        h_interp, axvalues, bxvalues, nbpoints, rhoinit, rho_values, adaptive=.false. &
      )
      call get_numerical_derivative(h_interp, rho_values, drho_values)
      ! these are no longer needed
      deallocate(axvalues, bxvalues)
    end if

    ! save to file if asked but don't load and save at the same time
    if (present(save_to) .and. .not. loaded) then
      call save_profile_to_file(save_to)
    end if

    ! set the various equilibrium attributes
    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      ! density
      rho_field % rho0(i) = lookup_table_value(x, h_interp, rho_values)
      ! temperature
      T_field % T0(i) = lookup_table_value(x, h_interp, T_interp)
      T_field % d_T0_dr(i) = lookup_table_value(x, h_interp, dT_interp)
    end do
    ! gravitational field
    grav_field % grav = gravity_prof(grid_gauss)
    ! magnetic fields
    B_field % B02 = b02_prof(grid_gauss)
    B_field % d_B02_dr = db02_prof(grid_gauss)
    B_field % B03 = b03_prof(grid_gauss)
    B_field % d_B03_dr = db03_prof(grid_gauss)
    B_field % B0 = sqrt(B_field % B02**2 + B_field % B03**2)

    ! set density derivative USING the differential equation. The difference between
    ! drho_interp and doing it like this is due to numerical accuracy, but the
    ! differences are quite acceptable (about 0.01% or lower). Setting this ensures
    ! that the equilibrium is satisfied up to machine precision.
    rho_field % d_rho0_dr = ( &
      - (T_field % d_T0_dr + grav_field % grav) * rho_field % rho0 &
      + (B_field % B02 * B_field % d_B02_dr + B_field % B03 * B_field % d_B03_dr) &
    ) / T_field % T0

    call log_message( &
      "rho, T, B and gravity attributes have been set", &
      level="info", &
      use_prefix=.false. &
    )

    deallocate(h_interp, T_interp, nh_interp, dT_interp)
    deallocate(rho_values, drho_values)
  end subroutine set_solar_atmosphere


  !> Interpolates the atmospheric tables to the desired resolution. The temperature
  !! derivative is obtained numerically.
  subroutine create_atmosphere_curves()
    use mod_atmosphere_curves, only: h_alc7, T_alc7, nh_alc7
    use mod_interpolation, only: interpolate_table, get_numerical_derivative
    use mod_units, only: unit_length, unit_temperature, unit_numberdensity
    use mod_global_variables, only: logging_level

    ! interpolate nh vs height
    call interpolate_table(nbpoints, h_alc7, nh_alc7, h_interp, nh_interp)
    ! interpolate T vs height
    call interpolate_table(nbpoints, h_alc7, T_alc7, h_interp, T_interp)

    ! rescale interpolated tables to actual values and normalise
    h_interp = h_interp * 1.0d5 / unit_length ! height is in km, so scale to cm first
    T_interp = T_interp / unit_temperature
    nh_interp = nh_interp / unit_numberdensity

    ! find temperature derivative
    call get_numerical_derivative(h_interp, T_interp, dT_interp)
    ! save these curves to a file if we're in debug mode
    if (logging_level >= 3) then ! LCOV_EXCL_START
      open( &
        unit=1002, &
        file="debug_atmocurves", &
        access="stream", &
        status="unknown", &
        action="write" &
      )
      write(1002) size(h_alc7)
      write(1002) h_alc7 * 1.0d5  ! save dimensionfull in cm
      write(1002) T_alc7
      write(1002) nh_alc7
      write(1002) size(h_interp)
      write(1002) h_interp * unit_length
      write(1002) T_interp * unit_temperature
      write(1002) nh_interp * unit_numberdensity
      write(1002) dT_interp * (unit_temperature / unit_length)
      close(1002)
      call log_message( &
        "atmo curves saved to file 'debug_atmocurves'", level="debug" &
      )
    end if ! LCOV_EXCL_STOP
  end subroutine create_atmosphere_curves


  !> Saves the density and density derivatives to the given filename. These can be
  !! used later on to set the values instead of solving the differential equation.
  subroutine save_profile_to_file(filename)
    use mod_units, only: unit_length, unit_temperature, unit_magneticfield, unit_density

    !> values are saved to this filename
    character(len=*), intent(in)  :: filename
    integer   :: unit

    unit = 1001
    open(unit=unit, file=filename, access="stream", status="unknown", action="write")
    write(unit) size(rho_values)
    ! write dimensions, so normalisations happen correctly when loading
    write(unit) unit_length, unit_temperature, unit_magneticfield, unit_density
    ! next we write the B and g values to file, these are checked when loading
    write(unit) h_interp
    write(unit) b02_prof(h_interp), db02_prof(h_interp)
    write(unit) b03_prof(h_interp), db03_prof(h_interp)
    write(unit) gravity_prof(h_interp)
    write(unit) rho_values, drho_values
    close(unit)
    call log_message( &
      "integrated density profiles saved to " // trim(adjustl(filename)), &
      level="info", &
      use_prefix=.false. &
    )
  end subroutine save_profile_to_file


  !> Loads a previously calculated profile and uses that to set the resolution,
  !! density and density derivatives.
  subroutine load_profile_from_file(filename, loading_ok)
    use mod_check_values, only: is_equal, is_constant
    use mod_units, only: unit_length, unit_temperature, unit_magneticfield, unit_density

    !> values are loaded from this file
    character(len=*), intent(in) :: filename
    logical, intent(out) :: loading_ok

    integer :: unit, resolution
    real(dp)  :: length_file, temperature_file, magneticfield_file, density_file
    logical   :: b02_cte, b03_cte, db02_zero, db03_zero
    character(len=:), allocatable :: prof_names
    real(dp), allocatable :: hfile(:), profile(:)

    loading_ok = .false.
    unit = 1002
    open(unit=unit, file=filename, access="stream", status="old", action="read")

    read(unit) resolution
    if (.not. resolution == nbpoints) then
      call log_message( &
        "set resolution (" // str(nbpoints) // &
        ") has been overriden by resolution from file (" // str(resolution) // ")", &
        level="warning" &
      )
    end if
    nbpoints = resolution
    call log_message( &
      "restoring density profiles from " // filename // " [" // str(nbpoints) // "]", &
      level="info", &
      use_prefix=.false. &
    )

    ! check normalisations
    read(unit) length_file, temperature_file, magneticfield_file, density_file
    if (.not. is_equal(length_file, unit_length)) then
      call log_message( &
        "profile inconsistency: length units do not match! Got " // &
        str(length_file) // " but expected " // str(unit_length), &
        level="warning" &
      )
      close(unit)
      return
    end if
    if (.not. is_equal(temperature_file, unit_temperature)) then
      call log_message( &
        "profile inconsistency: temperature units do not match! Got " // &
        str(temperature_file) // " but expected " // str(unit_temperature), &
        level="warning" &
      )
      close(unit)
      return
    end if
    if (.not. is_equal(magneticfield_file, unit_magneticfield)) then
      call log_message( &
        "profile inconsistency: magnetic units do not match! Got " // &
        str(magneticfield_file) // " but expected " // str(unit_magneticfield), &
        level="warning" &
      )
      close(unit)
      return
    end if
    if (.not. is_equal(density_file, unit_density)) then
      call log_message( &
        "profile inconsistency: density units do not match! Got " // &
        str(length_file) // " but expected " // str(unit_length), &
        level="warning" &
      )
      close(unit)
      return
    end if

    ! Here we check if the profiles that are provided correspond to the ones that
    ! were used to save the integrated profile. Since Fortran cannot save function
    ! statements themselves, the B and g values were saved to the file instead and we
    ! compare them at the same resolution here.
    allocate(hfile(nbpoints), profile(nbpoints))
    allocate(prof_names, mold="")

    ! check B02
    read(unit) hfile, profile
    if (.not. all(is_equal(profile, b02_prof(hfile)))) then
      prof_names = trim(prof_names // " B02")
    end if
    b02_cte = (is_constant(b02_prof(hfile)) .and. is_constant(profile))
    ! check dB02
    read(unit) profile
    if (.not. all(is_equal(profile, db02_prof(hfile)))) then
      prof_names = trim(prof_names // " dB02")
    end if
    db02_zero = ( &
      all(is_equal(db02_prof(hfile), 0.0d0)) .and. all(is_equal(profile, 0.0d0)) &
    )
    ! check B03
    read(unit) profile
    if (.not. all(is_equal(profile, b03_prof(hfile)))) then
      prof_names = trim(prof_names // " B03")
    end if
    b03_cte = (is_constant(b03_prof(hfile)) .and. is_constant(profile))
    ! check dB03
    read(unit) profile
    if (.not. all(is_equal(profile, db03_prof(hfile)))) then
      prof_names = trim(prof_names // " dB03")
    end if
    db03_zero = ( &
      all(is_equal(db03_prof(hfile), 0.0d0)) .and. all(is_equal(profile, 0.0d0)) &
    )

    ! if both B02 and B03 are constant then values do not matter for profile, so skip
    if (b02_cte .and. dB02_zero .and. b03_cte .and. dB03_zero) then
      deallocate(prof_names)
      allocate(prof_names, mold="")
      call log_message( &
        "load profile: B02 and B03 are constant, skipping B0 checks", level="debug" &
      )
    end if

    ! check gravity
    read(unit) profile
    if (.not. all(is_equal(profile, gravity_prof(hfile)))) then
      prof_names = trim(prof_names // " gravity")
    end if
    if (.not. prof_names == "") then
      call log_message( &
        "profile inconsistency in [" // prof_names // " ] when loading from file", &
        level="warning" &
      )
      close(unit)
      return
    end if
    deallocate(hfile, profile)

    ! if everything checks out, read in density profile
    allocate(rho_values(nbpoints), drho_values(nbpoints))
    read(unit) rho_values, drho_values
    close(unit)
    loading_ok = .true.
  end subroutine load_profile_from_file


  !> Sets the default profile for B02, taken to be zero.
  function default_b02_profile(x) result(b02)
    !> x-values to evaluate the profile in
    real(dp), intent(in)  :: x(:)
    !> resulting B02 profile
    real(dp)  :: b02(size(x))

    b02 = 0.0d0
  end function default_b02_profile


  !> Sets the default profile for dB02, taken to be zero.
  function default_db02_profile(x) result(db02)
    !> x-values to evaluate the profile in
    real(dp), intent(in)  :: x(:)
    !> resulting dB02 profile
    real(dp)  :: db02(size(x))

    db02 = 0.0d0
  end function default_db02_profile


  !> Sets the default profile for B03, taken to be a uniform field of 10 Gauss.
  function default_b03_profile(x) result(b03)
    use mod_units, only: unit_magneticfield

    !> x-values to evaluate the profile in
    real(dp), intent(in)  :: x(:)
    !> resulting B03 profile
    real(dp)  :: b03(size(x))

    b03 = 10.0d0 / unit_magneticfield
  end function default_b03_profile


  !> Sets the default profile for dB03, taken to be zero.
  function default_db03_profile(x) result(db03)
    !> x-values to evaluate the profile in
    real(dp), intent(in)  :: x(:)
    !> resulting dB03 profile
    real(dp)  :: db03(size(x))

    db03 = 0.0d0
  end function default_db03_profile


  !> Sets the default profile for the gravitational field, taken to be
  !! $$ g(x) = g_\odot \left(\frac{R_\odot}{(R_\odot + x)}\right)^2 $$
  function default_gravity_profile(x) result(grav)
    use mod_units, only: unit_length, unit_time
    use mod_physical_constants, only: gsun_cgs, Rsun_cgs

    !> x-values to evaluate the profile in
    real(dp), intent(in)  :: x(:)
    !> resulting gravitational profile
    real(dp)  :: grav(size(x))

    grav = ( &
      gsun_cgs &
      * (Rsun_cgs / (Rsun_cgs + x * unit_length))**2 &
      / (unit_length / unit_time**2) &
    )
  end function default_gravity_profile

end module mod_solar_atmosphere
