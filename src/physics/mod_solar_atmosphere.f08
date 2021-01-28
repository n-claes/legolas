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

  abstract interface
    function b_func(x)
      use mod_global_variables, only: dp
      real(dp), intent(in) :: x(:)
      real(dp)  :: b_func(size(x))
    end function b_func
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
    b02_func, db02_func, b03_func, db03_func, n_interp, load_from, save_to &
  )
    use mod_global_variables, only: ncool, geometry
    use mod_grid, only: grid_gauss
    use mod_interpolation, only: lookup_table_value, get_numerical_derivative
    use mod_integration, only: integrate_ode_rk
    use mod_physical_constants, only: gsun_cgs, Rsun_cgs
    use mod_units, only: unit_length, unit_time
    use mod_equilibrium, only: rho_field, T_field, B_field, grav_field

    !> function reference for calculation of B02
    procedure(b_func) :: b02_func
    !> function reference for calculation of B02'
    procedure(b_func) :: db02_func
    !> function reference for calculation of B03
    procedure(b_func) :: b03_func
    !> function reference for calculation of B03'
    procedure(b_func) :: db03_func
    !> points used for interpolation, defaults to <tt>ncool</tt> if not present
    integer, intent(in), optional :: n_interp
    !> if present loads the (previously) integrated density profile from this
    character(len=*), intent(in), optional  :: load_from
    !> if present, saves the integrated density profile to this
    character(len=*), intent(in), optional  :: save_to

    integer   :: i
    real(dp)  :: x, rhoinit
    real(dp), allocatable :: g_interp(:)
    real(dp), allocatable  :: axvalues(:), bxvalues(:)
    logical, save :: loaded

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
    loaded = .false.
    if (present(load_from)) then
      ! this allocates and sets rho_values and drho_values
      call load_profile_from_file(load_from)
      loaded = .true.
    else
      allocate(rho_values(nbpoints), drho_values(nbpoints))
    end if

    ! interpolate atmospheric tables
    allocate(h_interp(nbpoints), T_interp(nbpoints), dT_interp(nbpoints))
    allocate(nh_interp(nbpoints))
    ! curves are normalised on return
    call create_atmosphere_curves()

    ! create gravitational field
    allocate(g_interp(nbpoints))
    g_interp = ( &
      gsun_cgs &
      * (Rsun_cgs / (Rsun_cgs + h_interp * unit_length))**2 &
      / (unit_length / unit_time**2) &
    )

    ! only do integration if no profile was loaded
    if (.not. loaded) then
      ! fill ODE functions, use high resolution arrays
      allocate(axvalues(nbpoints), bxvalues(nbpoints))
      axvalues = -(dT_interp + g_interp) / (T_interp)
      bxvalues = ( &
        -( &
          b02_func(h_interp) * db02_func(h_interp) &
          + b03_func(h_interp) * db03_func(h_interp) &
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

    ! save to file if asked
    if (present(save_to)) then
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
      ! gravitational field
      grav_field % grav(i) = lookup_table_value(x, h_interp, g_interp)
    end do
    ! magnetic fields
    B_field % B02 = b02_func(grid_gauss)
    B_field % d_B02_dr = db02_func(grid_gauss)
    B_field % B03 = b03_func(grid_gauss)
    B_field % d_B03_dr = db03_func(grid_gauss)
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
    deallocate(g_interp)
    deallocate(rho_values, drho_values)
  end subroutine set_solar_atmosphere


  !> Interpolates the atmospheric tables to the desired resolution. The temperature
  !! derivative is obtained numerically.
  subroutine create_atmosphere_curves()
    use mod_atmosphere_curves, only: h_alc7, T_alc7, nh_alc7
    use mod_interpolation, only: interpolate_table, get_numerical_derivative
    use mod_units, only: unit_length, unit_temperature, unit_numberdensity

    ! interpolate T vs height
    call interpolate_table(nbpoints, h_alc7, T_alc7, h_interp, T_interp)
    ! interpolate nh vs height
    call interpolate_table(nbpoints, h_alc7, nh_alc7, h_interp, nh_interp)

    ! rescale interpolated tables to actual values and normalise
    h_interp = h_interp * 1.0d5 / unit_length ! height is in km, so scale to cm first
    T_interp = T_interp / unit_temperature
    nh_interp = nh_interp / unit_numberdensity

    ! find temperature derivative
    call get_numerical_derivative(h_interp, T_interp, dT_interp)
  end subroutine create_atmosphere_curves


  !> Saves the density and density derivatives to the given filename. These can be
  !! used later on to set the values instead of solving the differential equation.
  subroutine save_profile_to_file(filename)
    !> values are saved to this filename
    character(len=*), intent(in)  :: filename
    integer   :: unit

    unit = 1001
    open(unit=unit, file=filename, access="stream", status="unknown", action="write")
    write(unit) size(rho_values), rho_values, drho_values
    close(unit)
    call log_message( &
      "integrated density profiles saved to " // trim(adjustl(filename)), &
      level="info", &
      use_prefix=.false. &
    )
  end subroutine save_profile_to_file


  !> Loads a previously calculated profile and uses that to set the resolution,
  !! density and density derivatives.
  subroutine load_profile_from_file(filename)
    !> values are loaded from this file
    character(len=*), intent(in)  :: filename
    integer   :: unit, resolution

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
    allocate(rho_values(nbpoints), drho_values(nbpoints))
    read(unit) rho_values, drho_values
    close(unit)
    call log_message( &
      "restored density profiles from " // trim(adjustl(filename)) // &
      " (" // str(nbpoints) // " pts)", &
      level="info", &
      use_prefix=.false. &
    )
  end subroutine load_profile_from_file

end module mod_solar_atmosphere
