! =============================================================================
!> Module to set a realistic solar atmosphere, using tabulated density and
!! temperature profiles (see <tt>mod_atmosphere_curves</tt>), in Cartesian
!! geometries only.
module mod_solar_atmosphere
  use mod_global_variables, only: dp, gauss_gridpts
  use mod_logging, only: log_message, char_log, int_fmt
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

  private

  public :: set_solar_atmosphere

contains


  !> Sets the density, temperature and gravity attributes of the respective
  !! fields to a realistic solar atmosphere profile. This uses the (Gaussian) grid
  !! to interpolate the correct values from the tables, and should already be set.
  !! This routine first interpolates the temperature and numberdensity table at
  !! <tt>n_interp</tt> resolution, then solves the following ODE for the density:
  !! $$ \rho'(x) = -\frac{T'(x) + g(x)}{T(x)}\rho(x) -
  !!               \frac{B_{02}(x)B'_{02}(x) + B_{03}(x)B'_{03}(x)}{T(x)} $$
  !! using a fifth order Runge-Kutta method, assuming an initial density corresponding
  !! to the first (tabulated) density value in the grid. Integration is done using
  !! <tt>n_interp</tt> values.
  !! @warning   Throws an error if the geometry is not Cartesian.
  subroutine set_solar_atmosphere(rho_field, T_field, B_field, grav_field, n_interp)
    use mod_types, only: density_type, temperature_type, bfield_type, gravity_type
    use mod_global_variables, only: ncool, geometry
    use mod_grid, only: grid_gauss
    use mod_interpolation, only: lookup_table_value, get_numerical_derivative
    use mod_integration, only: integrate_ode_rk
    use mod_physical_constants, only: gsun_cgs, Rsun_cgs
    use mod_units, only: unit_length, unit_time

    !> type containing the density-related attributes
    type(density_type), intent(inout)     :: rho_field
    !> type containing the temperature-related attributes
    type(temperature_type), intent(inout) :: T_field
    !> type containing the magnetic field-related attributes
    type(bfield_type), intent(inout)      :: B_field
    !> type containing the gravity-related attributes
    type(gravity_type), intent(inout)     :: grav_field
    !> points used for interpolation, defaults to <tt>ncool</tt> if not present
    integer, intent(in), optional         :: n_interp

    real(dp)  :: x, gxi, rhoinit
    real(dp)  :: rhovalues(gauss_gridpts), drhovalues(gauss_gridpts)
    real(dp)  :: axvalues(gauss_gridpts), bxvalues(gauss_gridpts)
    integer   :: i

    nbpoints = ncool
    if (present(n_interp)) then
      nbpoints = n_interp
    end if

    if (geometry /= "Cartesian") then
      call log_message( &
        "solar atmosphere can only be set in Cartesian geometries!", &
        level="error" &
      )
      return
    end if

    allocate(h_interp(nbpoints))
    allocate(T_interp(nbpoints))
    allocate(nh_interp(nbpoints))
    allocate(dT_interp(nbpoints))

    call log_message("interpolating solar atmosphere", level="info")
    ! interpolate atmospheric tables
    call create_atmosphere_curves()

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      ! set temperature field with interpolated values
      T_field % T0(i) = lookup_table_value(x, h_interp, T_interp)
      T_field % d_T0_dr(i) = lookup_table_value(x, h_interp, dT_interp)
      ! set the gravitational field
      gxi = gsun_cgs * (Rsun_cgs / (Rsun_cgs + x * unit_length))**2
      grav_field % grav(i) = gxi / (unit_length / unit_time**2)
    end do
    ! set ODE functions
    axvalues = -(T_field % d_T0_dr + grav_field % grav) / (T_field % T0)
    bxvalues = -(B_field % B02 * B_field % d_B02_dr &
               + B_field % B03 * B_field % d_B03_dr) / (T_field % T0)
    ! find initial density value (numberdensity = density in normalised units)
    rhoinit = lookup_table_value(grid_gauss(1), h_interp, nh_interp)

    write(char_log, int_fmt) nbpoints
    call log_message( &
      "solving equilibrium ODE for density... (" &
      // trim(adjustl(char_log)) &
      // " points)", &
      level="info" &
    )
    ! solve differential equation, get derivative as well
    call integrate_ode_rk(grid_gauss, axvalues, bxvalues, rhovalues, &
                          rhoinit, nbpoints, dyvalues=drhovalues)

    ! set density attributes
    rho_field % rho0 = rhovalues
    rho_field % d_rho0_dr = drhovalues

    call log_message( &
      "solar atmosphere: rho, T and gravity attributes have been set", level="info" &
    )

    deallocate(h_interp)
    deallocate(T_interp)
    deallocate(nh_interp)
    deallocate(dT_interp)
  end subroutine set_solar_atmosphere


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

end module mod_solar_atmosphere