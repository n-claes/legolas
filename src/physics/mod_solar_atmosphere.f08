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
  !> derivative of interpolated numberdensity with respect to height
  real(dp), allocatable :: dnh_interp(:)
  !> amount of points used for interpolation, defaults to <tt>ncool</tt>
  integer   :: nbpoints

  private

  public :: set_solar_atmosphere

contains


  subroutine set_solar_atmosphere(rho_field, T_field, grav_field, ninterp)
    use mod_types, only: density_type, temperature_type, gravity_type
    use mod_global_variables, only: ncool
    use mod_grid, only: grid_gauss
    use mod_interpolation, only: lookup_table_value
    use mod_physical_constants, only: gsun_cgs
    use mod_units, only: unit_length, unit_time

    !> type containing the density-related attributes
    type(density_type), intent(inout)     :: rho_field
    !> type containing the temperature-related attributes
    type(temperature_type), intent(inout) :: T_field
    !> type containing the gravity-related attributes
    type(gravity_type), intent(inout)     :: grav_field
    !> points used for interpolation, defaults to <tt>ncool</tt> if not present
    integer, intent(in), optional         :: ninterp

    real(dp)  :: x
    integer   :: i

    nbpoints = ncool
    if (present(ninterp)) then
      nbpoints = ninterp
    end if

    allocate(h_interp(nbpoints))
    allocate(T_interp(nbpoints))
    allocate(nh_interp(nbpoints))
    allocate(dT_interp(nbpoints))
    allocate(dnh_interp(nbpoints))

    write(char_log, int_fmt) nbpoints
    call log_message( &
      "interpolating atmosphere using " // adjustl(trim(char_log)) // "points.", &
      level="info"&
    )
    call create_atmosphere_curves()

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      ! in normalised units numberdensity = density
      rho_field % rho0(i) = lookup_table_value(x, h_interp, nh_interp)
      rho_field % d_rho0_dr(i) = lookup_table_value(x, h_interp, dnh_interp)
      T_field % T0(i) = lookup_table_value(x, h_interp, T_interp)
      T_field % d_T0_dr(i) = lookup_table_value(x, h_interp, dT_interp)
    end do
    grav_field % grav = gsun_cgs / (unit_length / unit_time**2)

    call log_message( &
      "solar atmosphere: rho, T and gravity attributes have been set.", level="info" &
    )

    deallocate(h_interp)
    deallocate(T_interp)
    deallocate(nh_interp)
    deallocate(dT_interp)
    deallocate(dnh_interp)
  end subroutine set_solar_atmosphere


  subroutine create_atmosphere_curves()
    use mod_atmosphere_curves, only: h_alc7, T_alc7, nh_alc7
    use mod_interpolation, only: interpolate_table, get_numerical_derivative
    use mod_units, only: unit_length, unit_temperature, unit_numberdensity

    integer, parameter  :: ntable = size(h_alc7)
    real(dp)  :: log10_height(ntable)

    ! we interpolate on a log10 scale to reduce order of magnitude differences
    log10_height = dlog10(h_alc7)
    ! interpolate log(T) vs log(height)
    call interpolate_table(ntable, log10_height, dlog10(T_alc7), h_interp, T_interp)
    ! interpolate log(nh) vs log(height)
    call interpolate_table(ntable, log10_height, dlog10(nh_alc7), h_interp, nh_interp)

    ! rescale interpolated tables to actual values and normalise
    h_interp = 10.0d0**h_interp * 1.0d5 / unit_length   ! h given in km, scale to cm
    ! rescale other arrays
    T_interp = 10.0d0**T_interp / unit_temperature
    nh_interp = 10.0d0**nh_interp / unit_numberdensity

    ! fill derivatives
    call get_numerical_derivative(h_interp, T_interp, dT_interp)
    call get_numerical_derivative(h_interp, nh_interp, dnh_interp)
  end subroutine create_atmosphere_curves

end module mod_solar_atmosphere