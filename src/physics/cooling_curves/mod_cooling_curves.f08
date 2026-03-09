module mod_cooling_curves
  use mod_global_variables, only: dp, str_len
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_logging, only: logger
  use mod_interpolation, only: lookup_table_value, get_numerical_derivative

  use mod_cooling_curve_names
  use mod_data_rosner

  implicit none

  private

  real(dp), allocatable :: curve_T(:)
  real(dp), allocatable :: curve_lambda(:)
  real(dp), allocatable :: curve_dlambdadT(:)

  real(dp) :: minT, maxT

  public :: interpolate_cooling_curves
  public :: is_valid_cooling_curve
  public :: get_rosner_lambdaT, get_rosner_dlambdadT
  public :: get_interpolated_lambdaT, get_interpolated_dlambdadT
  public :: get_cooling_table
  public :: deallocate_cooling_curves

contains

  subroutine interpolate_cooling_curves(settings)
    use mod_interpolation, only: interpolate_table, get_numerical_derivative

    type(settings_t), intent(in) :: settings
    real(dp), allocatable :: table_T(:)
    real(dp), allocatable :: table_lambda(:)
    real(dp) :: unit_temperature, unit_lambdaT, unit_dlambdadT
    integer :: ncool

    call get_cooling_table( &
      name=settings%physics%cooling%get_cooling_curve(), &
      table_T=table_T, &
      table_lambda=table_lambda &
    )
    unit_temperature = settings%units%get_unit_temperature()
    unit_lambdaT = settings%units%get_unit_lambdaT()
    unit_dlambdadT = unit_lambdaT / unit_temperature

    !> @note
    !!     The cooling tables contain dimensionfull values on a logarithmic scale.
    !!     To avoid resampling the table on an unequally spaced temperature grid by doing
    !!     10**T, we interpolate the logarithmic table values on an equally spaced
    !!     T grid in log scale, so we get log10(lambda(T)) and log10(T) values.
    !! @endnote
    ncool = settings%physics%cooling%get_interpolation_points()
    allocate(curve_T(ncool), curve_lambda(ncool), curve_dlambdadT(ncool))
    call interpolate_table( &
      n_interp=ncool, &
      x_table=table_T, &
      y_table=table_lambda, &
      x_interp=curve_T, &
      y_interp=curve_lambda &
    )
    ! rescale to L(T) and normalise
    curve_lambda = 10.0_dp**curve_lambda / unit_lambdaT
    ! normalise logT values: log10(T/Tunit) = log10(T) - log10(Tunit)
    curve_T = curve_T - log10(unit_temperature)
    ! get dlambda(T)/dlogT, hence chain rule: dL(T)/dlogT = (dL(T)/dT) * T
    call get_numerical_derivative(x=curve_T, y=curve_lambda, dy=curve_dlambdadT)
    ! rescale to T (already normalised)
    curve_T = 10.0_dp**curve_T
    ! finally, dL(T)/dT = dL(T)/dlogT * (1 / T)
    curve_dlambdadT = curve_dlambdadT / curve_T

    ! set min/max temperature
    minT = minval(curve_T)
    maxT = maxval(curve_T)

    ! LCOV_EXCL_START
    ! save these curves to a file if we're in debug mode
    if (logger%get_logging_level() >= 3) then
      open( &
        unit=1002, &
        file="debug_coolingcurves", &
        access="stream", &
        status="unknown", &
        action="write" &
      )
      write(1002) size(table_T)
      write(1002) 10.0d0**table_T
      write(1002) 10.0d0**table_lambda
      write(1002) size(curve_T)
      write(1002) curve_T * unit_temperature
      write(1002) curve_lambda * unit_lambdaT
      write(1002) curve_dlambdadT * unit_dlambdadT
      close(1002)
      call logger%debug("cooling curves saved to file 'debug_coolingcurves'")
    end if
    ! LCOV_EXCL_STOP

    deallocate(table_T)
    deallocate(table_lambda)
  end subroutine interpolate_cooling_curves


  pure integer function get_rosner_index(logT0)
    !> dimensionfull log10(T0) value
    real(dp), intent(in) :: logT0
    integer :: j

    get_rosner_index = 1
    if (logT0 > logT_rosner(8)) then
      get_rosner_index = 9
    else
      do j = 1, size(logT_rosner)
        if (logT0 < logT_rosner(j)) then
          get_rosner_index = j
          exit
        end if
      end do
    end if
  end function get_rosner_index


  real(dp) function get_rosner_lambdaT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: logT0, logxi, alpha
    real(dp) :: unit_temperature, unit_lambdaT
    integer :: idx

    unit_temperature = settings%units%get_unit_temperature()
    unit_lambdaT = settings%units%get_unit_lambdaT()
    logT0 = log10(background%temperature%T0(x) * unit_temperature)
    idx = get_rosner_index(logT0)

    logxi = logxi_rosner(idx)
    alpha = alpha_rosner(idx)
    ! lambdaT = xi * T**alpha, so log10(lambdaT) = log10(xi) + alpha * log10(T)
    get_rosner_lambdaT = 10.0_dp**(logxi + alpha * logT0) / unit_lambdaT
  end function get_rosner_lambdaT


  real(dp) function get_rosner_dlambdadT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: logT0, logxi, alpha
    real(dp) :: unit_temperature, unit_lambdaT
    integer :: idx

    unit_temperature = settings%units%get_unit_temperature()
    unit_lambdaT = settings%units%get_unit_lambdaT()
    logT0 = log10(background%temperature%T0(x) * unit_temperature)
    idx = get_rosner_index(logT0)

    logxi = logxi_rosner(idx)
    alpha = alpha_rosner(idx)
    ! dlambdadT = alpha * xi * T**(alpha - 1), and so
    !           = alpha * 10**(logxi + (alpha - 1) * logT0)
    get_rosner_dlambdadT = ( &
      alpha * 10.0_dp**(logxi + (alpha - 1.0_dp) * logT0) &
    ) / (unit_lambdaT / unit_temperature)
  end function get_rosner_dlambdadT


  real(dp) function get_interpolated_lambdaT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: lambdaT, T0
    integer :: pts

    pts = size(curve_T)
    T0 = background%temperature%T0(x)
    if (T0 <= minT) then
      lambdaT = 0.0_dp
    else if (T0 >= maxT) then
      lambdaT = curve_lambda(pts) * sqrt(T0 / maxT)
    else
      lambdaT = lookup_table_value(T0, curve_T, curve_lambda)
    end if
    get_interpolated_lambdaT = lambdaT
  end function get_interpolated_lambdaT


  real(dp) function get_interpolated_dlambdadT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: dlambdadT, T0
    integer :: pts

    pts = size(curve_T)
    T0 = background%temperature%T0(x)
    if (T0 <= minT) then
      dlambdadT = 0.0_dp
    else if (T0 >= maxT) then
      dlambdadT = 0.5_dp * curve_lambda(pts) / sqrt(T0 * maxT)
    else
      dlambdadT = lookup_table_value(T0, curve_T, curve_dlambdadT)
    end if
    get_interpolated_dlambdadT = dlambdadT
  end function get_interpolated_dlambdadT


  logical function is_valid_cooling_curve(name)
    character(len=*), intent(in) :: name

    is_valid_cooling_curve = any(name == KNOWN_CURVES)
    if (.not. is_valid_cooling_curve) then
      call logger%error("unknown cooling curve: " // name)
    end if
  end function is_valid_cooling_curve


  subroutine get_cooling_table(name, table_T, table_lambda)
    use mod_cooling_curve_names
    use mod_data_dalgarno
    use mod_data_dalgarno2
    use mod_data_jccorona
    use mod_data_mlsolar
    use mod_data_rosner
    use mod_data_spex
    use mod_data_spex_enh
    use mod_data_colgan

    character(len=*), intent(in) :: name
    real(dp), intent(out), allocatable :: table_T(:)
    real(dp), intent(out), allocatable :: table_lambda(:)
    integer :: table_n

    select case(name)
    case(JC_CORONA)
      table_n = n_jccorona
      table_T = logT_jccorona
      table_lambda = logL_jccorona
    case(DALGARNO)
      table_n = n_dalgarno
      table_T = logT_dalgarno
      table_lambda = logL_dalgarno
    case(DALGARNO2)
      table_n = n_dalgarno2
      table_T = logT_dalgarno2
      table_lambda = logL_dalgarno2
    case(ML_SOLAR)
      table_n = n_mlsolar
      table_T = logT_mlsolar
      table_lambda = logL_mlsolar
    case(SPEX)
      table_n = n_spex
      table_T = logT_spex
      table_lambda = logL_spex
    case(SPEX_DALGARNO)
      table_n = n_spex + n_spex_enh_dalgarno - 6
      allocate(table_T(table_n))
      allocate(table_lambda(table_n))
      table_T(1:n_spex_enh_dalgarno - 1) = logT_spex_enh_dalgarno( &
        1:n_spex_enh_dalgarno - 1 &
      )
      table_T(n_spex_enh_dalgarno:) = logT_spex(6:n_spex)
      table_lambda(1:n_spex_enh_dalgarno - 1) = logL_spex_enh_dalgarno( &
        1:n_spex_enh_dalgarno - 1 &
      )
      table_lambda(n_spex_enh_dalgarno:) = ( &
        logL_spex(6:n_spex) + log10(L_spex_enh(6:n_spex)) &
      )
    case(COLGAN)
      table_n = n_colgan
      table_T = logT_colgan
      table_lambda = logL_colgan
    case(COLGAN_DM)
      table_n = n_colgan + n_dalgarno2
      allocate(table_T(table_n))
      allocate(table_lambda(table_n))
      table_T(1:n_dalgarno2) = logT_dalgarno2(1:n_dalgarno2)
      table_T(n_dalgarno2+1:) = logT_colgan(1:n_colgan)
      table_lambda(1:n_dalgarno2) = logL_dalgarno2(1:n_dalgarno2)
      table_lambda(n_dalgarno2+1:) = logL_colgan(1:n_colgan)
    end select
  end subroutine get_cooling_table


  subroutine deallocate_cooling_curves()
    if (allocated(curve_T)) deallocate(curve_T)
    if (allocated(curve_lambda)) deallocate(curve_lambda)
    if (allocated(curve_dlambdadT)) deallocate(curve_dlambdadT)
  end subroutine deallocate_cooling_curves

end module mod_cooling_curves
