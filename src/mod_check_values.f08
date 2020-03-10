module mod_check_values
  use mod_global_variables, only: dp, dp_LIMIT, gauss_gridpts
  implicit none

  private

  interface check_small_values
    module procedure small_values_real_array
    module procedure small_values_complex_array
    module procedure small_values_real_matrix
    module procedure small_values_complex_matrix
  end interface check_small_values

  interface check_nan_values
    module procedure check_nan_values_density
    module procedure check_nan_values_temperature
    module procedure check_nan_values_bfield
    module procedure check_nan_values_velocity
    module procedure check_nan_values_gravity
  end interface check_nan_values

  public :: check_small_values
  public :: check_negative_array
  public :: check_equilibrium_conditions
  public :: check_nan_values

contains

  subroutine small_values_real_array(array)
    real(dp), intent(inout) :: array(:)
    integer                 :: i

    do i = 1, size(array)
      if (abs(array(i)) < dp_LIMIT) then
        array(i) = 0.0d0
      end if
    end do
  end subroutine small_values_real_array

  subroutine small_values_complex_array(array)
    complex(dp), intent(inout) :: array(:)
    integer                    :: i
    real(dp)                   :: a_real, a_imag

    do i = 1, size(array)
      a_real = real(array(i))
      a_imag = aimag(array(i))

      if (abs(a_real) < dp_LIMIT) then
        a_real = 0.0d0
      end if
      if (abs(a_imag) < dp_LIMIT) then
        a_imag = 0.0d0
      end if

      array(i) = cmplx(a_real, a_imag, kind=dp)
    end do
  end subroutine small_values_complex_array

  subroutine small_values_real_matrix(matrix)
    real(dp), intent(inout)    :: matrix(:, :)
    integer                    :: i, j

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        if (abs(matrix(i, j)) < dp_LIMIT) then
          matrix(i, j) = 0.0d0
        end if
      end do
    end do
  end subroutine small_values_real_matrix

  subroutine small_values_complex_matrix(matrix)
    complex(dp), intent(inout) :: matrix(:, :)
    integer                    :: i, j
    real(dp)                   :: a_real, a_imag

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        a_real = real(matrix(i, j))
        a_imag = aimag(matrix(i, j))

        if (abs(a_real) < dp_LIMIT) then
          a_real = 0.0d0
        end if
        if (abs(a_imag) < dp_LIMIT) then
          a_imag = 0.0d0
        end if

        matrix(i, j) = cmplx(a_real, a_imag, kind=dp)
      end do
    end do
  end subroutine small_values_complex_matrix

  subroutine check_negative_array(array, variable_name)
    real(dp), intent(in)          :: array(:)
    character(len=*), intent(in)  :: variable_name
    integer                       :: i

    do i = 1, size(array)
      if (array(i) < 0.0d0) then
        write(*, *) "WARNING: ", trim(variable_name), " is negative somewhere!"
        error stop
      end if
    end do
  end subroutine check_negative_array


  subroutine check_equilibrium_conditions(rho_field, T_field, B_field, v_field, grav_field)
    use mod_types, only: density_type, temperature_type, bfield_type, velocity_type, gravity_type
    use mod_global_variables, only: geometry, dp_LIMIT, x_start
    use mod_grid, only: eps_grid, d_eps_grid_dr
    use mod_equilibrium_params, only: k2

    type(density_type), intent(in)      :: rho_field
    type(temperature_type), intent(in)  :: T_field
    type(bfield_type), intent(in)       :: B_field
    type(velocity_type), intent(in)     :: v_field
    type(gravity_type), intent(in)      :: grav_field

    real(dp)    :: rho, drho, B02, dB02, B03, dB03, T0, dT0, grav, v02, v03
    real(dp)    :: eps, d_eps, axis_limit
    real(dp)    :: eq_cond(gauss_gridpts)
    integer     :: i, k2_int

    k2_int = int(k2)
    axis_limit = 1.0d-2

    if (geometry == 'cylindrical') then
      ! in cylindrical geometry, m should be an integer
      if (abs(k2_int - k2) > dp_LIMIT) then
        write(*, *) "k2 value: ", k2
        error stop "cylindrical geometry defined but k2 (mode number m) is not an integer!"
      end if

      ! check if relevant values are zero on-axis, if applicable
      if (x_start <= 1.0d-5) then
        if (abs(B_field % B02(1)) > axis_limit) then
          write(*, *) "B_theta(0) value: ", B_field % B02(1)
          error stop "B_theta(0) is non-zero on axis!"
        end if
        if (abs(B_field % d_B03_dr(1)) > axis_limit) then
          write(*, *) "dBz/dr(0) value: ", B_field % d_B03_dr(1)
          error stop "dBz/dr(0) is non-zero on axis!"
        end if
        if (abs(v_field % v02(1)) > axis_limit) then
          write(*, *) "v_theta(0) value: ", v_field % v02(1)
          error stop "v_theta(0) is non-zero on axis!"
        end if
        if (abs(v_field % d_v03_dr(1)) > axis_limit) then
          write(*, *) "dvz/dr(0) value: ", v_field % d_v03_dr(1)
          error stop "dvz/dr(0) is non-zero on axis!"
        end if
      end if
    end if

    ! check if equilibrium conditions are met
    do i = 1, gauss_gridpts
      rho = rho_field % rho0(i)
      drho = rho_field % d_rho0_dr(i)
      B02 = B_field % B02(i)
      B03 = B_field % B03(i)
      dB02 = B_field % d_B02_dr(i)
      dB03 = B_field % d_B03_dr(i)
      T0 = T_field % T0(i)
      dT0 = T_field % d_T0_dr(i)
      grav = grav_field % grav(i)
      v02 = v_field % v02(i)
      v03 = v_field % v03(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i) = drho * T0 + rho * dT0 + B02 * dB02 + B03 * dB03 + rho * grav - (d_eps/eps) * (rho * v02**2 - B02**2)
      if (eq_cond(i) > dp_LIMIT) then
        error stop "equilibrium conditions not met!"
      end if
    end do
  end subroutine check_equilibrium_conditions


  subroutine stop_if_nan(array, array_name)
    use, intrinsic  :: ieee_arithmetic, only: ieee_is_nan

    real(dp), intent(in)  :: array(gauss_gridpts)
    character(len=*), intent(in)  :: array_name
    integer :: i

    do i = 1, gauss_gridpts
      if (ieee_is_nan(array(i))) then
        write(*, *) "Checking ", array_name
        error stop "NaN encountered!"
      end if
    end do
  end subroutine stop_if_nan



  subroutine check_nan_values_density(rho_field)
    use mod_types, only: density_type

    type(density_type), intent(in)  :: rho_field

    call stop_if_nan(rho_field % rho0, "rho0")
    call stop_if_nan(rho_field % d_rho0_dr, "drho0")
  end subroutine check_nan_values_density


  subroutine check_nan_values_temperature(T_field)
    use mod_types, only: temperature_type

    type(temperature_type), intent(in)  :: T_field

    call stop_if_nan(T_field % T0, "T0")
    call stop_if_nan(T_field % d_T0_dr, "dT0")
  end subroutine check_nan_values_temperature


  subroutine check_nan_values_bfield(B_field)
    use mod_types, only: bfield_type

    type(bfield_type), intent(in) :: B_field

    call stop_if_nan(B_field % B02, "B02")
    call stop_if_nan(B_field % B03, "B03")
    call stop_if_nan(B_field % B0, "B0")
    call stop_if_nan(B_field % d_B02_dr, "dB02")
    call stop_if_nan(B_field % d_B03_dr, "dB03")
  end subroutine check_nan_values_bfield


  subroutine check_nan_values_velocity(v_field)
    use mod_types, only: velocity_type

    type(velocity_type), intent(in) :: v_field

    call stop_if_nan(v_field % v02, "v02")
    call stop_if_nan(v_field % v03, "v03")
    call stop_if_nan(v_field % d_v02_dr, "dv02")
    call stop_if_nan(v_field % d_v03_dr, "dv03")
  end subroutine check_nan_values_velocity


  subroutine check_nan_values_gravity(grav_field)
    use mod_types, only: gravity_type

    type(gravity_type), intent(in)  :: grav_field

    call stop_if_nan(grav_field % grav, "g")
  end subroutine check_nan_values_gravity

end module mod_check_values
