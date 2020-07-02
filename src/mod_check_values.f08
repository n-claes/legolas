! =============================================================================
!> @brief   Module for sanity checks on values.
!! @details Performs sanity checks on given values. Methods are implemented to check for
!!          small, NaN or negative values in arrays. We also have single-value checks against
!!          zero or NaN, along with a double-precision equality check.
module mod_check_values
  use, intrinsic  :: ieee_arithmetic, only: ieee_is_nan
  use mod_global_variables, only: dp, dp_LIMIT, gauss_gridpts
  use mod_logging, only: log_message, int_fmt, exp_fmt, char_log
  implicit none

  private

  !> Interface to check small values in a real/complex array/matrix.
  interface check_small_values
    module procedure small_values_real_array
    module procedure small_values_complex_array
    module procedure small_values_real_matrix
    module procedure small_values_complex_matrix
  end interface check_small_values

  !> Interface to check NaN values in the equilibrium types.
  interface check_nan_values
    module procedure check_nan_values_density
    module procedure check_nan_values_temperature
    module procedure check_nan_values_bfield
    module procedure check_nan_values_velocity
    module procedure check_nan_values_gravity
  end interface check_nan_values

  !> Interface to check if a real/complex value is zero.
  interface value_is_zero
    module procedure real_is_zero
    module procedure complex_is_zero
  end interface value_is_zero

  public :: check_small_values
  public :: check_negative_array
  public :: check_nan_values
  public :: value_is_zero
  public :: value_is_equal
  public :: value_is_nan

contains


  !> @brief   Double-precision check for zero.
  !! @details Checks if a given double-precision real value equals zero
  !!          within a \p 1e-12 (value for \p DP_LIMIT).
  !! @param[in] value   the real value to check
  !! @return  \p True if \p value is zero, \p False otherwise
  function real_is_zero(value) result(is_zero)
    real(dp), intent(in)  :: value
    logical :: is_zero

    if (abs(value - 0.0d0) > dp_LIMIT) then
      is_zero = .false.
    else
      is_zero = .true.
    end if
  end function real_is_zero


  !> @brief   Double-precision complex check for zero.
  !! @details Checks if a given double-precision complex value equals zero
  !!          within \p 1e-12 (value for \p DP_LIMIT). For this to be true, both the real
  !!          and imaginary parts have to be zero.
  !! @param[in] value   the complex value to check
  !! @returns \p True if real and imaginary parts are zero, \p False otherwise
  function complex_is_zero(value) result(is_zero)
    complex(dp), intent(in) :: value
    logical :: is_zero

    if (abs(real(value) - 0.0d0) > dp_LIMIT .or. abs(aimag(value) - 0.0d0) > dp_LIMIT) then
      is_zero = .false.
    else
      is_zero = .true.
    end if
  end function complex_is_zero


  !> @brief   Equality check between variables.
  !! @details Checks if two real double-precision values are equal to
  !!          each other within \p 1e-12 (value for \p DP_LIMIT).
  !! @param[in] value   the first real value
  !! @param[in] value_base  the base value against which to check \p value
  !! @returns \p True if \p value == \p value_base, \p False otherwise
  function value_is_equal(value, value_base)  result(is_equal)
    real(dp), intent(in)  :: value, value_base
    logical :: is_equal

    if (abs(value - value_base) > dp_LIMIT) then
      is_equal = .false.
    else
      is_equal = .true.
    end if
  end function value_is_equal


  !> @brief   NaN check for real values.
  !! @details Checks if a real double-precision value is equal to NaN.
  !! @param[in] value   the real value to check
  !! @returns \p True if \p value equals \p NaN, \p False otherwise
  function value_is_nan(value)  result(is_nan)
    real(dp), intent(in)  :: value
    logical :: is_nan

    if (ieee_is_nan(value)) then
      is_nan = .true.
    else
      is_nan = .false.
    end if
  end function value_is_nan


  !> @brief   Small value checks on a real array.
  !! @details Checks a real double-precision array for small values.
  !!          Values that are smaller than the specified tolerance \p tol are set to zero.
  !!          If \p tol is not present, \p DP_LIMIT is used as tolerance.
  !! @param[in] array   the real array to check
  !! @param[in] tol     [optional] the tolerance to check against
  subroutine small_values_real_array(array, tol)
    real(dp), intent(inout)         :: array(:)
    real(dp), intent(in), optional  :: tol
    real(dp)    :: limit
    integer     :: i

    if (present(tol)) then
      limit = tol
    else
      limit = dp_LIMIT
    end if

    do i = 1, size(array)
      if (abs(array(i)) < limit) then
        array(i) = 0.0d0
      end if
    end do
  end subroutine small_values_real_array


  !> @brief   Small value checks on a complex array.
  !! @details Checks a complex double-precision array for small values.
  !!          Values that are smaller than the specified tolerance \p tol are set to zero.
  !!          The real and imaginary parts are checked (and set) separately.
  !!          If \p tol is not present, \p DP_LIMIT is used as tolerance.
  !! @param[in] array   the complex array to check
  !! @param[in] tol     [optional] the tolerance to check against
  subroutine small_values_complex_array(array, tol)
    complex(dp), intent(inout)      :: array(:)
    real(dp), intent(in), optional  :: tol
    real(dp)    :: a_real, a_imag
    real(dp)    :: limit
    integer     :: i

    if (present(tol)) then
      limit = tol
    else
      limit = dp_LIMIT
    end if

    do i = 1, size(array)
      a_real = real(array(i))
      a_imag = aimag(array(i))

      if (abs(a_real) < limit) then
        a_real = 0.0d0
      end if
      if (abs(a_imag) < limit) then
        a_imag = 0.0d0
      end if

      array(i) = cmplx(a_real, a_imag, kind=dp)
    end do
  end subroutine small_values_complex_array


  !> @brief   Small value checks on a real matrix.
  !! @details Checks a real double-precision matrix for small values.
  !!          Values that are smaller than the specified tolerance \p tol are set to zero.
  !!          If \p tol is not present, \p DP_LIMIT is used as tolerance.
  !! @param[in] matrix  the real matrix to check
  !! @param[in] tol     [optional] the tolerance to check against
  subroutine small_values_real_matrix(matrix, tol)
    real(dp), intent(inout)         :: matrix(:, :)
    real(dp), intent(in), optional  :: tol
    real(dp)    :: limit
    integer     :: i, j

    if (present(tol)) then
      limit = tol
    else
      limit = dp_LIMIT
    end if

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        if (abs(matrix(i, j)) < limit) then
          matrix(i, j) = 0.0d0
        end if
      end do
    end do
  end subroutine small_values_real_matrix


  !> @brief   Small value checks on a complex matrix.
  !! @details Checks a complex double-precision matrix for small values.
  !!          Values that are smaller than the specified tolerance \p tol are set to zero.
  !!          The real and imaginary parts are checked (and set) separately.
  !!          If \p tol is not present, \p DP_LIMIT is used as tolerance.
  !! @param[in] matrix  the complex matrix to check
  !! @param[in] tol     [optional] the tolerance to check against
  subroutine small_values_complex_matrix(matrix, tol)
    complex(dp), intent(inout)      :: matrix(:, :)
    real(dp), intent(in), optional  :: tol
    real(dp)    :: a_real, a_imag
    real(dp)    :: limit
    integer     :: i, j

    if (present(tol)) then
      limit = tol
    else
      limit = dp_LIMIT
    end if

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        a_real = real(matrix(i, j))
        a_imag = aimag(matrix(i, j))

        if (abs(a_real) < limit) then
          a_real = 0.0d0
        end if
        if (abs(a_imag) < limit) then
          a_imag = 0.0d0
        end if

        matrix(i, j) = cmplx(a_real, a_imag, kind=dp)
      end do
    end do
  end subroutine small_values_complex_matrix


  !> @brief   Negative value checks on a real array.
  !! @details Checks a real double-precision array for negative values.
  !!          If a negative value is encountered an error is thrown.
  !!          This is meant to check density, temperature, etc.
  !! @exception Error   If \p array contains a negative value.
  !! @param[in] array   the real array to check
  !! @param[in] variable_name   the name to include in the error message
  subroutine check_negative_array(array, variable_name)
    real(dp), intent(in)          :: array(:)
    character(len=*), intent(in)  :: variable_name
    integer                       :: i

    do i = 1, size(array)
      if (array(i) < 0.0d0) then
        call log_message("negative value encountered in " // trim(variable_name), level='error')
      end if
    end do
  end subroutine check_negative_array


  !> @brief   Throws an error if NaN is encountered in an array.
  !! @details Checks a real double-precision array for NaN values,
  !!          and throws an error if one is encountered.
  !! @exception Error If \p array contains a \p NaN value.
  !! @param[in] array   the real array to check
  !! @param[in] array_name  the name to include in the error message
  subroutine stop_if_nan(array, array_name)
    real(dp), intent(in)  :: array(gauss_gridpts)
    character(len=*), intent(in)  :: array_name
    integer :: i

    do i = 1, gauss_gridpts
      if (ieee_is_nan(array(i))) then
        call log_message("NaN encountered in " // trim(array_name), level='error')
      end if
    end do
  end subroutine stop_if_nan


  !> @brief   NaN checks for the density.
  !! @details Checks if one of the array attributes of the density type
  !!          contains a NaN value.
  !! @exception Error   If a \p NaN value is encountered.
  !! @param[in] rho_field   the type containing the density attributes
  subroutine check_nan_values_density(rho_field)
    use mod_types, only: density_type

    type(density_type), intent(in)  :: rho_field

    call stop_if_nan(rho_field % rho0, "rho0")
    call stop_if_nan(rho_field % d_rho0_dr, "drho0")
  end subroutine check_nan_values_density


  !> @brief   NaN checks for the temperature.
  !! @details Checks if one of the array attributes of
  !!          the temperature type contains a NaN value.
  !! @exception Error   If a \p NaN value is encountered.
  !! @param[in] T_field   the type containing the temperature attributes
  subroutine check_nan_values_temperature(T_field)
    use mod_types, only: temperature_type

    type(temperature_type), intent(in)  :: T_field

    call stop_if_nan(T_field % T0, "T0")
    call stop_if_nan(T_field % d_T0_dr, "dT0")
  end subroutine check_nan_values_temperature


  !> @brief   NaN checks for the magnetic field.
  !! @details Checks if one of the array attributes of
  !!          the magnetic field type contains a NaN value.
  !! @exception Error   If a \p NaN value is encountered.
  !! @param[in] B_field   the type containing the magnetic field attributes
  subroutine check_nan_values_bfield(B_field)
    use mod_types, only: bfield_type

    type(bfield_type), intent(in) :: B_field

    call stop_if_nan(B_field % B02, "B02")
    call stop_if_nan(B_field % B03, "B03")
    call stop_if_nan(B_field % B0, "B0")
    call stop_if_nan(B_field % d_B02_dr, "dB02")
    call stop_if_nan(B_field % d_B03_dr, "dB03")
  end subroutine check_nan_values_bfield


  !> @brief   NaN checks for the velocity.
  !! @details Checks if one of the array attributes of
  !!          the velocity type contains a NaN value.
  !! @exception Error   If a \p NaN value is encountered.
  !! @param[in] v_field   the type containing the velocity attributes
  subroutine check_nan_values_velocity(v_field)
    use mod_types, only: velocity_type

    type(velocity_type), intent(in) :: v_field

    call stop_if_nan(v_field % v02, "v02")
    call stop_if_nan(v_field % v03, "v03")
    call stop_if_nan(v_field % d_v02_dr, "dv02")
    call stop_if_nan(v_field % d_v03_dr, "dv03")
  end subroutine check_nan_values_velocity


  !> @brief   NaN checks for gravity.
  !! @details Checks if one of the array attributes of
  !!          the gravity type contains a NaN value.
  !! @exception Error   If a \p NaN value is encountered.
  !! @param[in] grav_field   the type containing the gravity attributes
  subroutine check_nan_values_gravity(grav_field)
    use mod_types, only: gravity_type

    type(gravity_type), intent(in)  :: grav_field

    call stop_if_nan(grav_field % grav, "g")
  end subroutine check_nan_values_gravity

end module mod_check_values
