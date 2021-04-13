! =============================================================================
!> This module performs sanity checks on given values. Methods are implemented to check for
!! small, NaN or negative values in arrays. We also have single-value checks against
!! zero or NaN, along with a double-precision equality check.
module mod_check_values
  use, intrinsic  :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  use mod_global_variables, only: dp, dp_LIMIT, gauss_gridpts
  use mod_logging, only: log_message, int_fmt, exp_fmt, char_log
  implicit none

  private

  !> Interface to check small values in a real or complex array/matrix.
  interface check_small_values
    module procedure small_values_real_array
    module procedure small_values_complex_array
    module procedure small_values_real_matrix
    module procedure small_values_complex_matrix
  end interface check_small_values

  !> Interface to check if a given matrix is square
  interface matrix_is_square
    module procedure matrix_real_is_square
    module procedure matrix_complex_is_square
  end interface matrix_is_square

  !> Interface to check NaN values in the equilibrium types.
  interface check_nan_values
    module procedure check_nan_values_density
    module procedure check_nan_values_temperature
    module procedure check_nan_values_bfield
    module procedure check_nan_values_velocity
    module procedure check_nan_values_gravity
  end interface check_nan_values

  !> Interface to check if a value is equal to another value.
  interface value_is_equal
    module procedure real_value_is_equal
    module procedure real_array_is_equal
    module procedure real_array_is_equal_to_value
  end interface value_is_equal

  !> Interface to check if a real or complex value is zero.
  interface value_is_zero
    module procedure real_is_zero
    module procedure complex_is_zero
  end interface value_is_zero

  !> Interface to check if a real or complex value is infinity.
  interface value_is_inf
    module procedure real_is_inf
    module procedure complex_is_inf
  end interface value_is_inf

  public :: check_small_values
  public :: check_negative_array
  public :: check_nan_values
  public :: matrix_is_square
  public :: value_is_zero
  public :: value_is_equal
  public :: value_is_nan
  public :: value_is_inf
  public :: check_coax

contains


  !> Double-precision check for zero.
  !! Checks if a given double-precision real value equals zero
  !! within a <tt>DP_LIMIT</tt> margin.
  !! Returns <tt>True</tt> if value equals zero, <tt>False</tt> otherwise.
  function real_is_zero(value) result(is_zero)
    !> the real value to check
    real(dp), intent(in)  :: value
    logical :: is_zero

    if (abs(value - 0.0d0) > dp_LIMIT) then
      is_zero = .false.
    else
      is_zero = .true.
    end if
  end function real_is_zero


  !> Double-precision complex check for zero.
  !! Checks if a given double-precision complex value equals zero
  !! within a <tt>DP_LIMIT</tt> margin.
  !! Returns <tt>True</tt> if value equals zero, <tt>False</tt> otherwise.
  !! @note  Both the real and imaginary parts have to
  !!        equal zero for this function to return <tt>True</tt>.
  function complex_is_zero(value) result(is_zero)
    !> the complex value to check
    complex(dp), intent(in) :: value
    logical :: is_zero

    if (abs(real(value) - 0.0d0) > dp_LIMIT .or. abs(aimag(value) - 0.0d0) > dp_LIMIT) then
      is_zero = .false.
    else
      is_zero = .true.
    end if
  end function complex_is_zero


  !> Equality check between variables.
  !! Checks if two real double-precision values are equal to
  !! each other within <tt>DP_LIMIT</tt>.
  !! Returns <tt>True</tt> if _value == value_base_, <tt>False</tt> otherwise.
  function real_value_is_equal(value, value_base, tol)  result(is_equal)
    !> the first real value
    real(dp), intent(in)  :: value
    !> the base value to check against
    real(dp), intent(in)  :: value_base
    !> optional tolerance used for comparison
    real(dp), intent(in), optional  :: tol
    logical   :: is_equal
    real(dp)  :: tolerance

    if (present(tol)) then
      tolerance = tol
    else
      tolerance = dp_LIMIT
    end if

    if (abs(value - value_base) > tolerance) then
      is_equal = .false.
    else
      is_equal = .true.
    end if
  end function real_value_is_equal


  !> Checks if all values of a given array are equal to a given value.
  function real_array_is_equal_to_value(array, value_base, tol)  result(is_equal)
    !> the first array of values
    real(dp), intent(in)  :: array(:)
    !> the base value to check against
    real(dp), intent(in)  :: value_base
    !> optional tolerance used for comparison
    real(dp), intent(in), optional  :: tol
    logical   :: is_equal
    real(dp)  :: tolerance

    if (present(tol)) then
      tolerance = tol
    else
      tolerance = dp_LIMIT
    end if

    if (all(abs(array - value_base) < dp_LIMIT)) then
      is_equal = .true.
    else
      is_equal = .false.
    end if
  end function real_array_is_equal_to_value


  !> Checks if two real arrays are equal.
  function real_array_is_equal(array, array_base, tol)  result(is_equal)
    !> the first array of values
    real(dp), intent(in)  :: array(:)
    !> the base values to check against
    real(dp), intent(in)  :: array_base(:)
    !> optional tolerance used for comparison
    real(dp), intent(in), optional  :: tol
    logical   :: is_equal
    real(dp)  :: tolerance

    if (present(tol)) then
      tolerance = tol
    else
      tolerance = dp_LIMIT
    end if

    if (all(abs(array - array_base) < dp_LIMIT)) then
      is_equal = .true.
    else
      is_equal = .false.
    end if
  end function real_array_is_equal


  !> NaN check for real values.
  !! Checks if a real double-precision value is equal to NaN.
  !! Returns <tt>True</tt> if value equals NaN, <tt>False</tt> otherwise.
  function value_is_nan(value)  result(is_nan)
    !> the real value to check
    real(dp), intent(in)  :: value
    logical :: is_nan

    if (ieee_is_nan(value)) then
      is_nan = .true.
    else
      is_nan = .false.
    end if
  end function value_is_nan


  !> Infinity check for real values.
  !! Returns <tt>True</tt> if value equals Infinity, <tt>False</tt> otherwise
  function real_is_inf(value) result(is_inf)
    !> the real value to check
    real(dp), intent(in)  :: value
    logical :: is_inf

    if (ieee_is_finite(value)) then
      is_inf = .false.
    else
      is_inf = .true.
    end if
  end function real_is_inf


  !> Infinity check for complex values.
  !! Returns <tt>True</tt> if value equals Infinity, <tt>False</tt> otherwise
  function complex_is_inf(value) result(is_inf)
    !> the real value to check
    complex(dp), intent(in)  :: value
    logical :: is_inf

    if (ieee_is_finite(real(value)) .and. ieee_is_finite(aimag(value))) then
      is_inf = .false.
    else
      is_inf = .true.
    end if
  end function complex_is_inf


  !> Small value checks on a real array.
  !! Checks a real double-precision array for small values.
  !! Values that are smaller than the specified tolerance <tt>tol</tt> are set to zero.
  !! If <tt>tol</tt> is not present, <tt>DP_LIMIT</tt> is used as tolerance.
  subroutine small_values_real_array(array, tol)
    !> the real array to check
    real(dp), intent(inout)         :: array(:)
    !> the tolerance to check against
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


  !> Small value checks on a complex array.
  !! Checks a complex double-precision array for small values.
  !! Values that are smaller than the specified tolerance <tt>tol</tt> are set to zero.
  !! The real and imaginary parts are checked (and set) separately.
  !! If <tt>tol</tt> is not present, <tt>DP_LIMIT</tt> is used as tolerance.
  subroutine small_values_complex_array(array, tol)
    !> the complex array to check
    complex(dp), intent(inout)      :: array(:)
    !> the tolerance to check against
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


  !> Small value checks on a real matrix.
  !! Checks a real double-precision matrix for small values.
  !! Values that are smaller than the specified tolerance <tt>tol</tt> are set to zero.
  !! If <tt>tol</tt> is not present, <tt>DP_LIMIT</tt> is used as tolerance.
  subroutine small_values_real_matrix(matrix, tol)
    !> the real matrix to check
    real(dp), intent(inout)         :: matrix(:, :)
    !> the tolerance to check against
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


  !> Small value checks on a complex matrix.
  !! Checks a complex double-precision matrix for small values.
  !! Values that are smaller than the specified tolerance <tt>tol</tt> are set to zero.
  !! The real and imaginary parts are checked (and set) separately.
  !! If <tt>tol</tt> is not present, <tt>DP_LIMIT</tt> is used as tolerance.
  subroutine small_values_complex_matrix(matrix, tol)
    !> the complex matrix to check
    complex(dp), intent(inout)      :: matrix(:, :)
    !> the tolerance to check against
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


  !> Checks if a given real matrix is square.
  !! Returns <tt>.true.</tt> if dimensions are equal, <tt>.false.</tt> otherwise.
  function matrix_real_is_square(matrix)  result(is_square)
    !> the real matrix to check
    real(dp), intent(in)  :: matrix(:, :)
    logical :: is_square

    if (size(matrix, dim=1) == size(matrix, dim=2)) then
      is_square = .true.
    else
      is_square = .false.
    end if
  end function matrix_real_is_square


  !> Checks if a given complex matrix is square.
  !! Returns <tt>.true.</tt> if dimensions are equal, <tt>.false.</tt> otherwise.
  function matrix_complex_is_square(matrix)  result(is_square)
    !> the real matrix to check
    complex(dp), intent(in)  :: matrix(:, :)
    logical :: is_square

    if (size(matrix, dim=1) == size(matrix, dim=2)) then
      is_square = .true.
    else
      is_square = .false.
    end if
  end function matrix_complex_is_square


  !> Negative value checks on a real array.
  !! Checks a real double-precision array for negative values.
  !! This is meant to check density, temperature, etc.
  !! @warning Throws an error if <tt>array</tt> contains a negative value.
  subroutine check_negative_array(array, variable_name)
    !> the real array to check
    real(dp), intent(in)          :: array(:)
    !> the name to include in the error message
    character(len=*), intent(in)  :: variable_name
    integer                       :: i

    do i = 1, size(array)
      if (array(i) < 0.0d0) then
        call log_message("negative value encountered in " // trim(variable_name), level='error')
        exit
      end if
    end do
  end subroutine check_negative_array


  !> Stops program if NaN is encountered in a real double-precision array.
  !! @warning If <tt>array</tt> contains a NaN value an error is thrown.
  subroutine stop_if_nan(array, array_name)
    !> the real array to check
    real(dp), intent(in)  :: array(gauss_gridpts)
    !> the name to include in the error message
    character(len=*), intent(in)  :: array_name
    integer :: i

    do i = 1, gauss_gridpts
      if (ieee_is_nan(array(i))) then
        call log_message("NaN encountered in " // trim(array_name), level='error')
      end if
    end do
  end subroutine stop_if_nan


  !> Checks if one of the array attributes of the density type contains NaN.
  !! @warning If NaN is encountered an error is thrown.
  subroutine check_nan_values_density(rho_field)
    use mod_types, only: density_type

    !> the type containing the density attributes
    type(density_type), intent(in)  :: rho_field

    call stop_if_nan(rho_field % rho0, "rho0")
    call stop_if_nan(rho_field % d_rho0_dr, "drho0")
  end subroutine check_nan_values_density


  !> Checks if one of the array attributes of the temperature type contains NaN.
  !! @warning If NaN is encountered an error is thrown.
  subroutine check_nan_values_temperature(T_field)
    use mod_types, only: temperature_type

    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field

    call stop_if_nan(T_field % T0, "T0")
    call stop_if_nan(T_field % d_T0_dr, "dT0")
  end subroutine check_nan_values_temperature


  !> Checks if one of the array attributes of the magnetic field type contains NaN.
  !! @warning If NaN is encountered an error is thrown.
  subroutine check_nan_values_bfield(B_field)
    use mod_types, only: bfield_type

    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in) :: B_field

    call stop_if_nan(B_field % B02, "B02")
    call stop_if_nan(B_field % B03, "B03")
    call stop_if_nan(B_field % B0, "B0")
    call stop_if_nan(B_field % d_B02_dr, "dB02")
    call stop_if_nan(B_field % d_B03_dr, "dB03")
  end subroutine check_nan_values_bfield


  !> Checks if one of the array attributes of the velocity type contains NaN.
  !! @warning If NaN is encountered an error is thrown.
  subroutine check_nan_values_velocity(v_field)
    use mod_types, only: velocity_type

    !> the type containing the velocity attributes
    type(velocity_type), intent(in) :: v_field

    call stop_if_nan(v_field % v02, "v02")
    call stop_if_nan(v_field % v03, "v03")
    call stop_if_nan(v_field % d_v02_dr, "dv02")
    call stop_if_nan(v_field % d_v03_dr, "dv03")
  end subroutine check_nan_values_velocity


  !> Checks if one of the array attributes of the gravity type contains NaN.
  !! @warning If NaN is encountered an error is thrown.
  subroutine check_nan_values_gravity(grav_field)
    use mod_types, only: gravity_type

    !> the type containing the gravity attributes
    type(gravity_type), intent(in)  :: grav_field

    call stop_if_nan(grav_field % grav, "g")
  end subroutine check_nan_values_gravity


  !> Checks if x_start is non-zero if an inner wall boundary is included for
  !! a cylindrical geometry
  subroutine check_coax()
    use mod_global_variables, only: dp_LIMIT, x_start, coaxial

    if (coaxial .and. x_start < dp_LIMIT) then
      call log_message("x_start must be larger than zero to introduce an inner wall boundary", level='error')
    end if
  end subroutine check_coax

end module mod_check_values
