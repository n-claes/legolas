module testmod_assert
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  implicit none

  private

  real(real64), parameter   :: tol = 1.0d-12

  interface assert_equal
    module procedure assert_int_equal
    module procedure assert_real_equal
    module procedure assert_complex_equal
  end interface assert_equal

  interface assert_smaller
    module procedure assert_int_smaller
    module procedure assert_real_smaller
    module procedure assert_complex_smaller
  end interface assert_smaller

  interface assert_larger
    module procedure assert_int_larger
    module procedure assert_real_larger
    module procedure assert_complex_larger
  end interface assert_larger

  interface assert_is_finite
    module procedure assert_real_is_finite
    module procedure assert_complex_is_finite
  end interface assert_is_finite

  interface assert_is_no_nan
    module procedure assert_real_is_no_nan
    module procedure assert_complex_is_no_nan
  end interface assert_is_no_nan

  public :: assert_true, assert_false
  public :: assert_equal
  public :: assert_smaller
  public :: assert_larger
  public :: assert_is_finite
  public :: assert_is_no_nan

contains

  !> Checks if boolean is true
  subroutine assert_true(bool_in, passed)
    logical, intent(in)   :: bool_in
    logical, intent(out)  :: passed

    if (bool_in) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_true

  !> Checks if boolean is false
  subroutine assert_false(bool_in, passed)
    logical, intent(in)   :: bool_in
    logical, intent(out)  :: passed

    if (bool_in) then
      passed = .false.
    else
      passed = .true.
    end if
  end subroutine assert_false





  !! Assert equal statements

  !> Integer assertion, checks if two integers are equal.
  !! @param[in] int1      First integer
  !! @param[in] int2      Second integer
  !! @param[out] passed   True if int1 = int2, false otherwise
  subroutine assert_int_equal(int1, int2, passed)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: passed

    if (int1 == int2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_int_equal

  !> Double precision real assertion, checks if two reals are equal within
  !! a given tolerance.
  !! @param[in] dble1     First double precision real
  !! @param[in] dble2     Second double precision real
  !! @param[out] passed   True if dble1 = dble2, false otherwise
  subroutine assert_real_equal(dble1, dble2, passed)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: passed

    if (abs(dble1 - dble2) < tol) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_real_equal

  !> Double precision complex assertion, checks if two complex numbers are
  !! equal within a given tolerance. Both the real and imaginary parts must
  !! be equal.
  !! @param[in] z1      First double precision complex
  !! @param[in] z2      Second double precision complex
  !! @param[out] passed True if z1 = z2, false otherwise
  subroutine assert_complex_equal(z1, z2, passed)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: passed

    if (abs(real(z1) - real(z2)) < tol .and. &
        abs(aimag(z1) - aimag(z2)) < tol) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_complex_equal





  !! Assert smaller statements

  !> Integer assertion, checks if one integer is smaller than the other.
  !! @param[in]   int1    First integer
  !! @param[in]   int2    Second integer
  !! @param[out]  passed  True if int1 < int2, false otherwise
  subroutine assert_int_smaller(int1, int2, passed)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: passed

    if (int1 < int2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_int_smaller

  !> Double precision real assertion. Checks if one real is smaller than the
  !! other one within a given tolerance.
  !! @param[in]  dble1  First double precision real
  !! @param[in]  dble2  Second double precision real
  !! @param[out] passed true if dble1 < dble2, false otherwise
  subroutine assert_real_smaller(dble1, dble2, passed)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: passed

    if (dble1 + tol < dble2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_real_smaller

  !> Double precision complex assertion. Checks if one complex number is
  !! smaller than the other one within a given tolerance. Both the real and
  !! imaginary parts must be smaller in order to be true.
  !! @param[in]  z1     First double precision complex
  !! @param[in]  z2     Second double precision complex
  !! @param[out] passed True if z1 < z2 (both real and imaginary parts),
  !!                    false otherwise
  subroutine assert_complex_smaller(z1, z2, passed)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: passed

    if (real(z1) + tol < real(z2) .and. &
        aimag(z1) + tol < aimag(z2)) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_complex_smaller




  !! Assert larger statements

  !> Integer assertion, checks if one integer is larger than the other.
  !! @param[in]   int1    First integer
  !! @param[in]   int2    Second integer
  !! @param[out]  passed  True if int1 > int2, false otherwise
  subroutine assert_int_larger(int1, int2, passed)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: passed

    if (int1 > int2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_int_larger

  !> Double precision real assertion. Checks if one real is larger than the
  !! other one within a given tolerance.
  !! @param[in]  dble1  First double precision real
  !! @param[in]  dble2  Second double precision real
  !! @param[out] passed true if dble1 > dble2, false otherwise
  subroutine assert_real_larger(dble1, dble2, passed)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: passed

    if (dble1 - tol > dble2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_real_larger

  !> Double precision complex assertion. Checks if one complex number is
  !! larger than the other one within a given tolerance. Both the real and
  !! imaginary parts must be larger in order to be true.
  !! @param[in]  z1     First double precision complex
  !! @param[in]  z2     Second double precision complex
  !! @param[out] passed True if z1 > z2 (both real and imaginary parts),
  !!                    false otherwise
  subroutine assert_complex_larger(z1, z2, passed)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: passed

    if (real(z1) - tol > real(z2) .and. &
        aimag(z1) - tol > aimag(z2)) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_complex_larger




  !! Finite and NaN assertions

  !> Double precision real assertion, checks if real is infinite.
  !! @param[in] dble1   Double precision real to check
  !! @param[out] passed True if dble1 is finite, false if dlbe1 = inf
  subroutine assert_real_is_finite(dble1, passed)
    real(real64), intent(in)  :: dble1
    logical, intent(out)      :: passed

    passed = ieee_is_finite(dble1)
  end subroutine assert_real_is_finite

  !> Double precision complex assertion, checks if complex number is infinite.
  !! If either the real or imaginary part is infinite, assertion will fail.
  !! @param[in]  z1      Double precision complex to check
  !! @param[out] passed  True if z1 is finite, false if z1 = inf
  subroutine assert_complex_is_finite(z1, passed)
    complex(real64), intent(in) :: z1
    logical, intent(out)        :: passed
    logical                     :: bool1, bool2

    bool1 = ieee_is_finite(real(z1))
    bool2 = ieee_is_finite(aimag(z1))
    if (bool1 .and. bool2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_complex_is_finite


  !> Double precision real assertion, checks if real is NaN.
  !! @param[in]  dble1    Integer to check
  !! @param[out] passed   True if dble1 is a number, false if dble1 = NaN
  subroutine assert_real_is_no_nan(dble1, passed)
    real(real64), intent(in)  :: dble1
    logical, intent(out)      :: passed

    passed = .not. ieee_is_nan(dble1)
  end subroutine assert_real_is_no_nan


  !> Double precision complex assertion, checks if complex number is NaN.
  !! If either the real or imaginary part is NaN, assertion will fail.
  !! @param[in]  z1      Double precision complex to check
  !! @param[out] passed  True if z1 is a number, false if z1 = NaN
  subroutine assert_complex_is_no_nan(z1, passed)
    complex(real64), intent(in) :: z1
    logical, intent(out)        :: passed

    logical     :: bool1, bool2

    bool1 = .not. ieee_is_nan(real(z1))
    bool2 = .not. ieee_is_nan(aimag(z1))
    if (bool1 .and. bool2) then
      passed = .true.
    else
      passed = .false.
    end if
  end subroutine assert_complex_is_no_nan

end module testmod_assert
