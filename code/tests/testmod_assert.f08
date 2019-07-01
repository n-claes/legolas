module testmod_assert
  use, intrinsic :: iso_fortran_env
  implicit none

  real(real64), parameter   :: tol = 1.0d-8

contains

  !> Checks if boolean is true
  subroutine assert_true(bool1, bool)
    logical, intent(in)   :: bool1
    logical, intent(out)  :: bool

    if (bool1) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_true

  !> Checks if boolean is false
  subroutine assert_false(bool1, bool)
    logical, intent(in)   :: bool1
    logical, intent(out)  :: bool

    if (bool1) then
      bool = .false.
      return
    else
      bool = .true.
      return
    end if
  end subroutine assert_false

  !> Checks if int1 = int2
  subroutine assert_int_equal(int1, int2, bool)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: bool

    if (int1 == int2) then
      bool = .true.
    else
      bool = .false.
    end if
  end subroutine assert_int_equal

  !> Checks if int1 > int2
  subroutine assert_int_larger(int1, int2, bool)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: bool

    if (int1 > int2) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_int_larger

  !> Checks if int1 < int2
  subroutine assert_int_smaller(int1, int2, bool)
    integer, intent(in)   :: int1, int2
    logical, intent(out)  :: bool

    if (int1 < int2) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_int_smaller

  !> Checks if dble1 = dble2
  subroutine assert_real_equal(dble1, dble2, bool)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: bool

    if (abs(dble1 - dble2) < tol) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_real_equal

  !> Checks if dble1 > dble2
  subroutine assert_real_larger(dble1, dble2, bool)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: bool

    if (dble1 - tol > dble2) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_real_larger

  !> Checks if dble1 < dble2
  subroutine assert_real_smaller(dble1, dble2, bool)
    real(real64), intent(in)  :: dble1, dble2
    logical, intent(out)      :: bool

    if (dble1 + tol < dble2) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_real_smaller

  !> Checks if z1 = z2
  subroutine assert_complex_equal(z1, z2, bool)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: bool

    if (abs(real(z1) - real(z2)) < tol .and. &
        abs(aimag(z1) - aimag(z2)) < tol) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_complex_equal


  !> Checks if z1 > z2
  subroutine assert_complex_larger(z1, z2, bool)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: bool

    if (real(z1) - tol < real(z2) .and. &
        aimag(z1) - tol > aimag(z2)) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_complex_larger


  !> Checks if z1 < z2
  subroutine assert_complex_smaller(z1, z2, bool)
    complex(real64), intent(in) :: z1, z2
    logical, intent(out)        :: bool

    if (real(z1) + tol < real(z2) .and. &
        aimag(z1) + tol < aimag(z2)) then
      bool = .true.
      return
    else
      bool = .false.
      return
    end if
  end subroutine assert_complex_smaller















end module testmod_assert
