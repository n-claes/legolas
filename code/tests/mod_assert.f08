module mod_assert
  implicit none

  double precision, parameter       :: tol = 1.0d-5
  integer                           :: test_success, test_fail, test_total

contains

  subroutine assert_init()
    test_success = 0
    test_fail = 0
    test_total = 0
  end subroutine assert_init


  subroutine assert_equal(a, b)
    double precision, intent(in) :: a, b

    if (abs(a - b) < tol) then
      call increment_success
    else
      call increment_failure
      print*,"FAIL: numbers not equal within", tol,". Actual value", &
             abs(a-b)
    end if
  end subroutine assert_equal

  subroutine assert_larger(a, b)
    double precision, intent(in) :: a, b

    if (a > b - tol) then
      call increment_success
    else
      call increment_failure
      print*,"FAIL: a not larger than b witin", tol,". Actual values a=",a, &
             "b=",b
    end if
  end subroutine assert_larger

  subroutine assert_smaller(a, b)
    double precision, intent(in) :: a, b

    if (a < b + tol) then
      call increment_success
    else
      call increment_failure
      print*,"FAIL: number not smaller within", tol,"."
    end if
  end subroutine assert_smaller


  subroutine increment_success()
    test_success = test_success + 1
    test_total = test_total + 1
    print*,"--> OK"
  end subroutine increment_success

  subroutine increment_failure()
    test_fail = test_fail + 1
    test_total = test_total + 1
    print*,"--> FAIL"
  end subroutine increment_failure

  subroutine get_test_results()
    print*,
    print*,
    print*,"===================="
    print*,"TESTS TOTAL    : ", test_total
    print*,"===================="
    print*,"TESTS SUCCEEDED: ", test_success
    print*,"TESTS FAILED   : ", test_fail

  end subroutine get_test_results


end module mod_assert
