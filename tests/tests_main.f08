program tests_main
  use testmod_core_tests
  implicit none

  ! Core tests
  call init_core()
  call run_core_tests()
  call finish_core()

end program tests_main
