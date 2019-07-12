program tests_main
  use testmod_core_tests
  use testmod_homogeneous
  implicit none

  ! Core tests
  call init_core()
  call run_core_tests()
  call finish_core()

  ! Homogeneous adiabatic medium
  call init_homogeneous()
  call run_homogeneous_test()
  call finish_homogeneous()
  call execute_command_line("python python/test_homogeneous.py")



end program tests_main
