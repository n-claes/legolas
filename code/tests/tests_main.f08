program tests_main
  use testmod_core_tests
  use testmod_homo_adiabatic
  use testmod_homo_gravity
  implicit none

  ! Core tests
  call init_core()
  call run_core_tests()
  call finish_core()

  ! Homogeneous adiabatic medium
  call init_homo_adiabatic_test()
  call run_homo_adiabatic_test()
  call finish_homo_adiabatic_test()

  ! Homogeneous medium with gravity
  call init_homo_gravity_test()
  call run_homo_gravity_test()
  call finish_homo_gravity_test()


  call execute_command_line("python python/test_homogeneous.py")



end program tests_main
