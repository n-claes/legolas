!
! SUBMODULE: smod_test_beta0
! 
! DESCRIPTION: 
!> This submodule tests the case for beta = 0.
submodule (mod_equilibrium) smod_test_beta0
  implicit none
  
contains
  
  module subroutine beta0_test_eq()
    
    call adiabatic_homo_eq()
    
    write(*, *) "This test displays the adiabatic homogeneous case in the beta=0 limit."
    
    T_field % T0 = 0.0d0 
    T_field % d_T0_dr = 0.0d0
    
  end subroutine beta0_test_eq

end submodule smod_test_beta0
