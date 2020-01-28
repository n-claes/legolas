!
! SUBMODULE: smod_test_hydro
! 
! DESCRIPTION:
!> This submodule tests the hydrodynamic case (no magnetic fields).
submodule (mod_equilibrium) smod_test_hydro
  implicit none
  
contains
  
  module subroutine hydro_test_eq()
    write(*, *) "This test displays the adiabatic homogeneous case "
    write(*, *) "in the hydrodynamics limit (B=0)."
    
    call adiabatic_homo_eq()
    
    B_field % B02 = 0.0d0
    B_field % B03 = 0.0d0
    B_field % B0 = 0.0d0
  end subroutine hydro_test_eq
  
end submodule smod_test_hydro
