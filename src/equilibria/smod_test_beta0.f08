!
! SUBMODULE: smod_test_beta0
!
! DESCRIPTION:
!> This submodule tests the case for beta = 0.
!! Note that some predefined equilibria no longer satisfy the equilibrium
!! condition in this limit.
submodule (mod_equilibrium) smod_test_beta0
  implicit none

contains

  module subroutine beta0_test_eq()

    !! Call to the desired equilibrium to test in the beta = 0 limit
    call adiabatic_homo_eq()
    write(*, *) "This test displays the adiabatic homogeneous case in the beta=0 limit."

    T_field % T0      = 0.0d0
    T_field % d_T0_dr = 0.0d0

  end subroutine beta0_test_eq

end submodule smod_test_beta0
