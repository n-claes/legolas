! =============================================================================
!> @brief   Module containing all physical constants
!! @details All physical constants used in the code are defined in this module.
!!          We include values both in SI units and in cgs units for convenience.
module mod_physical_constants
  use mod_global_variables, only: dp
  implicit none

  public

  !> value for pi
  real(dp), parameter   :: dpi = 3.141592653589793238462643383279d0
  !> coulomb logarithm
  real(dp), parameter   :: coulomb_log = 22.0d0
  !> proton mass in g (cgs)
  real(dp), parameter   :: mp_cgs = 1.672621777d-24
  !> proton mass in kg (SI)
  real(dp), parameter   :: mp_si  = 1.672621777d-27
  !> hydrogen mass in g (cgs)
  real(dp), parameter   :: mH_cgs = 1.6733d-24
  !> hydrogen mass in kg (SI)
  real(dp), parameter   :: mH_si  = 1.6733d-27
  !> electron mass in g (cgs)
  real(dp), parameter   :: me_cgs = 9.1094d-28
  !> electron mass in kg (SI)
  real(dp), parameter   :: me_si  = 9.1094d-31
  !> elementary charge in statcoul (cgs)
  real(dp), parameter   :: ec_cgs = 4.8032d-10
  !> elementary charge in C (SI)
  real(dp), parameter   :: ec_si  = 1.6022d-19
  !> Boltzmann constant in erg/K (cgs)
  real(dp), parameter   :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K (SI)
  real(dp), parameter   :: kB_si  = 1.3806488d-23
  !> magnetic constant in H/m (SI)
  real(dp), parameter   :: mu0_si = 1.2566370614d-6
  !> magnetic constant (cgs)
  real(dp), parameter   :: mu0_cgs = 4.0d0*dpi
  !> permittivity of free space in F/m (SI)
  real(dp), parameter   :: e0_si  = 8.8542d-12
  !> degree of ionization (assumed fully ionized)
  real(dp), parameter   :: Z_ion = 1.0d0
  !> gas constant in J/K (SI)
  real(dp), parameter   :: R_si = 8.3145d0
  !> gas constant in erg/deg
  real(dp), parameter   :: R_cgs = 8.3145d7

end module mod_physical_constants
