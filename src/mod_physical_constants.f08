!
! MODULE: mod_physical_constants
!
! DESCRIPTION:
!> Module containing all relevant physical constants.
!
module mod_physical_constants
  use mod_global_variables, only: dp
  implicit none

  public

  !> Value for pi
  real(dp), parameter   :: dpi = 3.141592653589793238462643383279d0
  !> Coulomb logarithm
  real(dp), parameter   :: coulomb_log = 22.0d0
  !> Proton mass in g (cgs)
  real(dp), parameter   :: mp_cgs = 1.672621777d-24
  !> Proton mass in kg (SI)
  real(dp), parameter   :: mp_si  = 1.672621777d-27
  !> Hydrogen mass in g (cgs)
  real(dp), parameter   :: mH_cgs = 1.6733d-24
  !> Hydrogen mass in kg (SI)
  real(dp), parameter   :: mH_si  = 1.6733d-27
  !> Electron mass in g (cgs)
  real(dp), parameter   :: me_cgs = 9.1094d-28
  !> Electron mass in kg (SI)
  real(dp), parameter   :: me_si  = 9.1094d-31
  !> Elementary charge in statcoul (cgs)
  real(dp), parameter   :: ec_cgs = 4.8032d-10
  !> Elementary charge in C (SI)
  real(dp), parameter   :: ec_si  = 1.6022d-19
  !> Boltzmann constant in erg/K (cgs)
  real(dp), parameter   :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K (SI)
  real(dp), parameter   :: kB_si  = 1.3806488d-23
  !> Magnetic constant in H/m (SI)
  real(dp), parameter   :: mu0_si = 1.2566370614d-6
  !> Magnetic constant (cgs)
  real(dp), parameter   :: mu0_cgs = 4.0d0*dpi
  !> Permittivity of free space in F/m (SI)
  real(dp), parameter   :: e0_si  = 8.8542d-12
  !> Degree of ionization (assumed fully ionized)
  real(dp), parameter   :: Z_ion = 1.0d0
  !> Gas constant in J/K (SI)
  real(dp), parameter   :: R_si = 8.3145d0
  !> Gas constant in erg/deg
  real(dp), parameter   :: R_cgs = 8.3145d7

end module mod_physical_constants
