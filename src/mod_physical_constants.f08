! =============================================================================
!> All physical constants used in the code are defined in this module.
!! We include values both in SI units and in cgs units for convenience.
!! All values are taken from the
!! [NRL Plasma Formulary](https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary).
module mod_physical_constants
  use mod_global_variables, only: dp
  implicit none

  public

  !> value for pi
  real(dp), parameter   :: dpi = 3.141592653589793238462643383279d0
  !> coulomb logarithm
  real(dp), parameter   :: coulomb_log = 22.0d0
  !> proton mass in g
  real(dp), parameter   :: mp_cgs = 1.672621777d-24
  !> proton mass in kg
  real(dp), parameter   :: mp_si  = 1.672621777d-27
  !> hydrogen mass in g
  real(dp), parameter   :: mH_cgs = 1.6733d-24
  !> hydrogen mass in kg
  real(dp), parameter   :: mH_si  = 1.6733d-27
  !> electron mass in g
  real(dp), parameter   :: me_cgs = 9.1094d-28
  !> electron mass in kg
  real(dp), parameter   :: me_si  = 9.1094d-31
  !> elementary charge in statcoul
  real(dp), parameter   :: ec_cgs = 4.8032d-10
  !> elementary charge in C
  real(dp), parameter   :: ec_si  = 1.6022d-19
  !> Boltzmann constant in erg/K
  real(dp), parameter   :: kB_cgs = 1.3806488d-16
  !> Boltzmann constant in J/K
  real(dp), parameter   :: kB_si  = 1.3806488d-23
  !> magnetic constant in H/m
  real(dp), parameter   :: mu0_si = 1.2566370614d-6
  !> magnetic constant
  real(dp), parameter   :: mu0_cgs = 4.0d0*dpi
  !> permittivity of free space in F/m
  real(dp), parameter   :: e0_si  = 8.8542d-12
  !> degree of ionization (assumed fully ionized)
  real(dp), parameter   :: Z_ion = 1.0d0
  !> gas constant in J/K
  real(dp), parameter   :: R_si = 8.3145d0
  !> gas constant in erg/deg
  real(dp), parameter   :: R_cgs = 8.3145d7

  !! Solar physics parameters
  !> solar mass in g
  real(dp), parameter   :: msun_cgs = 1.99d33
  !> solar radius in cm
  real(dp), parameter   :: Rsun_cgs = 6.96d10
  !> solar gravity in cm/s^2
  real(dp), parameter   :: gsun_cgs = 2.74d4
  !> solar escape velocity in cm/s
  real(dp), parameter   :: Vesc_sun = 6.18d7
  !> solar effective temperature in K
  real(dp), parameter   :: Teff_sun = 5770.0d0
  !> solar luminosity in erg/s
  real(dp), parameter   :: Lsun_cgs = 3.83d33
  !> astronomical unit in cm
  real(dp), parameter   :: AU_cgs = 1.50d13
  !> solar constant (intensity at 1 AU) in erg/cm^2/s
  real(dp), parameter   :: fsun_cgs = 1.36d6

  !> coefficient of parallel thermal conduction in erg*s/cm/K
  !! the coulomb logarithm ln(Lambda) has not yet been included
  real(dp), parameter :: tc_pf_kappa_para = 1.8d-5
  !> coefficient of perpendicular thermal conduction in erg*s/cm/K
  !! the coulomb logarithm ln(Lambda) has not yet been included
  real(dp), parameter :: tc_pf_kappa_perp = 8.2d-13

end module mod_physical_constants
