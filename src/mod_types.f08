! =============================================================================
!> Module containing the different types used in the code.
!! All types are defined here, this includes all types where the
!! equilibrium arrays are set as attributes, as well as a type
!! to handle eigenfunctions.
module mod_types
  use mod_global_variables, only: dp, str_len, str_len_arr
  use mod_settings, only: settings_t
  implicit none

  private

  !> type containing all gravity related equilibrium variables
  type gravity_type
    !> equilibrium gravity
    real(dp), allocatable   :: grav(:)
  end type gravity_type

  !> type containing all resistivity related equilibrium variables
  type resistivity_type
    !> equilibrium resistivity
    real(dp), allocatable   :: eta(:)
    !> derivative of eta with respect to temperature
    real(dp), allocatable   :: d_eta_dT(:)
    !> derivative of eta with respect to radius
    real(dp), allocatable   :: d_eta_dr(:)
    !> second derivative equilibrium B02
    real(dp), allocatable   :: dd_B02_dr(:)
    !> second derivative equilibrium B03
    real(dp), allocatable   :: dd_B03_dr(:)
  end type resistivity_type

  !> type containing all radiative cooling related equilibrium variables
  type cooling_type
    !> equilibrium heat-loss function L0 (losses - gains)
    real(dp), allocatable :: L0(:)
    !> derivative dL/dT
    real(dp), allocatable :: dL_dT(:)
    !> derivative dL/drho
    real(dp), allocatable :: dL_drho(:)
    real(dp), allocatable :: H0(:)
    real(dp), allocatable :: dH_drho(:)
    real(dp), allocatable :: dH_dT(:)
    real(dp), allocatable :: lambdaT(:)
  end type cooling_type

  !> Type containing all conduction related equilibrium variables
  type conduction_type
    !> equilibrium parallel thermal conduction
    real(dp), allocatable   :: kappa_para(:)
    !> derivative d(kappa_para)/dT
    real(dp), allocatable   :: d_kappa_para_dT(:)
    !> equilibrium perpendicular thermal conduction
    real(dp), allocatable   :: kappa_perp(:)
    !> derivative d(kappa_perp)/drho
    real(dp), allocatable   :: d_kappa_perp_drho(:)
    !> derivative d(kappa_perp)/dT
    real(dp), allocatable   :: d_kappa_perp_dT(:)
    !> derivative d(kappa_perp)/d(B**2)
    real(dp), allocatable   :: d_kappa_perp_dB2(:)
    !> derivative d(kappa_perp)/dr
    real(dp), allocatable   :: d_kappa_perp_dr(:)
    !> conduction prefactor (kappa_para - kappa_perp) / B**2
    real(dp), allocatable   :: prefactor(:)
    !> derivative conduction prefactor with respect to r
    real(dp), allocatable   :: d_prefactor_dr(:)
  end type conduction_type

  !> type containing Hall related variables
  type hall_type
    !> Hall parameter
    real(dp), allocatable    :: hallfactor(:)
    !> electron inertia parameter
    real(dp), allocatable    :: inertiafactor(:)
  end type hall_type


  !> interface to initialise all the different types
  interface initialise_type
    module procedure initialise_gravity_type
    module procedure initialise_resistivity_type
    module procedure initialise_cooling_type
    module procedure initialise_conduction_type
    module procedure initialise_hall_type
  end interface initialise_type

  !> interface to deallocate all the different types
  interface deallocate_type
    module procedure deallocate_gravity_type
    module procedure deallocate_resistivity_type
    module procedure deallocate_cooling_type
    module procedure deallocate_conduction_type
    module procedure deallocate_hall_type
  end interface deallocate_type

  public :: gravity_type
  public :: resistivity_type
  public :: cooling_type
  public :: conduction_type
  public :: hall_type

  public :: initialise_type
  public :: deallocate_type

contains


  !> Allocates the gravity type and initialises all values to zero.
  subroutine initialise_gravity_type(settings, type_grav)
    type(settings_t), intent(in) :: settings
    !> the type containing the gravity attributes
    type (gravity_type), intent(inout)  :: type_grav
    integer :: gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()

    allocate(type_grav % grav(gauss_gridpts))

    type_grav % grav = 0.0d0
  end subroutine initialise_gravity_type


  !> Allocates the resistivity type and initialises all values to zero.
  !! @note The second derivatives of the magnetic field components
  !!       are also included in this type, since they are only used
  !!       when resistivity is included.
  subroutine initialise_resistivity_type(settings, type_eta)
    type(settings_t), intent(in) :: settings
    !> the type containing the resistivity attributes
    type (resistivity_type), intent(inout)  :: type_eta
    integer :: gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()

    allocate(type_eta % eta(gauss_gridpts))
    allocate(type_eta % d_eta_dT(gauss_gridpts))
    allocate(type_eta % d_eta_dr(gauss_gridpts))
    allocate(type_eta % dd_B02_dr(gauss_gridpts))
    allocate(type_eta % dd_B03_dr(gauss_gridpts))

    type_eta % eta = 0.0d0
    type_eta % d_eta_dT = 0.0d0
    type_eta % d_eta_dr = 0.0d0
    type_eta % dd_B02_dr = 0.0d0
    type_eta % dd_B03_dr = 0.0d0
  end subroutine initialise_resistivity_type


  !> Allocates the radiative cooling type and initialises all values to zero.
  subroutine initialise_cooling_type(settings, type_rc)
    type(settings_t), intent(in) :: settings
    !> the type containing the radiative cooling attributes
    type (cooling_type), intent(inout)  :: type_rc
    integer :: gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()

    allocate(type_rc % L0(gauss_gridpts))
    allocate(type_rc % dL_dT(gauss_gridpts))
    allocate(type_rc % dL_drho(gauss_gridpts))
    allocate(type_rc % H0(gauss_gridpts))
    allocate(type_rc % dH_drho(gauss_gridpts))
    allocate(type_rc % dH_dT(gauss_gridpts))
    allocate(type_rc % lambdaT(gauss_gridpts))

    type_rc % L0 = 0.0d0
    type_rc % dL_dT = 0.0d0
    type_rc % dL_drho = 0.0d0
    type_rc % H0 = 0.0d0
    type_rc % dH_drho = 0.0d0
    type_rc % dH_dT = 0.0d0
    type_rc % lambdaT = 0.0d0
  end subroutine initialise_cooling_type


  !> Allocates the Hall type and initialises all values to zero.
  subroutine initialise_hall_type(settings, type_hall)
    type(settings_t), intent(in) :: settings
    !> the type containing the density attributes
    type (hall_type), intent(inout)  :: type_hall
    integer :: gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()

    allocate(type_hall % hallfactor(gauss_gridpts))
    allocate(type_hall % inertiafactor(gauss_gridpts))

    type_hall % hallfactor = 0.0d0
    type_hall % inertiafactor = 0.0d0
  end subroutine initialise_hall_type


  !> Allocates the thermal conduction type and initialises all values to zero.
  subroutine initialise_conduction_type(settings, type_kappa)
    type(settings_t), intent(in) :: settings
    !> the type containing the thermal conduction attributes
    type (conduction_type), intent(inout) :: type_kappa
    integer :: gauss_gridpts

    gauss_gridpts = settings%grid%get_gauss_gridpts()

    allocate(type_kappa % kappa_para(gauss_gridpts))
    allocate(type_kappa % d_kappa_para_dT(gauss_gridpts))
    allocate(type_kappa % kappa_perp(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_drho(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_dT(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_dB2(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_dr(gauss_gridpts))
    allocate(type_kappa % prefactor(gauss_gridpts))
    allocate(type_kappa % d_prefactor_dr(gauss_gridpts))

    type_kappa % kappa_para = 0.0d0
    type_kappa % d_kappa_para_dT = 0.0d0
    type_kappa % kappa_perp = 0.0d0
    type_kappa % d_kappa_perp_drho = 0.0d0
    type_kappa % d_kappa_perp_dT = 0.0d0
    type_kappa % d_kappa_perp_dB2 = 0.0d0
    type_kappa % d_kappa_perp_dr = 0.0d0
    type_kappa % prefactor = 0.0d0
    type_kappa % d_prefactor_dr = 0.0d0
  end subroutine initialise_conduction_type


  !> Deallocates all attributes contained in the gravity type.
  subroutine deallocate_gravity_type(type_grav)
    !> the type containing the gravity attributes
    type (gravity_type), intent(inout) :: type_grav

    deallocate(type_grav % grav)
  end subroutine deallocate_gravity_type


  !> Deallocates all attributes contained in the resistivity type.
  subroutine deallocate_resistivity_type(type_eta)
    !> the type containing the resistivity attributes
    type (resistivity_type), intent(inout) :: type_eta

    deallocate(type_eta % eta)
    deallocate(type_eta % d_eta_dT)
    deallocate(type_eta % d_eta_dr)
    deallocate(type_eta % dd_B02_dr)
    deallocate(type_eta % dd_B03_dr)
  end subroutine deallocate_resistivity_type


  !> Deallocates all attributes contained in the radiative cooling type.
  subroutine deallocate_cooling_type(type_rc)
    !> the type containing the radiative cooling attributes
    type (cooling_type), intent(inout) :: type_rc

    deallocate(type_rc % L0)
    deallocate(type_rc % dL_dT)
    deallocate(type_rc % dL_drho)
    deallocate(type_rc % H0)
    deallocate(type_rc % dH_drho)
    deallocate(type_rc % dH_dT)
    deallocate(type_rc % lambdaT)
  end subroutine deallocate_cooling_type


  !> Deallocates all attributes contained in the thermal conduction type.
  subroutine deallocate_conduction_type(type_kappa)
    !> the type containing the thermal conduction attributes
    type (conduction_type), intent(inout)  :: type_kappa

    deallocate(type_kappa % kappa_para)
    deallocate(type_kappa % d_kappa_para_dT)
    deallocate(type_kappa % kappa_perp)
    deallocate(type_kappa % d_kappa_perp_drho)
    deallocate(type_kappa % d_kappa_perp_dT)
    deallocate(type_kappa % d_kappa_perp_dB2)
    deallocate(type_kappa % d_kappa_perp_dr)
    deallocate(type_kappa % prefactor)
    deallocate(type_kappa % d_prefactor_dr)
  end subroutine deallocate_conduction_type


  !> Deallocates all attributes contained in the Hall type.
  subroutine deallocate_hall_type(type_hall)
    !> the type containing the density attributes
    type (hall_type), intent(inout)  :: type_hall

    deallocate(type_hall % hallfactor)
    deallocate(type_hall % inertiafactor)
  end subroutine deallocate_hall_type

end module mod_types
