! =============================================================================
!> Module containing the different types used in the code.
!! All types are defined here, this includes all types where the
!! equilibrium arrays are set as attributes, as well as a type
!! to handle eigenfunctions.
module mod_types
  use mod_global_variables, only: dp, str_len, str_len_arr, gauss_gridpts
  implicit none

  private

  !> type containing all density related equilibrium variables
  type density_type
    !> equilibrium density
    real(dp), allocatable   :: rho0(:)
    !> derivative of the equilibrium density
    real(dp), allocatable   :: d_rho0_dr(:)
  end type density_type

  !> type containing all temperature related equilibrium variables
  type temperature_type
    !> equilibrium temperature
    real(dp), allocatable   :: T0(:)
    !> derivative of the equilibrium temperature
    real(dp), allocatable   :: d_T0_dr(:)
    !> second derivative of the equilibrium temperature
    real(dp), allocatable   :: dd_T0_dr(:)
  end type temperature_type

  !> type containing all magnetic field related equilibrium variables
  type bfield_type
    !> constant equilibrium magnetic field in x or r direction
    real(dp)                :: B01
    !> equilibrium magnetic field in y or theta direction
    real(dp), allocatable   :: B02(:)
    !> equilibrium magnetic field in z direction
    real(dp), allocatable   :: B03(:)
    !> equilibrium magnetic field (total)
    real(dp), allocatable   :: B0(:)
    !> derivative of equilibrium B02
    real(dp), allocatable   :: d_B02_dr(:)
    !> derivative of equilibrium B03
    real(dp), allocatable   :: d_B03_dr(:)
  end type bfield_type

  !> type containing all flow related equilibrium variables
  type velocity_type
    !> equilibrium velocity in the x or r direction
    real(dp), allocatable   :: v01(:)
    !> derivative of equilibrium v01
    real(dp), allocatable   :: d_v01_dr(:)
    !> second derivative of equilibrium v01
    real(dp), allocatable   :: dd_v01_dr(:)
    !> equilibrium velocity in the y or theta-direction
    real(dp), allocatable   :: v02(:)
    !> derivative of equilibrium v02
    real(dp), allocatable   :: d_v02_dr(:)
    !> second derivative of equilibrium v02
    real(dp), allocatable   :: dd_v02_dr(:)
    !> equilibrium velocity in the z direction
    real(dp), allocatable   :: v03(:)
    !> derivative of equilibrium v03
    real(dp), allocatable   :: d_v03_dr(:)
    !> second derivative of equilibrium v03
    real(dp), allocatable   :: dd_v03_dr(:)
  end type velocity_type

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
    real(dp), allocatable   :: heat_loss(:)
    !> derivative dL/dT
    real(dp), allocatable   :: d_L_dT(:)
    !> derivative dL/drho
    real(dp), allocatable   :: d_L_drho(:)
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

  !> type containing all eigenfuction related variables
  type ef_type
    !> index of the variable in the state vector
    integer                  :: index
    !> name of the eigenfunction (rho, v1, etc.)
    character(str_len_arr)   :: name
    !> array containing all eigenfunctions for this index
    complex(dp), allocatable :: eigenfunctions(:, :)
  end type ef_type

  !> type containing Hall related variables
  type hall_type
    !> Hall parameter
    real(dp), allocatable    :: hallfactor(:)
    !> electron inertia parameter
    real(dp), allocatable    :: inertiafactor(:)
  end type hall_type


  !> interface to initialise all the different types
  interface initialise_type
    module procedure initialise_density_type
    module procedure initialise_temperature_type
    module procedure initialise_bfield_type
    module procedure initialise_velocity_type
    module procedure initialise_gravity_type
    module procedure initialise_resistivity_type
    module procedure initialise_cooling_type
    module procedure initialise_conduction_type
    module procedure initialise_hall_type
  end interface initialise_type

  !> interface to deallocate all the different types
  interface deallocate_type
    module procedure deallocate_density_type
    module procedure deallocate_temperature_type
    module procedure deallocate_bfield_type
    module procedure deallocate_velocity_type
    module procedure deallocate_gravity_type
    module procedure deallocate_resistivity_type
    module procedure deallocate_cooling_type
    module procedure deallocate_conduction_type
    module procedure deallocate_hall_type
  end interface deallocate_type

  public :: density_type
  public :: temperature_type
  public :: bfield_type
  public :: velocity_type
  public :: gravity_type
  public :: resistivity_type
  public :: cooling_type
  public :: conduction_type
  public :: ef_type
  public :: hall_type

  public :: initialise_type
  public :: deallocate_type

contains


  !> Allocates the density type and initialises all values to zero.
  subroutine initialise_density_type(type_rho)
    !> the type containing the density attributes
    type (density_type), intent(inout)  :: type_rho

    allocate(type_rho % rho0(gauss_gridpts))
    allocate(type_rho % d_rho0_dr(gauss_gridpts))

    type_rho % rho0 = 0.0d0
    type_rho % d_rho0_dr = 0.0d0
  end subroutine initialise_density_type


  !> Allocates the temperature type and initialises all values to zero.
  subroutine initialise_temperature_type(type_T)
    !> the type containing the temperature attributes
    type (temperature_type), intent(inout) :: type_T

    allocate(type_T % T0(gauss_gridpts))
    allocate(type_T % d_T0_dr(gauss_gridpts))
    allocate(type_T % dd_T0_dr(gauss_gridpts))

    type_T % T0 = 0.0d0
    type_T % d_T0_dr = 0.0d0
    type_T % dd_T0_dr = 0.0d0
  end subroutine initialise_temperature_type


  !> Allocates the magnetic field type and initialises all values to zero.
  subroutine initialise_bfield_type(type_B)
    !> the type containing the magnetic field attributes
    type (bfield_type), intent(inout) :: type_B

    allocate(type_B % B02(gauss_gridpts))
    allocate(type_B % B03(gauss_gridpts))
    allocate(type_B % B0(gauss_gridpts))
    allocate(type_B % d_B02_dr(gauss_gridpts))
    allocate(type_B % d_B03_dr(gauss_gridpts))

    type_B % B01 = 0.0d0
    type_B % B02 = 0.0d0
    type_B % B03 = 0.0d0
    type_B % B0 = 0.0d0
    type_B % d_B02_dr = 0.0d0
    type_B % d_B03_dr = 0.0d0
  end subroutine initialise_bfield_type


  !> Allocates the velocity type and initialises all values to zero.
  subroutine initialise_velocity_type(type_v)
    !> the type containing the velocity attributes
    type (velocity_type), intent(inout) :: type_v

    allocate(type_v % v01(gauss_gridpts))
    allocate(type_v % d_v01_dr(gauss_gridpts))
    allocate(type_v % dd_v01_dr(gauss_gridpts))
    allocate(type_v % v02(gauss_gridpts))
    allocate(type_v % d_v02_dr(gauss_gridpts))
    allocate(type_v % dd_v02_dr(gauss_gridpts))
    allocate(type_v % v03(gauss_gridpts))
    allocate(type_v % d_v03_dr(gauss_gridpts))
    allocate(type_v % dd_v03_dr(gauss_gridpts))

    type_v % v01 = 0.0d0
    type_v % d_v01_dr = 0.0d0
    type_v % dd_v01_dr = 0.0d0
    type_v % v02 = 0.0d0
    type_v % d_v02_dr = 0.0d0
    type_v % dd_v02_dr = 0.0d0
    type_v % v03 = 0.0d0
    type_v % d_v03_dr = 0.0d0
    type_v % dd_v03_dr = 0.0d0
  end subroutine initialise_velocity_type


  !> Allocates the gravity type and initialises all values to zero.
  subroutine initialise_gravity_type(type_grav)
    !> the type containing the gravity attributes
    type (gravity_type), intent(inout)  :: type_grav

    allocate(type_grav % grav(gauss_gridpts))

    type_grav % grav = 0.0d0
  end subroutine initialise_gravity_type


  !> Allocates the resistivity type and initialises all values to zero.
  !! @note The second derivatives of the magnetic field components
  !!       are also included in this type, since they are only used
  !!       when resistivity is included.
  subroutine initialise_resistivity_type(type_eta)
    !> the type containing the resistivity attributes
    type (resistivity_type), intent(inout)  :: type_eta

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
  subroutine initialise_cooling_type(type_rc)
    !> the type containing the radiative cooling attributes
    type (cooling_type), intent(inout)  :: type_rc

    allocate(type_rc % heat_loss(gauss_gridpts))
    allocate(type_rc % d_L_dT(gauss_gridpts))
    allocate(type_rc % d_L_drho(gauss_gridpts))

    type_rc % heat_loss = 0.0d0
    type_rc % d_L_dT = 0.0d0
    type_rc % d_L_drho = 0.0d0
  end subroutine initialise_cooling_type


  !> Allocates the Hall type and initialises all values to zero.
  subroutine initialise_hall_type(type_hall)
    !> the type containing the density attributes
    type (hall_type), intent(inout)  :: type_hall

    allocate(type_hall % hallfactor(gauss_gridpts))
    allocate(type_hall % inertiafactor(gauss_gridpts))

    type_hall % hallfactor = 0.0d0
    type_hall % inertiafactor = 0.0d0
  end subroutine initialise_hall_type


  !> Allocates the thermal conduction type and initialises all values to zero.
  subroutine initialise_conduction_type(type_kappa)
    !> the type containing the thermal conduction attributes
    type (conduction_type), intent(inout) :: type_kappa

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


  !> Deallocates all attributes contained in the  density type.
  subroutine deallocate_density_type(type_rho)
    !> the type containing the density attributes
    type (density_type), intent(inout)  :: type_rho

    deallocate(type_rho % rho0)
    deallocate(type_rho % d_rho0_dr)
  end subroutine deallocate_density_type


  !> Deallocates all attributes contained in the temperature type.
  subroutine deallocate_temperature_type(type_T)
    !> the type containing the temperature attributes
    type (temperature_type), intent(inout)  :: type_T

    deallocate(type_T % T0)
    deallocate(type_T % d_T0_dr)
    deallocate(type_T % dd_T0_dr)
  end subroutine deallocate_temperature_type


  !> Deallocates all attributes contained in the magnetic field type.
  subroutine deallocate_bfield_type(type_B)
    !> the type containing the magnetic field attributes
    type (bfield_type), intent(inout) :: type_B

    deallocate(type_B % B02)
    deallocate(type_B % B03)
    deallocate(type_B % B0)
    deallocate(type_B % d_B02_dr)
    deallocate(type_B % d_B03_dr)
  end subroutine deallocate_bfield_type


  !> Deallocates all attributes contained in the velocity type.
  subroutine deallocate_velocity_type(type_v)
    !> the type containing the velocity attributes
    type (velocity_type), intent(inout) :: type_v

    deallocate(type_v % v01)
    deallocate(type_v % d_v01_dr)
    deallocate(type_v % dd_v01_dr)
    deallocate(type_v % v02)
    deallocate(type_v % d_v02_dr)
    deallocate(type_v % dd_v02_dr)
    deallocate(type_v % v03)
    deallocate(type_v % d_v03_dr)
    deallocate(type_v % dd_v03_dr)
  end subroutine deallocate_velocity_type


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

    deallocate(type_rc % heat_loss)
    deallocate(type_rc % d_L_dT)
    deallocate(type_rc % d_L_drho)
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
