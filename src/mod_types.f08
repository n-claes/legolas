module mod_types
  use mod_global_variables, only: dp, str_len, str_len_arr, gauss_gridpts
  implicit none

  private

  !> Type containing all density related equilibrium variables
  type density_type
    !> Equilibrium density
    real(dp), allocatable   :: rho0(:)
    !> Derivative rho0
    real(dp), allocatable   :: d_rho0_dr(:)
  end type density_type

  !> Type containing all temperature related equilibrium variables
  type temperature_type
    !> Equilibrium temperature
    real(dp), allocatable   :: T0(:)
    !> Derivative T0
    real(dp), allocatable   :: d_T0_dr(:)
  end type temperature_type

  !> Type containing all magnetic field related equilibrium
  !! variables
  type bfield_type
    !> Equilibrium magnetic field in y or theta-direction
    real(dp), allocatable   :: B02(:)
    !> Equilibrium magnetic field in z direction
    real(dp), allocatable   :: B03(:)
    !> Equilibrium magnetic field (total)
    real(dp), allocatable   :: B0(:)
    !> Derivative B02
    real(dp), allocatable   :: d_B02_dr(:)
    !> Derivative B03
    real(dp), allocatable   :: d_B03_dr(:)
  end type bfield_type

  !> Type containing all flow related equilibrium variables
  type velocity_type
    !> Equilibrium velocity in the y or theta-direction
    real(dp), allocatable   :: v02(:)
    !> Equilibrium velocity in the z direction
    real(dp), allocatable   :: v03(:)
    !> Derivative v02
    real(dp), allocatable   :: d_v02_dr(:)
    !> Derivative v03
    real(dp), allocatable   :: d_v03_dr(:)
  end type velocity_type

  !> Type containing all gravity related equilibrium
  !! variables
  type gravity_type
    !> Equilibrium gravity
    real(dp), allocatable   :: grav(:)
  end type gravity_type

  !> Type containing all resistivity related equilibrium
  !! variables
  type resistivity_type
    !> Equilibrium resistivity
    real(dp), allocatable   :: eta(:)
    !> Derivative of eta with respect to temperature
    real(dp), allocatable   :: d_eta_dT(:)
    !> Second derivative B02
    real(dp), allocatable   :: dd_B02_dr(:)
    !> Second derivative B03
    real(dp), allocatable   :: dd_B03_dr(:)
  end type resistivity_type

  !> Type containing all radiative cooling related
  !! equilibrium variables
  type cooling_type
    !> Equilibrium heat-loss function (losses - gains)
    real(dp), allocatable   :: heat_loss(:)
    !> Derivative dL/dT
    real(dp), allocatable   :: d_L_dT(:)
    !> Derivative dL/drho
    real(dp), allocatable   :: d_L_drho(:)
  end type cooling_type

  !> Type containing all conduction related equilibrium
  !! variables
  type conduction_type
    !> Equilibrium parallel thermal conduction
    real(dp), allocatable   :: kappa_para(:)
    !> Equilibrium perpendicular thermal conduction
    real(dp), allocatable   :: kappa_perp(:)
    !> Derivative d(kappa_perp)/drho
    real(dp), allocatable   :: d_kappa_perp_drho(:)
    !> Derivative d(kappa_perp)/dT
    real(dp), allocatable   :: d_kappa_perp_dT(:)
    !> Derivative d(kappa_perp)/d(B**2)
    real(dp), allocatable   :: d_kappa_perp_dB2(:)
  end type conduction_type

  !> Type containing all eigenfuction related variables
  type ef_type
    !> Index of the eigenfunction (1 -> 8)
    integer                  :: index
    !> Name of the eigenfunction
    character(str_len_arr)   :: name
    !> Array containing all eigenfunctions for this index
    complex(dp), allocatable :: eigenfunctions(:, :)
  end type ef_type


  interface initialise_type
    module procedure initialise_density_type
    module procedure initialise_temperature_type
    module procedure initialise_bfield_type
    module procedure initialise_velocity_type
    module procedure initialise_gravity_type
    module procedure initialise_resistivity_type
    module procedure initialise_cooling_type
    module procedure initialise_conduction_type
  end interface initialise_type

  interface deallocate_type
    module procedure deallocate_density_type
    module procedure deallocate_temperature_type
    module procedure deallocate_bfield_type
    module procedure deallocate_velocity_type
    module procedure deallocate_gravity_type
    module procedure deallocate_resistivity_type
    module procedure deallocate_cooling_type
    module procedure deallocate_conduction_type
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

  public :: initialise_type
  public :: deallocate_type

contains

  subroutine initialise_density_type(type_rho)
    type (density_type), intent(inout)  :: type_rho

    allocate(type_rho % rho0(gauss_gridpts))
    allocate(type_rho % d_rho0_dr(gauss_gridpts))

    type_rho % rho0 = 0.0d0
    type_rho % d_rho0_dr = 0.0d0
  end subroutine initialise_density_type


  subroutine initialise_temperature_type(type_T)
    type (temperature_type), intent(inout) :: type_T

    allocate(type_T % T0(gauss_gridpts))
    allocate(type_T % d_T0_dr(gauss_gridpts))

    type_T % T0 = 0.0d0
    type_T % d_T0_dr = 0.0d0
  end subroutine initialise_temperature_type


  subroutine initialise_bfield_type(type_B)
    type (bfield_type), intent(inout) :: type_B

    allocate(type_B % B02(gauss_gridpts))
    allocate(type_B % B03(gauss_gridpts))
    allocate(type_B % B0(gauss_gridpts))
    allocate(type_B % d_B02_dr(gauss_gridpts))
    allocate(type_B % d_B03_dr(gauss_gridpts))

    type_B % B02 = 0.0d0
    type_B % B03 = 0.0d0
    type_B % B0 = 0.0d0
    type_B % d_B02_dr = 0.0d0
    type_B % d_B03_dr = 0.0d0
  end subroutine initialise_bfield_type


  subroutine initialise_velocity_type(type_v)
    type (velocity_type), intent(inout) :: type_v

    allocate(type_v % v02(gauss_gridpts))
    allocate(type_v % v03(gauss_gridpts))
    allocate(type_v % d_v02_dr(gauss_gridpts))
    allocate(type_v % d_v03_dr(gauss_gridpts))

    type_v % v02 = 0.0d0
    type_v % v03 = 0.0d0
    type_v % d_v02_dr = 0.0d0
    type_v % d_v03_dr = 0.0d0
  end subroutine initialise_velocity_type


  subroutine initialise_gravity_type(type_grav)
    type (gravity_type), intent(inout)  :: type_grav

    allocate(type_grav % grav(gauss_gridpts))

    type_grav % grav = 0.0d0
  end subroutine initialise_gravity_type


  subroutine initialise_resistivity_type(type_eta)
    type (resistivity_type), intent(inout)  :: type_eta

    allocate(type_eta % eta(gauss_gridpts))
    allocate(type_eta % d_eta_dT(gauss_gridpts))
    allocate(type_eta % dd_B02_dr(gauss_gridpts))
    allocate(type_eta % dd_B03_dr(gauss_gridpts))

    type_eta % eta = 0.0d0
    type_eta % d_eta_dT = 0.0d0
    type_eta % dd_B02_dr = 0.0d0
    type_eta % dd_B03_dr = 0.0d0
  end subroutine initialise_resistivity_type


  subroutine initialise_cooling_type(type_rc)
    type (cooling_type), intent(inout)  :: type_rc

    allocate(type_rc % heat_loss(gauss_gridpts))
    allocate(type_rc % d_L_dT(gauss_gridpts))
    allocate(type_rc % d_L_drho(gauss_gridpts))

    type_rc % heat_loss = 0.0d0
    type_rc % d_L_dT = 0.0d0
    type_rc % d_L_drho = 0.0d0
  end subroutine initialise_cooling_type


  subroutine initialise_conduction_type(type_kappa)
    type (conduction_type), intent(inout) :: type_kappa

    allocate(type_kappa % kappa_para(gauss_gridpts))
    allocate(type_kappa % kappa_perp(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_drho(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_dT(gauss_gridpts))
    allocate(type_kappa % d_kappa_perp_dB2(gauss_gridpts))

    type_kappa % kappa_para = 0.0d0
    type_kappa % kappa_perp = 0.0d0
    type_kappa % d_kappa_perp_drho = 0.0d0
    type_kappa % d_kappa_perp_dT = 0.0d0
    type_kappa % d_kappa_perp_dB2 = 0.0d0
  end subroutine initialise_conduction_type


  subroutine deallocate_density_type(type_rho)
    type (density_type), intent(inout)  :: type_rho

    deallocate(type_rho % rho0)
    deallocate(type_rho % d_rho0_dr)
  end subroutine deallocate_density_type


  subroutine deallocate_temperature_type(type_T)
    type (temperature_type), intent(inout)  :: type_T

    deallocate(type_T % T0)
    deallocate(type_T % d_T0_dr)
  end subroutine deallocate_temperature_type


  subroutine deallocate_bfield_type(type_B)
    type (bfield_type), intent(inout) :: type_B

    deallocate(type_B % B02)
    deallocate(type_B % B03)
    deallocate(type_B % B0)
    deallocate(type_B % d_B02_dr)
    deallocate(type_B % d_B03_dr)
  end subroutine deallocate_bfield_type


  subroutine deallocate_velocity_type(type_v)
    type (velocity_type), intent(inout) :: type_v

    deallocate(type_v % v02)
    deallocate(type_v % v03)
    deallocate(type_v % d_v02_dr)
    deallocate(type_v % d_v03_dr)
  end subroutine deallocate_velocity_type


  subroutine deallocate_gravity_type(type_grav)
    type (gravity_type), intent(inout) :: type_grav

    deallocate(type_grav % grav)
  end subroutine deallocate_gravity_type


  subroutine deallocate_resistivity_type(type_eta)
    type (resistivity_type), intent(inout) :: type_eta

    deallocate(type_eta % eta)
    deallocate(type_eta % d_eta_dT)
    deallocate(type_eta % dd_B02_dr)
    deallocate(type_eta % dd_B03_dr)
  end subroutine deallocate_resistivity_type


  subroutine deallocate_cooling_type(type_rc)
    type (cooling_type), intent(inout) :: type_rc

    deallocate(type_rc % heat_loss)
    deallocate(type_rc % d_L_dT)
    deallocate(type_rc % d_L_drho)
  end subroutine deallocate_cooling_type


  subroutine deallocate_conduction_type(type_kappa)
    type (conduction_type), intent(inout)  :: type_kappa

    deallocate(type_kappa % kappa_para)
    deallocate(type_kappa % kappa_perp)
    deallocate(type_kappa % d_kappa_perp_drho)
    deallocate(type_kappa % d_kappa_perp_dT)
    deallocate(type_kappa % d_kappa_perp_dB2)
  end subroutine deallocate_conduction_type

end module mod_types
