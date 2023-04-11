module mod_heatloss
  use mod_global_variables, only: dp
  use mod_logging, only: logger
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_radiative_cooling, only: cooling_t
  use mod_heating, only: heating_t
  use mod_thermal_conduction, only: conduction_t
  use mod_grid, only: grid_t
  implicit none
  private

  type(settings_t), pointer :: settings => null()
  type(background_t), pointer :: background => null()
  type(cooling_t), pointer :: cooling =>  null()
  type(heating_t), pointer :: heating => null()
  type(conduction_t), pointer :: conduction => null()
  type(grid_t), pointer :: grid => null()

  type, public :: heatloss_t
    procedure(real(dp)), pointer, nopass :: L0
    procedure(real(dp)), pointer, nopass :: dLdT
    procedure(real(dp)), pointer, nopass :: dLdrho

  contains
    procedure, public :: check_if_thermal_balance_needs_enforcing
    procedure, public :: delete
  end type heatloss_t

  public :: new_heatloss

contains

  function new_heatloss( &
    settings_tgt, background_tgt, cooling_tgt, heating_tgt &
  ) result(heatloss)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(cooling_t), target, intent(in) :: cooling_tgt
    type(heating_t), target, intent(in) :: heating_tgt
    type(heatloss_t) :: heatloss

    settings => settings_tgt
    background => background_tgt
    cooling => cooling_tgt
    heating => heating_tgt

    heatloss%L0 => get_L0
    heatloss%dLdT => get_dLdT
    heatloss%dLdrho => get_dLdrho
  end function new_heatloss


  real(dp) function get_L0(x)
    real(dp), intent(in) :: x
    get_L0 = background%density%rho0(x) * cooling%lambdaT(x) - heating%H(x)
  end function get_L0


  real(dp) function get_dLdT(x)
    real(dp), intent(in) :: x
    get_dLdT = background%density%rho0(x) * cooling%dlambdadT(x) - heating%dHdT(x)
  end function get_dLdT


  real(dp) function get_dLdrho(x)
    real(dp), intent(in) :: x
    get_dLdrho = cooling%lambdaT(x) - heating%dHdrho(x)
  end function get_dLdrho


  subroutine check_if_thermal_balance_needs_enforcing(this, conduction_tgt, grid_tgt)
    class(heatloss_t), intent(inout) :: this
    type(conduction_t), target, intent(in) :: conduction_tgt
    type(grid_t), target, intent(in) :: grid_tgt

    if (.not. settings%physics%heating%is_enabled()) return
    if (.not. settings%physics%heating%force_thermal_balance) return

    ! needed for eps and deps + physics in function below
    grid => grid_tgt
    conduction => conduction_tgt

    call logger%info("enforcing thermal balance by setting a constant heating term")
    heating%H => H_for_thermal_balance
  end subroutine check_if_thermal_balance_needs_enforcing


  real(dp) function H_for_thermal_balance(x)
    real(dp), intent(in) :: x
    real(dp) :: rho0, T0, dT0, ddT0, v01, dv01, B01
    real(dp) :: eps, deps
    real(dp) :: Kp, dKp, tcperp, dtcperpdr

    rho0 = background%density%rho0(x)
    T0 = background%temperature%T0(x)
    dT0 = background%temperature%dT0(x)
    ddT0 = background%temperature%ddT0(x)
    v01 = background%velocity%v01(x)
    dv01 = background%velocity%dv01(x)
    B01 = background%magnetic%B01(x)
    eps = grid%get_eps(x)
    deps = grid%get_deps()
    Kp = conduction%tcprefactor(x)
    dKp = conduction%dtcprefactordr(x)
    tcperp = conduction%tcperp(x)
    dtcperpdr = conduction%dtcperpdr(x)

    H_for_thermal_balance = rho0 * cooling%lambdaT(x) + (1.0_dp / rho0) * ( &
      T0 * rho0 * (deps * v01 + eps * dv01) / eps &
      - B01**2 * (dKp * dT0 + Kp * ddT0) &
      - (1.0_dp / eps) * ( &
          deps * tcperp * dT0 &
          + eps * dtcperpdr * dT0 &
          + eps * tcperp * ddT0 &
        ) &
      + (1.0_dp / settings%physics%get_gamma_1()) * dT0 * rho0 * v01 &
    )
  end function H_for_thermal_balance


  subroutine delete(this)
    class(heatloss_t), intent(inout) :: this
    nullify(settings)
    nullify(background)
    nullify(cooling)
    nullify(heating)
    nullify(grid)
    nullify(conduction)
    nullify(this%L0)
    nullify(this%dLdT)
    nullify(this%dLdrho)
  end subroutine delete

end module mod_heatloss
