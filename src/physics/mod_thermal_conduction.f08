! =============================================================================
!> This module is responsible for calculating and setting the
!! thermal conduction values based on the equilibrium configuration.
module mod_thermal_conduction
  use mod_global_variables, only: dp
  use mod_physics_utils, only: physics_i
  use mod_physical_constants, only: dpi, coulomb_log, tc_pf_kappa_para, tc_pf_kappa_perp
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  implicit none

  private

  type, public :: conduction_t
    procedure(physics_i), pointer, nopass :: tcpara
    procedure(physics_i), pointer, nopass :: dtcparadT
    procedure(physics_i), pointer, nopass :: dtcparadr

    procedure(physics_i), pointer, nopass :: tcperp
    procedure(physics_i), pointer, nopass :: dtcperpdrho
    procedure(physics_i), pointer, nopass :: dtcperpdT
    procedure(physics_i), pointer, nopass :: dtcperpdB2
    procedure(physics_i), pointer, nopass :: dtcperpdr
  contains
    procedure, public :: tcprefactor
    procedure, public :: dtcprefactordr
    procedure, public :: delete
  end type conduction_t

  public :: new_conduction

contains


  function new_conduction() result(conduction)
    type(conduction_t) :: conduction

    conduction%tcpara => get_tcpara
    conduction%dtcparadT => get_dtcparadT
    conduction%dtcparadr => get_dtcparadr

    conduction%tcperp => get_tcperp
    conduction%dtcperpdrho => get_dtcperpdrho
    conduction%dtcperpdT => get_dtcperpdT
    conduction%dtcperpdB2 => get_dtcperpdB2
    conduction%dtcperpdr => get_dtcperpdr
  end function new_conduction


  real(dp) function get_tcpara(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: T0

    get_tcpara = 0.0_dp
    if (.not. settings%physics%conduction%has_parallel_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_para()) then
      get_tcpara = settings%physics%conduction%get_fixed_tc_para()
      return
    end if

    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    get_tcpara = ( &
      tc_pf_kappa_para * T0**2.5_dp / coulomb_log &
    ) / settings%units%get_unit_conduction()
  end function get_tcpara


  real(dp) function get_dtcparadT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: unit_temperature, unit_conduction
    real(dp) :: T0

    get_dtcparadT = 0.0_dp
    if (.not. settings%physics%conduction%has_parallel_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_para()) return

    unit_temperature = settings%units%get_unit_temperature()
    unit_conduction = settings%units%get_unit_conduction()
    T0 = background%temperature%T0(x) * unit_temperature
    get_dtcparadT = ( &
      tc_pf_kappa_para * 2.5_dp * T0**1.5_dp / coulomb_log &
    ) / (unit_conduction / unit_temperature)
  end function get_dtcparadT


  real(dp) function get_dtcparadr(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: dT0

    get_dtcparadr = 0.0_dp
    if (.not. settings%physics%conduction%has_parallel_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_para()) return

    ! position derivative is dtcpara/dT * T0'
    dT0 = background%temperature%dT0(x)
    get_dtcparadr = get_dtcparadT(x, settings, background) * dT0
  end function get_dtcparadr



  real(dp) function get_tcperp(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: T0, nh0, B0

    get_tcperp = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) then
      get_tcperp = settings%physics%conduction%get_fixed_tc_perp()
      return
    end if
    if (.not. settings%has_bfield()) then
      get_tcperp = get_tcpara(x, settings, background)
      return
    end if

    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    nh0 = background%density%rho0(x) * settings%units%get_unit_numberdensity()
    B0 = background%magnetic%get_B0(x) * settings%units%get_unit_magneticfield()
    get_tcperp = ( &
      tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nh0**2 / (B0**2 * sqrt(T0)) &
    ) / settings%units%get_unit_conduction()
  end function get_tcperp


  real(dp) function get_dtcperpdrho(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: T0, nh0, B0

    get_dtcperpdrho = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return
    if (.not. settings%has_bfield()) return

    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    nh0 = background%density%rho0(x) * settings%units%get_unit_numberdensity()
    B0 = background%magnetic%get_B0(x) * settings%units%get_unit_magneticfield()
    get_dtcperpdrho = ( &
      2.0_dp * tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nh0 &
      / (B0**2 * sqrt(T0)) &
    ) / (settings%units%get_unit_conduction() / settings%units%get_unit_numberdensity())
  end function get_dtcperpdrho


  real(dp) function get_dtcperpdT(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: T0, nh0, B0

    get_dtcperpdT = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return
    if (.not. settings%has_bfield()) then
      get_dtcperpdT = get_dtcparadT(x, settings, background)
      return
    end if

    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    nh0 = background%density%rho0(x) * settings%units%get_unit_numberdensity()
    B0 = background%magnetic%get_B0(x) * settings%units%get_unit_magneticfield()
    get_dtcperpdT = ( &
      -0.5_dp * tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nh0**2 &
      / (B0**2 * T0**(3.0_dp / 2.0_dp)) &
    ) / (settings%units%get_unit_conduction() / settings%units%get_unit_temperature())
  end function get_dtcperpdT


  real(dp) function get_dtcperpdB2(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: T0, nh0, B0
    real(dp) :: unit_magneticfield

    get_dtcperpdB2 = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return
    if (.not. settings%has_bfield()) return

    unit_magneticfield = settings%units%get_unit_magneticfield()
    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    nh0 = background%density%rho0(x) * settings%units%get_unit_numberdensity()
    B0 = background%magnetic%get_B0(x) * unit_magneticfield
    get_dtcperpdB2 = ( &
      -tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nh0**2 / (B0**4 * sqrt(T0)) &
    ) / (settings%units%get_unit_conduction() / unit_magneticfield**2)
  end function get_dtcperpdB2


  real(dp) function get_dtcperpdr(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: drho0, dT0, B0, dB0

    get_dtcperpdr = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return

    if (.not. settings%has_bfield()) then
      get_dtcperpdr = get_dtcparadr(x, settings, background)
      return
    end if
    drho0 = background%density%drho0(x)
    dT0 = background%temperature%dT0(x)
    B0 = background%magnetic%get_B0(x)
    dB0 = background%magnetic%get_dB0(x)
    get_dtcperpdr = ( &
      get_dtcperpdrho(x, settings, background) * drho0 &
      + get_dtcperpdT(x, settings, background) * dT0 &
      + get_dtcperpdB2(x, settings, background) * 2.0_dp * B0 * dB0 &
    )
  end function get_dtcperpdr


  real(dp) function tcprefactor(this, x, settings, background)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: B0

    tcprefactor = 0.0_dp
    if (.not. settings%physics%conduction%is_enabled()) return
    if (.not. settings%has_bfield()) return

    B0 = background%magnetic%get_B0(x)
    tcprefactor = ( &
      this%tcpara(x, settings, background) - this%tcperp(x, settings, background) &
    ) / B0**2
  end function tcprefactor


  real(dp) function dtcprefactordr(this, x, settings, background)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    real(dp) :: B0, dB0

    dtcprefactordr = 0.0_dp
    if (.not. settings%physics%conduction%is_enabled()) return
    if (.not. settings%has_bfield()) return

    B0 = background%magnetic%get_B0(x)
    dB0 = background%magnetic%get_dB0(x)
    dtcprefactordr = ( &
      ( &
        this%dtcparadr(x, settings, background) &
        - this%dtcperpdr(x, settings, background) &
      ) * B0 &
      - 2.0_dp * ( &
        this%tcpara(x, settings, background) &
        - this%tcperp(x, settings, background) &
      ) * dB0 &
    ) / B0**3
  end function dtcprefactordr


  pure subroutine delete(this)
    class(conduction_t), intent(inout) :: this
    nullify(this%tcpara)
    nullify(this%dtcparadT)
    nullify(this%dtcparadr)

    nullify(this%tcperp)
    nullify(this%dtcperpdrho)
    nullify(this%dtcperpdT)
    nullify(this%dtcperpdB2)
    nullify(this%dtcperpdr)
  end subroutine delete

end module mod_thermal_conduction
