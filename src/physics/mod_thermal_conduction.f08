! =============================================================================
!> This module is responsible for calculating and setting the
!! thermal conduction values based on the equilibrium configuration.
module mod_thermal_conduction
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi, coulomb_log, tc_pf_kappa_para, tc_pf_kappa_perp
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  implicit none

  private

  type(settings_t), pointer :: settings
  type(background_t), pointer :: background

  type, public :: conduction_t
    procedure(real(dp)), pointer, nopass :: tcpara
    procedure(real(dp)), pointer, nopass :: dtcparadT

    procedure(real(dp)), pointer, nopass :: tcperp
    procedure(real(dp)), pointer, nopass :: dtcperpdrho
    procedure(real(dp)), pointer, nopass :: dtcperpdT
    procedure(real(dp)), pointer, nopass :: dtcperpdB2
  contains
    procedure, public :: get_dtcparadr
    procedure, public :: get_dtcperpdr
    procedure, public :: get_tcprefactor
    procedure, public :: get_dtcprefactordr
    procedure, public :: delete
  end type conduction_t

  public :: new_conduction

contains


  function new_conduction(settings_tgt, background_tgt) result(conduction)
    type(settings_t), target, intent(in) :: settings_tgt
    type(background_t), target, intent(in) :: background_tgt
    type(conduction_t) :: conduction

    settings => settings_tgt
    background => background_tgt

    conduction%tcpara => get_tcpara
    conduction%dtcparadT => get_dtcparadT

    conduction%tcperp => get_tcperp
    conduction%dtcperpdrho => get_dtcperpdrho
    conduction%dtcperpdT => get_dtcperpdT
    conduction%dtcperpdB2 => get_dtcperpdB2
  end function new_conduction


  real(dp) function get_tcpara(x)
    real(dp), intent(in) :: x
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


  real(dp) function get_dtcparadT(x)
    real(dp), intent(in) :: x
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


  impure elemental real(dp) function get_dtcparadr(this, x)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: dT0

    get_dtcparadr = 0.0_dp
    if (.not. settings%physics%conduction%has_parallel_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_para()) return

    ! position derivative is dtcpara/dT * T0'
    dT0 = background%temperature%dT0(x)
    get_dtcparadr = this%dtcparadT(x) * dT0
  end function get_dtcparadr



  real(dp) function get_tcperp(x)
    real(dp), intent(in) :: x
    real(dp) :: T0, nh0, B0

    get_tcperp = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) then
      get_tcperp = settings%physics%conduction%get_fixed_tc_perp()
      return
    end if
    if (.not. settings%has_bfield()) then
      get_tcperp = get_tcpara(x)
      return
    end if

    T0 = background%temperature%T0(x) * settings%units%get_unit_temperature()
    nh0 = background%density%rho0(x) * settings%units%get_unit_numberdensity()
    B0 = background%magnetic%get_B0(x) * settings%units%get_unit_magneticfield()
    get_tcperp = ( &
      tc_pf_kappa_para * tc_pf_kappa_perp * coulomb_log * nh0**2 / (B0**2 * sqrt(T0)) &
    ) / settings%units%get_unit_conduction()
  end function get_tcperp


  real(dp) function get_dtcperpdrho(x)
    real(dp), intent(in) :: x
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


  real(dp) function get_dtcperpdT(x)
    real(dp), intent(in) :: x
    real(dp) :: T0, nh0, B0

    get_dtcperpdT = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return
    if (.not. settings%has_bfield()) then
      get_dtcperpdT = get_dtcparadT(x)
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


  real(dp) function get_dtcperpdB2(x)
    real(dp), intent(in) :: x
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


  impure elemental real(dp) function get_dtcperpdr(this, x)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: drho0, dT0, B0, dB0

    get_dtcperpdr = 0.0_dp
    if (.not. settings%physics%conduction%has_perpendicular_conduction()) return
    if (settings%physics%conduction%has_fixed_tc_perp()) return

    if (.not. settings%has_bfield()) then
      get_dtcperpdr = this%get_dtcparadr(x)
      return
    end if
    drho0 = background%density%drho0(x)
    dT0 = background%temperature%dT0(x)
    B0 = background%magnetic%get_B0(x)
    dB0 = background%magnetic%get_dB0(x)
    get_dtcperpdr = ( &
      this%dtcperpdrho(x) * drho0 &
      + this%dtcperpdT(x) * dT0 &
      + this%dtcperpdB2(x) * 2.0_dp * B0 * dB0 &
    )
  end function get_dtcperpdr


  impure elemental real(dp) function get_tcprefactor(this, x)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: B0

    get_tcprefactor = 0.0_dp
    if (.not. settings%physics%conduction%is_enabled()) return
    if (.not. settings%has_bfield()) return

    B0 = background%magnetic%get_B0(x)
    get_tcprefactor = (this%tcpara(x) - this%tcperp(x)) / B0**2
  end function get_tcprefactor


  impure elemental real(dp) function get_dtcprefactordr(this, x)
    class(conduction_t), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: B0, dB0

    get_dtcprefactordr = 0.0_dp
    if (.not. settings%physics%conduction%is_enabled()) return
    if (.not. settings%has_bfield()) return

    B0 = background%magnetic%get_B0(x)
    dB0 = background%magnetic%get_dB0(x)
    get_dtcprefactordr = ( &
      (this%get_dtcparadr(x) - this%get_dtcperpdr(x)) * B0 &
      - 2.0_dp * (this%tcpara(x) - this%tcperp(x)) * dB0 &
    ) / B0**3
  end function get_dtcprefactordr


  subroutine delete(this)
    class(conduction_t), intent(inout) :: this

    nullify(settings)
    nullify(background)

    nullify(this%tcpara)
    nullify(this%dtcparadT)

    nullify(this%tcperp)
    nullify(this%dtcperpdrho)
    nullify(this%dtcperpdT)
    nullify(this%dtcperpdB2)
  end subroutine delete

end module mod_thermal_conduction
