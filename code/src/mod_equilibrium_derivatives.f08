!
! MODULE: mod_equilibrium_derivatives
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module containing the derivatives of the equilibrium arrays.
!
module mod_equilibrium_derivatives
  use mod_global_variables
  implicit none

  public

  !! Derivatives of the equilibrium variable arrays.
  !! d_xxx_dr     means d(xxx)/dr
  !! d_xxx_yyy_dr means d(xxx/yyy)/dr

  !! Default derivatives
  !> Derivative rho0
  real(dp), allocatable       :: d_rho0_dr(:)
  !> Derivative r*B02
  real(dp), allocatable       :: d_rB02_dr(:)
  !> Derivative B02/r
  real(dp), allocatable       :: d_B02_r_dr(:)
  !> Derivative B03
  real(dp), allocatable       :: d_B03_dr(:)
  !> Derivative T0
  real(dp), allocatable       :: d_T0_dr(:)
  !> Derivative B02
  real(dp), allocatable       :: d_B02_dr(:)

  !! Flow derivatives
  !> Derivative r*v02
  real(dp), allocatable       :: d_rv02_dr(:)
  !> Derivative v03
  real(dp), allocatable       :: d_v03_dr(:)

  !! Thermal conduction derivatives
  !> Derivative d(kappa_perp)/d(rho)
  real(dp), allocatable       :: d_tc_perp_eq_drho(:)
  !> Derivative d(kappa_perp)/d(T)
  real(dp), allocatable       :: d_tc_perp_eq_dT(:)
  !> Derivative d(kappa_perp)/(dB^2)
  real(dp), allocatable       :: d_tc_perp_eq_dB2(:)

  !! Radiative cooling derivatives
  !> Derivative dL/d(T)
  real(dp), allocatable       :: d_L_dT(:)
  !> Derivative dL/d(rho)
  real(dp), allocatable       :: d_L_drho(:)

  !! Resistivity derivatives
  !> derivative d(eta)/d(T)
  real(dp), allocatable       :: d_eta_dT(:)
  !> Double derivative B03
  real(dp), allocatable       :: dd_B03_dr(:)
  !> Double derivative B02
  real(dp), allocatable       :: dd_B02_dr(:)

contains

  !> Allocates and initialises all equilibrium derivatives arrays.
  subroutine initialise_equilibrium_derivatives()
    allocate(d_rho0_dr(gauss_gridpts))
    allocate(d_rB02_dr(gauss_gridpts))
    allocate(d_B02_r_dr(gauss_gridpts))
    allocate(d_B03_dr(gauss_gridpts))
    allocate(d_T0_dr(gauss_gridpts))
    allocate(d_B02_dr(gauss_gridpts))

    allocate(d_rv02_dr(gauss_gridpts))
    allocate(d_v03_dr(gauss_gridpts))

    allocate(d_tc_perp_eq_drho(gauss_gridpts))
    allocate(d_tc_perp_eq_dT(gauss_gridpts))
    allocate(d_tc_perp_eq_dB2(gauss_gridpts))

    allocate(d_L_dT(gauss_gridpts))
    allocate(d_L_drho(gauss_gridpts))

    allocate(d_eta_dT(gauss_gridpts))
    allocate(dd_B03_dr(gauss_gridpts))
    allocate(dd_B02_dr(gauss_gridpts))


    d_rho0_dr  = 0.0d0
    d_rB02_dr  = 0.0d0
    d_B02_r_dr = 0.0d0
    d_B03_dr   = 0.0d0
    d_T0_dr    = 0.0d0
    d_B02_dr   = 0.0d0

    d_rv02_dr  = 0.0d0
    d_v03_dr   = 0.0d0

    d_tc_perp_eq_drho = 0.0d0
    d_tc_perp_eq_dT   = 0.0d0
    d_tc_perp_eq_dB2  = 0.0d0

    d_L_dT     = 0.0d0
    d_L_drho   = 0.0d0

    d_eta_dT   = 0.0d0
    dd_B03_dr  = 0.0d0
    dd_B02_dr  = 0.0d0

  end subroutine initialise_equilibrium_derivatives

  !> Calculates the derivative thermal conduction arrays,
  !! if conduction is enabled.
  !! @param[in] T0_eq   Equilibrium temperature
  !! @param[in] rho0_eq Equilibrium density
  !! @param[in] B0_eq   Equilibrium magnetic field
  subroutine set_conduction_derivatives(T0_eq, rho0_eq, B0_eq)
    use mod_thermal_conduction

    real(dp), intent(in)    :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts), &
                               B0_eq(gauss_gridpts)

    call get_dkappa_perp_drho(T0_eq, rho0_eq, B0_eq, d_tc_perp_eq_drho)
    call get_dkappa_perp_dT(T0_eq, rho0_eq, B0_eq, d_tc_perp_eq_dT)
    call get_dkappa_perp_dB2(T0_eq, rho0_eq, B0_eq, d_tc_perp_eq_dB2)
  end subroutine set_conduction_derivatives

  !> Calculates the derivative radiative cooling arrays, if cooling is enabled.
  !! @param[in] T0_eq   Equilibrium temperature
  !! @param[in] rho0_eq Equilibrium density
  subroutine set_cooling_derivatives(T0_eq, rho0_eq)
    use mod_radiative_cooling

    real(dp), intent(in) :: T0_eq(gauss_gridpts), rho0_eq(gauss_gridpts)
    real(dp)             :: d_lambda_dT(gauss_gridpts)

    !! dL/dT = rho0 * d_lambda_dT (where lambda(T) = cooling curve)
    call get_dLambdadT(T0_eq, d_lambda_dT)
    d_L_dT = rho0_eq * d_lambda_dT

    !! dL/drho = lambda(T)
    call get_Lambda(T0_eq, d_L_drho)

  end subroutine set_cooling_derivatives

  !> Calculates the derivative resistivity arrays, if resistivity is enabled.
  !! @param[in] T0_eq   Equilibrium temperature
  subroutine set_resistivity_derivatives(T0_eq)
    use mod_resistivity

    real(dp), intent(in)  :: T0_eq(gauss_gridpts)

    call get_deta_dT(T0_eq, d_eta_dT)

  end subroutine set_resistivity_derivatives


  !> Cleaning routine, deallocates all arrays in this module.
  subroutine equilibrium_derivatives_clean()
    if (allocated(d_rho0_dr)) then
      deallocate(d_rho0_dr)
    end if
    if (allocated(d_rB02_dr)) then
      deallocate(d_rB02_dr)
    end if
    if (allocated(d_B02_r_dr)) then
      deallocate(d_B02_r_dr)
    end if
    if (allocated(d_B03_dr)) then
      deallocate(d_B03_dr)
    end if
    if (allocated(d_T0_dr)) then
      deallocate(d_T0_dr)
    end if
    if (allocated(d_B02_dr)) then
      deallocate(d_B02_dr)
    end if

    if (allocated(d_rv02_dr)) then
      deallocate(d_rv02_dr)
    end if
    if (allocated(d_v03_dr)) then
      deallocate(d_v03_dr)
    end if

    if (allocated(d_tc_perp_eq_drho)) then
      deallocate(d_tc_perp_eq_drho)
    end if
    if (allocated(d_tc_perp_eq_dT)) then
      deallocate(d_tc_perp_eq_dT)
    end if
    if (allocated(d_tc_perp_eq_dB2)) then
      deallocate(d_tc_perp_eq_dB2)
    end if

    if (allocated(d_L_dT)) then
      deallocate(d_L_dT)
    end if
    if (allocated(d_L_drho)) then
      deallocate(d_L_drho)
    end if

    if (allocated(d_eta_dT)) then
      deallocate(d_eta_dT)
    end if
    if (allocated(dd_B03_dr)) then
      deallocate(dd_B03_dr)
    end if
    if (allocated(dd_B02_dr)) then
      deallocate(dd_B02_dr)
    end if
  end subroutine equilibrium_derivatives_clean

end module mod_equilibrium_derivatives
