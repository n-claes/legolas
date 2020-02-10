module mod_equilibrium_params
  use mod_global_variables, only: dp
  implicit none

  public

  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)  :: k2
  !> Wavenumber in z-direction
  real(dp)  :: k3

  !! General equilibrium values
  real(dp)  :: cte_rho0
  real(dp)  :: cte_T0
  real(dp)  :: cte_B02
  real(dp)  :: cte_B03
  real(dp)  :: cte_v02
  real(dp)  :: cte_v03
  real(dp)  :: cte_grav
  real(dp)  :: cte_p0

  !! General constants, used in various equilibrium configurations
  real(dp)  :: p1, p2, p3, p4, p5, p6, p7, p8
  real(dp)  :: alpha, beta, delta, theta, tau, lambda, nu
  real(dp)  :: r0, rc, rj
  real(dp)  :: Bth, Bz, V

contains

  subroutine init_equilibrium_params()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    real(dp)  :: nan

    nan = ieee_value(nan, ieee_quiet_nan)

    k2 = nan
    k3 = nan

    cte_rho0 = nan
    cte_T0 = nan
    cte_B02 = nan
    cte_B03 = nan
    cte_v02 = nan
    cte_v03 = nan
    cte_grav = nan
    cte_p0 = nan

    p1 = nan
    p2 = nan
    p3 = nan
    p4 = nan
    p5 = nan
    p6 = nan
    p7 = nan
    p8 = nan

    alpha = nan
    beta = nan
    delta = nan
    theta = nan
    tau = nan
    lambda = nan
    nu = nan

    r0 = nan
    rc = nan
    rj = nan
    Bth = nan
    Bz = nan
    V = nan
  end subroutine init_equilibrium_params

end module mod_equilibrium_params
