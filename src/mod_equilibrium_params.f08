module mod_equilibrium_params
  use mod_global_variables, only: dp
  implicit none

  public

  !> Wavenumber in y-direction (Cartesian) or theta-direction (cylindrical)
  real(dp)  :: k2
  !> Wavenumber in z-direction
  real(dp)  :: k3

  !! General equilibrium values, either values themselves or prefactors
  real(dp)  :: cte_rho0
  real(dp)  :: cte_T0
  real(dp)  :: cte_B02
  real(dp)  :: cte_B03
  real(dp)  :: cte_v02
  real(dp)  :: cte_v03
  real(dp)  :: cte_p0

  !! General constants, used in various equilibrium configurations
  real(dp)  :: p1, p2, p3, p4, p5, p6, p7, p8
  real(dp)  :: alpha, beta, delta, theta, tau, lambda, nu
  real(dp)  :: r0, rc, rj
  real(dp)  :: Bth0, Bz0, V, j0, g

contains

  subroutine init_equilibrium_params()
    use mod_global_variables, only: NaN

    k2 = NaN
    k3 = NaN

    cte_rho0 = NaN
    cte_T0 = NaN
    cte_B02 = NaN
    cte_B03 = NaN
    cte_v02 = NaN
    cte_v03 = NaN
    cte_p0 = NaN

    p1 = NaN
    p2 = NaN
    p3 = NaN
    p4 = NaN
    p5 = NaN
    p6 = NaN
    p7 = NaN
    p8 = NaN

    alpha = NaN
    beta = NaN
    delta = NaN
    theta = NaN
    tau = NaN
    lambda = NaN
    nu = NaN

    r0 = NaN
    rc = NaN
    rj = NaN
    Bth0 = NaN
    Bz0 = NaN
    V = NaN
    j0 = NaN
    g = NaN
  end subroutine init_equilibrium_params

end module mod_equilibrium_params
