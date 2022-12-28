module mod_matrix_shortcuts
  use mod_global_variables, only: dp, NaN
  use mod_equilibrium, only: B_field
  use mod_equilibrium_params, only: k2, k3
  use mod_grid, only: eps_grid, d_eps_grid_dr
  use mod_logging, only: logger, str
  implicit none

  private

  public :: get_F_operator
  public :: get_diffF_operator
  public :: get_G_operator
  public :: get_wv_operator
  public :: get_Kp_operator

contains

  !> Calculates the $$\boldsymbol{\mathcal{F}}$$ operator, given as
  !! $$
  !! \boldsymbol{\mathcal{F}} =
  !!      \left(\frac{k_2}{\varepsilon}B_{02} \pm k_3B_{03}\right)
  !! $$
  function get_F_operator(gauss_idx, which) result(Foperator)
    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> which operator to calculate, <tt>"plus"</tt> or <tt>"minus"</tt>
    character(len=*), intent(in)  :: which
    !> the F operator on return
    real(dp)  :: Foperator

    if (which == "plus") then
      Foperator = ( &
        k2 * B_field % B02(gauss_idx) / eps_grid(gauss_idx) &
        + k3 * B_field % B03(gauss_idx) &
      )
    else if (which == "minus") then
      Foperator = ( &
        k2 * B_field % B02(gauss_idx) / eps_grid(gauss_idx) &
        - k3 * B_field % B03(gauss_idx) &
      )
    else
      Foperator = NaN
      call logger%error("requesting invalid F-operator sign: " // trim(which))
    end if
  end function get_F_operator


  !> Calculates the derivative of the $$\boldsymbol{\mathcal{F}}$$ operator, given as
  !! $$
  !! \boldsymbol{\mathcal{F}}' = \left[\frac{k_2}{\varepsilon}\left(
  !!    B_{02}' - \frac{\varepsilon'}{\varepsilon}B_{02}
  !! \right) + k_3B_{02}\right]
  !! $$
  function get_diffF_operator(gauss_idx, which) result(dFoperator)
    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> which operator to calculate, <tt>"plus"</tt> or <tt>"minus"</tt>
    character(len=*), intent(in)  :: which
    !> the F operator on return
    real(dp)  :: dFoperator
    real(dp)  :: eps, deps

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)


    if (which == "plus") then
      dFoperator = ( &
        (k2 / eps) * ( &
          B_field % d_B02_dr(gauss_idx) - deps * B_field % B02(gauss_idx) / eps &
        ) + k3 * B_field % d_B03_dr(gauss_idx) &
      )
    else if (which == "minus") then
      dFoperator = ( &
        (k2 / eps) * ( &
          B_field % d_B02_dr(gauss_idx) - deps * B_field % B02(gauss_idx) / eps &
        ) - k3 * B_field % d_B03_dr(gauss_idx) &
      )
    else
      dFoperator = NaN
      call logger%error("requesting invalid dF-operator sign: " // trim(which))
    end if
  end function get_diffF_operator


  !> Calculates the $$\boldsymbol{\mathcal{G}}$$ operator, given as
  !! $$
  !! \boldsymbol{\mathcal{G}} =
  !!      \left(k_3 B_{02} \mp \frac{k_2}{\varepsilon}B_{03}\right)
  !! $$
  function get_G_operator(gauss_idx, which) result(Goperator)
    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> which operator to calculate, <tt>"plus"</tt> or <tt>"minus"</tt>
    character(len=*), intent(in)  :: which
    !> the G operator on return
    real(dp)  :: Goperator

    if (which == "minus") then
      Goperator = ( &
        k3 * B_field % B02(gauss_idx) &
        - k2 * B_field % B03(gauss_idx) / eps_grid(gauss_idx) &
      )
    else if (which == "plus") then
      Goperator = ( &
        k3 * B_field % B02(gauss_idx) &
        + k2 * B_field % B03(gauss_idx) / eps_grid(gauss_idx) &
      )
    else
      Goperator = NaN
      call logger%error("requesting invalid G-operator sign: " // trim(which))
    end if
  end function get_G_operator


  !> Calculates the wave vector operator $$\boldsymbol{\mathcal{K}}$$, given as
  !! $$
  !! \boldsymbol{\mathcal{K}} =
  !!      \left(\frac{k_2^2}{\varepsilon} + \varepsilon k_3^2\right)
  !! $$
  function get_wv_operator(gauss_idx) result(wvoperator)
    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> the wave vector operator on return
    real(dp)  :: wvoperator

    wvoperator = k2**2 / eps_grid(gauss_idx) + eps_grid(gauss_idx) * k3**2
  end function get_wv_operator


  !> Calculates the (modified) conduction prefactor, given as
  !! \boldsymbol{K_p^+} =
  !!      \left(\boldsymbol{K_p} + \frac{\partial \kappa_\perp}{\partial(B^2)}\right)
  !! $$
  !! $$ \boldsymbol{K_p^{++}} = \left(
  !!    \frac{\partial \kappa_\perp}{\partial(B^2)}
  !!    - \frac{B_{01}^2}{B_0^2}\boldsymbol{K_p^+}
  !! \right)
  !! $$
  function get_Kp_operator(gauss_idx, which) result(Kp_operator)
    use mod_equilibrium, only: B_field, kappa_field

    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> which operator to calculate, <tt>"+", "++"</tt>
    character(len=*), intent(in)  :: which
    !> the (modified) $K_p$ operator on return
    real(dp)  :: Kp_operator

    real(dp)  :: Kp, Kp_plus, Kp_plusplus

    Kp = kappa_field % prefactor(gauss_idx)
    Kp_plus = Kp + kappa_field % d_kappa_perp_dB2(gauss_idx)
    Kp_plusplus = ( &
      kappa_field % d_kappa_perp_dB2(gauss_idx) &
      - (B_field % B01**2 * Kp_plus / B_field % B0(gauss_idx)**2) &
    )

    if (which == "+") then
      Kp_operator = Kp_plus
    else if (which == "++") then
      Kp_operator = Kp_plusplus
    else
      Kp_operator = NaN
      call logger%error("requesting invalid Kp-operator: " // trim(which))
    end if
  end function get_Kp_operator
end module mod_matrix_shortcuts
