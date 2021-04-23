module mod_matrix_shortcuts
  use mod_global_variables, only: dp
  use mod_equilibrium, only: B_field
  use mod_equilibrium_params, only: k2, k3
  use mod_grid, only: eps_grid, d_eps_grid_dr
  use mod_logging, only: log_message, str
  implicit none

  private

  public :: get_F_operator
  public :: get_G_operator
  public :: get_wv_operator

contains

  !> Calculates the $$\boldsymbol{\mathcal{F}}$$ operator, given as
  !! $$
  !! \boldsymbol{\mathcal{F}} =
  !!      \left(\frac{k_2}{\varepsilon}B_{03} \pm k_3B_{02}\right)
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
      call log_message( &
        "requesting invalid F-operator sign: " // trim(which), level="error" &
      )
    end if
  end function get_F_operator


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
      call log_message( &
        "requesting invalid G-operator sign: " // trim(which), level="error" &
      )
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

end module mod_matrix_shortcuts
