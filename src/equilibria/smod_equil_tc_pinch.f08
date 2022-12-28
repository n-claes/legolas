! =============================================================================
!> This submodule defines a steady Taylor-Couette flow in a cylindrical geometry
!! where a plasma is confined between two (rotating) coaxial cylinders
!! with an azimuthal magnetic field and constant resistivity.
!!
!! This equilibrium is taken from
!! _Shalybkov, Dima.
!! "Rotational stabilization of pinch instabilities in Taylor-Couette flow.",
!! Physical Review E 75, 047302 (2007)_.
!!
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : density (constant)
!! - <tt>alpha</tt> = 1 : rotational speed of the inner cylinder
!! - <tt>beta</tt> = 2 : rotational speed of the outer cylinder
!! - <tt>cte_B02</tt> = 1 : azimuthal magnetic field at inner cylinder
!! - <tt>tau</tt> = 10 : azimuthal magnetic field at outer cylinder
!!
!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_tc_pinch
submodule (mod_equilibrium) smod_equil_tc_pinch
  implicit none

contains

  module procedure tc_pinch_eq
    use mod_equilibrium_params, only: cte_rho0, cte_B02, alpha, beta, tau

    real(dp) :: r, h, A, B, A2, B2, Bc1, Bc2, Bc3, Tstart, Tend, Ta, Ha, Re, Pm
    real(dp) :: x_start, x_end
    real(dp) :: fixed_eta_value, viscosity_value
    integer :: i, gauss_gridpts

    call settings%physics%enable_flow()

    settings%grid%coaxial = .true.

    if (settings%equilibrium%use_defaults) then
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(1.0_dp, 2.0_dp)
      cte_rho0 = 1.0d3
      alpha = 1.0d-6
      beta = 1.5d-6
      cte_B02 = 1.0d-3
      tau = 4.0d-3

      k2 = 0.0d0
      k3 = 1.0d0

      call settings%physics%enable_resistivity(fixed_resistivity_value=1.0e-4_dp)
      call settings%physics%enable_viscosity(viscosity_value=1.0e-6_dp)
    end if
    x_start = settings%grid%get_grid_start()
    x_end = settings%grid%get_grid_end()
    gauss_gridpts = settings%grid%get_gauss_gridpts()
    call initialise_grid(settings)


    fixed_eta_value = settings%physics%resistivity%get_fixed_resistivity()
    viscosity_value = settings%physics%viscosity%get_viscosity_value()

    rho_field % rho0 = cte_rho0

    h = x_end - x_start
    A = (beta * x_end**2 - alpha * x_start**2) / (x_end**2 - x_start**2)
    B = (x_start * x_end)**2 * (alpha - beta) / (x_end**2 - x_start**2)
    A2 = (x_end * tau - x_start * cte_B02) / (x_end**2 - x_start**2)
    B2 = x_start * x_end * (x_end * cte_B02 - x_start * tau) / (x_end**2 - x_start**2)

    Bc1 = cte_rho0 * A**2 - A2**2
    Bc2 = cte_rho0 * A * B - A2 * B2
    Bc3 = cte_rho0 * B**2 - B2**2

    Tstart = ( 0.5d0 * Bc1 * x_start**2 &
              + 2.0d0 * Bc2 * log(x_start) - 0.5d0 * Bc3 / x_start**2 &
              - 0.5d0 * (A2 * x_start + B2 / x_start)**2 ) / cte_rho0
    Tend = ( 0.5d0 * Bc1 * x_end**2 &
              + 2.0d0 * Bc2 * log(x_end) - 0.5d0 * Bc3 / x_end**2 &
              - 0.5d0 * (A2 * x_end + B2 / x_end)**2 ) / cte_rho0

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      v_field % v02(i) = A * r + B / r
      v_field % d_v02_dr(i) = A - B / r**2
      v_field % dd_v02_dr(i) = 2.0d0 * B / r**3

      B_field % B02(i)         = A2 * r + B2 / r
      B_field % B0(i)          = B_field % B02(i)
      B_field % d_B02_dr(i)    = A2 - B2 / r**2
      eta_field % dd_B02_dr(i) = 2.0d0 * B2 / r**3

      if ((Tstart > 0) .and. (Tend > 0)) then
        T_field % T0(i) = ( 0.5d0 * Bc1 * r**2 &
                          + 2.0d0 * Bc2 * log(r) - 0.5d0 * Bc3 / r**2 &
                          - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
      else
        T_field % T0(i) = 2.0d0 * max(abs(Tstart), abs(Tend)) + ( 0.5d0 * Bc1 * r**2 &
                          + 2.0d0 * Bc2 * log(r) - 0.5d0 * Bc3 / r**2 &
                          - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
      end if

      T_field % d_T0_dr(i) = ( (cte_rho0 * (v_field % v02(i))**2 &
                                - (B_field % B02(i))**2) / r &
                                - (B_field % B02(i)) * (B_field % d_B02_dr(i)) &
                              ) / cte_rho0
    end do

    Ta = (cte_rho0 * (v_field % v02(int(gauss_gridpts/2))) * h / viscosity_value)**2 &
          * 2.0d0 * h / (x_start + x_end)
    call logger%info("Taylor number:  " // str(Ta, fmt=exp_fmt))
    Pm = viscosity_value / (cte_rho0 * fixed_eta_value)
    call logger%info("Prandtl number: " // str(Pm, fmt=exp_fmt))
    Ha = cte_B02 * sqrt(x_start * (x_end - x_start)) &
          / sqrt(viscosity_value * fixed_eta_value)
    call logger%info("Hartmann number " // str(Ha, fmt=exp_fmt))
    Re = alpha * cte_rho0 * x_start * (x_end - x_start) / viscosity_value
    call logger%info("Reynolds number " // str(Re, fmt=exp_fmt))

  end procedure tc_pinch_eq

end submodule smod_equil_tc_pinch
