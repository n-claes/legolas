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

  module subroutine tc_pinch_eq()
    use mod_global_variables, only: coaxial
    use mod_equilibrium_params, only: cte_rho0, cte_B02, alpha, beta, tau
    use mod_global_variables, only: viscosity_value, use_fixed_resistivity, fixed_eta_value

    real(dp)    :: r, h, Rrat, A, B, A2, B2, Bc1, Bc2, Bc3, Tstart, Tend, Ta, Ha, Re, Pm
    integer     :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=1.0d0, default_x_end=2.0d0)
    call initialise_grid()

    flow = .true.
    resistivity = .true.
    use_fixed_resistivity = .true.
    viscosity = .true.
    coaxial = .true.

    if (use_defaults) then
      cte_rho0 = 1.0d0
      alpha = 1.0d0
      beta = 2.0d0
      cte_B02 = 1.0d0
      tau = 10.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      viscosity_value = 1.0d-3
      fixed_eta_value = 1.0d-4
    end if

    rho_field % rho0 = cte_rho0

    h = x_end - x_start
    Rrat = x_start / x_end
    A = (alpha * Rrat**2 - beta) / (Rrat**2 - 1.0d0)
    B = x_start**2 * (alpha - beta) / (1.0d0 - Rrat**2)
    A2 = (x_end * tau - x_start * cte_B02) / (x_end**2 - x_start**2)
    B2 = x_start * (cte_B02 - x_start * tau / x_end) / (1.0d0 - (x_start / x_end)**2)

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
    call log_message('Taylor number:    ' // str(Ta, fmt='(e8.2)'), level='info')
    Pm = viscosity_value / (cte_rho0 * fixed_eta_value)
    call log_message('Prandtl number:   ' // str(Pm, fmt='(e8.2)'), level='info')
    Ha = cte_B02 * sqrt(x_start * (x_end - x_start)) &
          / sqrt(viscosity_value * fixed_eta_value)
    call log_message('Hartmann number:  ' // str(Ha, fmt='(e8.2)'), level='info')
    Re = alpha * cte_rho0 * x_start * (x_end - x_start) / viscosity_value
    call log_message('Reynolds number:  ' // str(Re, fmt='(e8.2)'), level='info')

  end subroutine tc_pinch_eq

end submodule smod_equil_tc_pinch
