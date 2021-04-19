! =============================================================================
!> This submodule defines a steady Taylor-Couette flow in a cylindrical geometry
!! where the plasma is confined between two (rotating) coaxial cylinders
!! (without a magnetic field).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : density (constant)
!! - <tt>alpha</tt> = 1 : rotational speed (cycles per second) of the inner cylinder
!! - <tt>beta</tt> = 2 : rotational speed (cycles per second) of the outer cylinder
!! - <tt>eq_bool</tt> = <tt>False</tt> : to include magnetic field, set to <tt>True</tt>
!! - <tt>cte_B02</tt> = 0 : azimuthal magnetic field at inner cylinder, undefined if <tt>eq_bool</tt> = <tt>False</tt>
!! - <tt>tau</tt> = 0 : azimuthal magnetic field at outer cylinder, undefined if <tt>eq_bool</tt> = <tt>False</tt>
!!
!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_taylor_couette
submodule (mod_equilibrium) smod_equil_taylor_couette
  implicit none

contains

  module subroutine taylor_couette_eq()
    use mod_global_variables, only: coaxial, dp_LIMIT
    use mod_equilibrium_params, only: cte_rho0, cte_B02, alpha, beta, tau, eq_bool
    use mod_global_variables, only: viscosity_value

    real(dp)    :: r, h, Rrat, A, B, Ta, Tstart, A2, B2, K1, K2, K3
    integer     :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=1.0d0, default_x_end=2.0d0)
    call initialise_grid()

    coaxial = .true.

    if (use_defaults) then
      cte_rho0 = 1.0d0
      alpha = 1.0d0
      beta = 2.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      viscosity = .true.
      viscosity_value = 1.0d-3

      eq_bool = .false.
      cte_B02 = 0.0d0
      tau = 0.0d0
    end if

    rho_field % rho0 = cte_rho0
    h = x_end - x_start

    if ((abs(alpha) > dp_LIMIT) .or. (abs(beta) > dp_LIMIT)) then
      flow = .true.
      Rrat = x_start / x_end
      A = (alpha * Rrat**2 - beta) / (Rrat**2 - 1.0d0)
      B = x_start**2 * (alpha - beta) / (1.0d0 - Rrat**2)

      if (eq_bool) then
        A2 = (x_end * tau - x_start * cte_B02) / (x_end**2 - x_start**2)
        B2 = x_start * (cte_B02 - x_start * tau / x_end) / (1.0d0 - (x_start / x_end)**2)

        K1 = cte_rho0 * A**2 - A2**2
        K2 = cte_rho0 * A * B - A2 * B2
        K3 = cte_rho0 * B**2 - B2**2

        Tstart = ( 0.5d0 * K1 * x_start**2 &
                  + 2.0d0 * K2 * log(x_start) - 0.5d0 * K3 / x_start**2 &
                  - 0.5d0 * (A2 * x_start + B2 / x_start)**2 ) / cte_rho0
      else
        Tstart = 0.5d0 * ((A * x_start)**2 + 4.0d0 * A * B * log(x_start) - (B / x_start)**2)
      end if

      do i = 1, gauss_gridpts
        r = grid_gauss(i)

        v_field % v02(i) = A * r + B / r
        v_field % d_v02_dr(i) = A - B / r**2
        v_field % dd_v02_dr(i) = 2.0d0 * B / r**3

        if (eq_bool) then
          B_field % B02(i)         = A2 * r + B2 / r
          B_field % d_B02_dr(i)    = A2 - B2 / r**2
          eta_field % dd_B02_dr(i) = 2.0d0 * B2 / r**3

          if (Tstart > 0) then
            T_field % T0(i) = ( 0.5d0 * K1 * r**2 &
                              + 2.0d0 * K2 * log(r) - 0.5d0 * K3 / r**2 &
                              - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
          else
            T_field % T0(i) = 2.0d0 * abs(Tstart) + ( 0.5d0 * K1 * r**2 &
                              + 2.0d0 * K2 * log(r) - 0.5d0 * K3 / r**2 &
                              - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
          end if

          T_field % d_T0_dr(i) = ( (cte_rho0 * (v_field % v02(i))**2 &
                                    - (B_field % B02(i))**2) / r &
                                    - (B_field % B02(i)) * (B_field % d_B02_dr(i)) &
                                  ) / cte_rho0
        else
          if (Tstart > 0) then
            T_field % T0(i) = 0.5d0 * ((A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2)
          else
            T_field % T0(i) = 2.0d0 * abs(Tstart) + 0.5d0 * ((A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2)
          end if

          T_field % d_T0_dr(i) = (v_field % v02(i))**2 / r
        end if
      end do
    else if (eq_bool) then
      A2 = (x_end * tau - x_start * cte_B02) / (x_end**2 - x_start**2)
      B2 = x_start * (cte_B02 - x_start * tau / x_end) / (1.0d0 - (x_start / x_end)**2)

      K1 = - A2**2
      K2 = - A2 * B2
      K3 = - B2**2

      Tstart = ( 0.5d0 * K1 * x_start**2 &
                + 2.0d0 * K2 * log(x_start) - 0.5d0 * K3 / x_start**2 &
                - 0.5d0 * (A2 * x_start + B2 / x_start)**2 ) / cte_rho0

      do i = 1, gauss_gridpts
        r = grid_gauss(i)

        B_field % B02(i)         = A2 * r + B2 / r
        B_field % d_B02_dr(i)    = A2 - B2 / r**2
        eta_field % dd_B02_dr(i) = 2.0d0 * B2 / r**3

        if (Tstart > 0) then
          T_field % T0(i) = ( 0.5d0 * K1 * r**2 &
                            + 2.0d0 * K2 * log(r) - 0.5d0 * K3 / r**2 &
                            - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
        else
          T_field % T0(i) = 2.0d0 * abs(Tstart) + ( 0.5d0 * K1 * r**2 &
                            + 2.0d0 * K2 * log(r) - 0.5d0 * K3 / r**2 &
                            - 0.5d0 * (A2 * r + B2 / r)**2 ) / cte_rho0
        end if

        T_field % d_T0_dr(i) = ((- (B_field % B02(i))**2) / r &
                                  - (B_field % B02(i)) * (B_field % d_B02_dr(i)) &
                                ) / cte_rho0
      end do
    else
      T_field % T0 = 1.0d0
    end if

    Ta = (cte_rho0 * (v_field % v02(int(gauss_gridpts/2))) * h / viscosity_value)**2 &
          * 2.0d0 * h / (x_start + x_end)
    call log_message('Taylor number is ' // str(int(Ta)), level='info')

  end subroutine taylor_couette_eq

end submodule smod_equil_taylor_couette
