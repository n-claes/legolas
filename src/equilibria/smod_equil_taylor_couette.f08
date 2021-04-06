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
!!
!! and can all be changed in the parfile. @endnote
!! @note The implementation does not allow for an immobile inner cylinder
!! (<tt>alpha</tt> = 0). @endnote
! SUBMODULE: smod_equil_taylor_couette
submodule (mod_equilibrium) smod_equil_taylor_couette
  implicit none

contains

  module subroutine taylor_couette_eq()
    use mod_global_variables, only: coaxial, dp_LIMIT
    use mod_equilibrium_params, only: cte_rho0, alpha, beta
    use mod_global_variables, only: viscosity_value

    real(dp)    :: r, h, Rrat, Orat, A, B, Ta, Tstart
    integer     :: i

    call allow_geometry_override(default_geometry='cylindrical', default_x_start=1.0d0, default_x_end=2.0d0)
    call initialise_grid()

    flow = .true.
    coaxial = .true.

    if (use_defaults) then
      cte_rho0 = 1.0d0
      alpha = 1.0d0
      beta = 2.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      viscosity = .true.
      viscosity_value = 1.0d-3
    end if

    h = x_end - x_start
    Rrat = x_start / x_end
    if (alpha > dp_LIMIT) then
      Orat = beta / alpha
      A = alpha * (Rrat**2 - Orat) / (Rrat**2 - 1.0d0)
      B = alpha * x_start**2 * (1.0d0 - Orat) / (1.0d0 - Rrat**2)
    else
      A = beta / (1.0d0 - Rrat**2)
      B = -beta * x_start**2 / (1.0d0 - Rrat**2)
    end if

    Tstart = 0.5d0 * ((A * x_start)**2 + 4.0d0 * A * B * log(x_start) - (B / x_start)**2)
    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      v_field % v02(i)    = A * r + B / r
      if (Tstart > 0) then
        T_field % T0(i)     = 0.5d0 * ((A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2)
      else
        T_field % T0(i) = 2.0d0*abs(Tstart) + 0.5d0 * ((A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2)
      end if

      v_field % d_v02_dr(i) = A - B / r**2
      T_field % d_T0_dr(i)  = (v_field % v02(i))**2 / r
    end do

    Ta = (cte_rho0 * (v_field % v02(int(gauss_gridpts/2))) * h / viscosity_value)**2 &
          * 2.0d0 * h / (x_start + x_end)
    call log_message('Taylor number is ' // str(int(Ta)), level='info')

  end subroutine taylor_couette_eq

end submodule smod_equil_taylor_couette
