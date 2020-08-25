! =============================================================================
!> Module to inspect if certain conditions are fulfilled by doing
!! additional sanity checks on the equilibrium configuration.
!! For cylindrical geometries we check if \(k_2\) is an integer and if the
!! on-axis values obey regularity conditions. Equilibrium balance
!! for both the Cartesian and cylindrical cases is checked.
module mod_inspections
  use mod_global_variables, only: dp
  use mod_types, only: density_type, temperature_type, bfield_type, velocity_type, &
                       gravity_type, cooling_type, conduction_type
  use mod_logging, only: log_message, dp_fmt, exp_fmt, int_fmt, char_log
  implicit none

  private

  public :: perform_sanity_checks
  public :: handle_spurious_eigenvalues

contains


  !> General routine to do sanity checks on the different equilibrium types.
  !! We check the wavenumbers and on-axis values, as well as standard
  !! and non-adiabatic equilibrium force balance.
  subroutine perform_sanity_checks(rho_field, T_field, B_field, v_field, grav_field, rc_field, kappa_field)
    !> the type containing the density attributes
    type(density_type), intent(in)      :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)       :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)     :: v_field
    !> the type containing the gravity attributes
    type(gravity_type), intent(in)      :: grav_field
    !> the type containing the radiative cooling attributes
    type(cooling_type), intent(in)      :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in)   :: kappa_field

    call check_wavenumbers()
    call check_on_axis_values(B_field, v_field)
    call standard_equil_conditions(rho_field, T_field, B_field, v_field, grav_field)
    call nonadiab_equil_conditions(rho_field, T_field, rc_field, kappa_field)
  end subroutine perform_sanity_checks


  !> Handles spurious eigenvalue through removal.
  !! If requested, this can remove spurious eigenvalues on the edges
  !! of the spectrum. This usually only occurs in cylindrical geometries
  !! with forced on-axis conditions. The amount of eigenvalues to handle on every side
  !! of the imaginary axis is specified in the parfile.
  !! Example: <tt>nb_spurious_eigenvalues = 1</tt> removes the outermost
  !! eigenvalue on each side of the imaginary axis (so two in total).
  !! @note We don't actually _remove_ the spurious eigenvalues, but replace their value
  !!       with a large number so they can be filtered out during post-processing. @endnote
  !! @warning This routine should **ONLY** be used if **ABSOLUTELY** necessary. A better
  !!          on-axis treatment (e.g. <tt>r = 0.025</tt> instead of <tt>r = 0</tt>)
  !!          usually does a better job. A warning will be thrown if eigenvalues are removed.
  subroutine handle_spurious_eigenvalues(eigenvalues)
    use mod_global_variables, only: matrix_gridpts, remove_spurious_eigenvalues, nb_spurious_eigenvalues

    !> the eigenvalues with spurious modes replaced on exit
    complex(dp), intent(inout)  :: eigenvalues(matrix_gridpts)
    integer                     :: i, idx
    complex(dp)                 :: replacement

    if (.not. remove_spurious_eigenvalues) then
      return
    end if

    call log_message("handling spurious eigenvalues", level='debug')

    ! For now, the largest real eigenvalues are set to a large number so they
    ! do not appear on the plots.
    ! Do NOT sort the eigenvalues, otherwise the order is messed up for the eigenfunctions
    replacement = (1.0d20, 0.0d0)

    do i = 1, nb_spurious_eigenvalues
      ! handle real values, take large values from boundaries into account
      idx = maxloc(real(eigenvalues), dim=1, mask=(real(eigenvalues) < 1.0d15))
      eigenvalues(idx) = replacement
      idx = minloc(real(eigenvalues), dim=1, mask=(real(eigenvalues) < 1.0d15))
      eigenvalues(idx) = replacement
    end do

    write(char_log, int_fmt) nb_spurious_eigenvalues
    call log_message("spurious eigenvalues removed on every side: " // adjustl(char_log), level='warning')
  end subroutine handle_spurious_eigenvalues


  !> Sanity check on the wavenumbers.
  !! Checks if k2 is an integer in cylindrical geometry.
  !! @warning An error if thrown if the geometry is cylindrical and k2 is not
  !!          an integer.
  subroutine check_wavenumbers()
    use mod_global_variables, only: geometry, dp_LIMIT
    use mod_equilibrium_params, only: k2

    integer   :: k2_int

    k2_int = int(k2)

    if (geometry == 'cylindrical') then
      ! in cylindrical geometry k2 should be an integer
      if (abs(k2_int - k2) > dp_LIMIT) then
        write(char_log, dp_fmt) k2
        call log_message("cylindrical geometry but k2 is not an integer! Value: " // trim(char_log), level='error')
      end if
    end if
  end subroutine check_wavenumbers


  !> Checks if on-axis regularity conditions are satisfied in cylindrical geometry.
  !! We check if \(B_\theta, B_z', v_\theta\) and \(v_z'\) are smaller than <tt>1e-3</tt> on-axis.
  !! Nothing is checked if the geometry is Cartesian.
  !! @warning Throws a warning if:
  !!
  !! - \(B_\theta\) is not zero on-axis.
  !! - \(B_z'\) is not zero on-axis.
  !! - \(v_\theta\) is not zero on-axis.
  !! - \(v_z'\) is not zero on-axis. @endwarning
  subroutine check_on_axis_values(B_field, v_field)
    use mod_global_variables, only: geometry, x_start

    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)   :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in) :: v_field

    real(dp)  :: on_axis_limit

    if (geometry == 'Cartesian') then
      return
    end if

    on_axis_limit = 1.0d-3
    if (x_start > on_axis_limit) then
      return
    end if

    if (abs(B_field % B02(1)) > on_axis_limit) then
      write(char_log, exp_fmt) B_field % B02(1)
      call log_message("B_theta non-zero on axis! Value: " // trim(char_log), level='warning')
    end if
    if (abs(B_field % d_B03_dr(1)) > on_axis_limit) then
      write(char_log, exp_fmt) B_field % d_B03_dr(1)
      call log_message("dBz/dr non-zero on axis! Value: " // trim(char_log), level='warning')
    end if
    if (abs(v_field % v02(1)) > on_axis_limit) then
      write(char_log, exp_fmt) v_field % v02(1)
      call log_message("v_theta non-zero on axis! Value: " // trim(char_log), level='warning')
    end if
    if (abs(v_field % d_v03_dr(1)) > on_axis_limit) then
      write(char_log, exp_fmt) v_field % d_v03_dr(1)
      call log_message("dvz_dr non-zero on axis! Value: " // trim(char_log), level='warning')
    end if
  end subroutine check_on_axis_values


  !> Checks the standard force-balance equation for the equilibrium state. This is given by
  !! $$ \Bigl(p_0 + \frac{1}{2}B_0^2\Bigr)' + \rho_0 g
  !!    - \frac{\varepsilon'}{\varepsilon}\bigl(\rho_0 v_{02}^2 - B_{02}^2\bigr) $$
  !! and should be fulfilled.
  !! @warning   Throws a warning if force-balance is not satisfied.
  subroutine standard_equil_conditions(rho_field, T_field, B_field, v_field, grav_field)
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    !> the type containing the density attributes
    type(density_type), intent(in)      :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)       :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)     :: v_field
    !> the type containing the gravity attributes
    type(gravity_type), intent(in)      :: grav_field

    real(dp)  :: rho, drho, B02, dB02, B03, dB03, T0, dT0, grav, v02, v03
    real(dp)  :: eps, d_eps, r, discrepancy
    real(dp)  :: eq_cond(gauss_gridpts)
    integer   :: i, counter
    logical   :: satisfied

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, gauss_gridpts
      rho = rho_field % rho0(i)
      drho = rho_field % d_rho0_dr(i)
      B02 = B_field % B02(i)
      B03 = B_field % B03(i)
      dB02 = B_field % d_B02_dr(i)
      dB03 = B_field % d_B03_dr(i)
      T0 = T_field % T0(i)
      dT0 = T_field % d_T0_dr(i)
      grav = grav_field % grav(i)
      v02 = v_field % v02(i)
      v03 = v_field % v03(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i) = drho * T0 + rho * dT0 + B02 * dB02 + B03 * dB03 + rho * grav - (d_eps/eps) * (rho * v02**2 - B02**2)
      if (abs(eq_cond(i)) > dp_LIMIT) then
        counter = counter + 1
        satisfied = .false.
        if (abs(eq_cond(i)) > discrepancy) then
          discrepancy = abs(eq_cond(i))
          r = grid_gauss(i)
        end if
      end if
    end do

    if (.not. satisfied) then
      call log_message("standard equilibrium conditions not satisfied!", level='warning')
      write(char_log, dp_fmt) r
      call log_message("location of largest discrepancy: x = " // adjustl(trim(char_log)), level='warning')
      write(char_log, exp_fmt) discrepancy
      call log_message("value of largest discrepancy: " // adjustl(trim(char_log)), level='warning')
      write(char_log, int_fmt) counter
      call log_message("amount of nodes not satisfying criterion: " // adjustl(trim(char_log)), level='warning')
    end if
  end subroutine standard_equil_conditions


  !> Checks the non-adiabatic force-balance equation for the equilibrium state.
  !! This is given by
  !! $$ \frac{1}{\varepsilon}\bigl(\varepsilon \kappa_\bot T_0'\bigr)'
  !!    - \rho_0\mathscr{L}_0 = 0 $$
  !! The second derivative of the equilibrium temperature is evaluated numerically
  !! and does not have to be explicitly specified.
  !! @warning   Throws a warning if force-balance is not satisfied.
  subroutine nonadiab_equil_conditions(rho_field, T_field, rc_field, kappa_field)
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    !> the type containing the density attributes
    type(density_type), intent(in)      :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the radiative cooling attributes
    type(cooling_type), intent(in)      :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in)   :: kappa_field

    real(dp)  :: rho, dT0, ddT0, L0, kperp, dkperpdT
    real(dp)  :: eps, d_eps, r, discrepancy
    real(dp)  :: eq_cond(gauss_gridpts)
    integer   :: i, counter
    logical   :: satisfied

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, gauss_gridpts-1
      rho = rho_field % rho0(i)
      dT0 = T_field % d_T0_dr(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)
      kperp = kappa_field % kappa_perp(i)
      dkperpdT = kappa_field % d_kappa_perp_dT(i)
      L0 = rc_field % heat_loss(i)

      ! Do numerical differentiation for second T0 derivative, as it is only used here.
      ! This prevents having to calculate it every time in the submodules, 'approximately' equal here is fine.
      ddT0 = (T_field % d_T0_dr(i + 1) - T_field % d_T0_dr(i)) / (grid_gauss(i + 1) - grid_gauss(i))

      eq_cond(i) = d_eps / eps * kperp * dT0 + dkperpdT * dT0**2 + kperp * ddT0 - rho * L0
      if (abs(eq_cond(i)) > dp_LIMIT) then
        counter = counter + 1
        satisfied = .false.
        if (abs(eq_cond(i)) > discrepancy) then
          discrepancy = abs(eq_cond(i))
          r = grid_gauss(i)
        end if
      end if
    end do

    if (.not. satisfied) then
      call log_message("non-adiabatic equilibrium conditions not satisfied!", level='warning')
      write(char_log, dp_fmt) r
      call log_message("location of largest discrepancy: x = " // adjustl(trim(char_log)), level='warning')
      write(char_log, exp_fmt) discrepancy
      call log_message("value of largest discrepancy: " // adjustl(trim(char_log)), level='warning')
      write(char_log, int_fmt) counter
      call log_message("amount of nodes not satisfying criterion: " // adjustl(trim(char_log)), level='warning')
    end if
  end subroutine nonadiab_equil_conditions

end module mod_inspections
