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
  use mod_logging, only: log_message, dp_fmt, exp_fmt, int_fmt, char_log, char_log2
  implicit none

  private

  public :: perform_sanity_checks
  public :: handle_spurious_eigenvalues

contains


  !> General routine to do sanity checks on the different equilibrium types.
  !! We check the wavenumbers and on-axis values, as well as standard
  !! and non-adiabatic equilibrium force balance.
  subroutine perform_sanity_checks( &
    rho_field, T_field, B_field, v_field, grav_field, rc_field, kappa_field &
  )
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
    type(cooling_type), intent(inout)   :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in)   :: kappa_field

    call check_wavenumbers()
    call check_on_axis_values(B_field, v_field)
    call standard_equil_conditions(rho_field, T_field, B_field, v_field, grav_field)
    call continuity_equil_conditions(rho_field, v_field)
    call induction_equil_conditions(B_field, v_field)
    ! set the energy balance based on the equilibrium conditions
    call set_energy_balance(rho_field, T_field, B_field, v_field, rc_field, kappa_field)
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
    call log_message( &
      "spurious eigenvalues removed on every side: " // adjustl(char_log), &
      level='warning' &
    )
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
        call log_message( &
          "cylindrical geometry but k2 is not an integer! Value: " &
            // adjustl(trim(char_log)), &
          level='error' &
        )
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
      call log_message( &
        "B_theta non-zero on axis! Value: " // trim(char_log), &
        level='warning' &
      )
    end if
    if (abs(B_field % d_B03_dr(1)) > on_axis_limit) then
      write(char_log, exp_fmt) B_field % d_B03_dr(1)
      call log_message( &
        "dBz/dr non-zero on axis! Value: " // trim(char_log), &
        level='warning' &
      )
    end if
    if (abs(v_field % v02(1)) > on_axis_limit) then
      write(char_log, exp_fmt) v_field % v02(1)
      call log_message( &
        "v_theta non-zero on axis! Value: " // trim(char_log), &
        level='warning' &
      )
    end if
    if (abs(v_field % d_v03_dr(1)) > on_axis_limit) then
      write(char_log, exp_fmt) v_field % d_v03_dr(1)
      call log_message( &
        "dvz_dr non-zero on axis! Value: " // trim(char_log), &
        level='warning' &
      )
    end if
  end subroutine check_on_axis_values


  !> Checks the standard force-balance equation for the equilibrium state. This
  !! results in three expressions,
  !! $$ \Bigl(p_0 + \frac{1}{2}B_0^2\Bigr)' + \rho_0 g + \rho_0 v_{01} v_{01}'
  !!    - \frac{\varepsilon'}{\varepsilon}\bigl(\rho_0 v_{02}^2 - B_{02}^2\bigr) = 0, $$
  !! $$ \rho_0 v_{01} \bigl( v_{02}' + \frac{\varepsilon'}{\varepsilon} v_{02} \bigr)
  !!    - \frac{B_{01}}{\varepsilon} (\varepsilon B_{02})' = 0, $$
  !! $$ \rho_0 v_{01} v_{03}' - B_{01} B_{03}' = 0, $$
  !! and they should all be fulfilled.
  !! @warning   Throws a warning if force-balance is not satisfied.
  subroutine standard_equil_conditions(rho_field, T_field, B_field, v_field, grav_field)
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT, geometry
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

    real(dp)  :: rho, drho, B01, B02, dB02, B03, dB03, T0, dT0, grav
    real(dp)  :: v01, v02, v03, dv01, dv02, dv03
    real(dp)  :: eps, d_eps, r(3), discrepancy(3)
    real(dp)  :: eq_cond(gauss_gridpts, 3)
    integer   :: i, j, counter(3)
    logical   :: satisfied(3)

    B01 = B_field % B01
    if ((geometry == 'cylindrical') .and. (abs(B01) > dp_LIMIT)) then
      call log_message('B01 component currently not supported for cylindrical &
                        &geometries !', level='error')
    end if

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
      v01 = v_field % v01(i)
      v02 = v_field % v02(i)
      v03 = v_field % v03(i)
      dv01 = v_field % d_v01_dr(i)
      dv02 = v_field % d_v02_dr(i)
      dv03 = v_field % d_v03_dr(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i, 1) = drho * T0 + rho * dT0 + B02 * dB02 + B03 * dB03 + rho * grav &
                    - (d_eps/eps) * (rho * v02**2 - B02**2) + rho * v01 * dv01
      eq_cond(i, 2) = rho * v01 * (dv02 + v02 * d_eps / eps) - B01 * (dB02 + B02 * d_eps / eps)
      eq_cond(i, 3) = rho * v01 * dv03 - B01 * dB03

      do j = 1, 3
        if (abs(eq_cond(i, j)) > dp_LIMIT) then
          counter(j) = counter(j) + 1
          satisfied(j) = .false.
          if (abs(eq_cond(i, j)) > discrepancy(j)) then
            discrepancy(j) = abs(eq_cond(i, j))
            r(j) = grid_gauss(i)
          end if
        end if
      end do
    end do

    do j = 1, 3
      if (satisfied(j)) then
        cycle
      end if
      call log_message( &
        "standard equilibrium conditions not satisfied!", &
        level="warning" &
      )
      write(char_log2, int_fmt) j
      write(char_log, dp_fmt) r(j)
      call log_message( &
        "location of largest discrepancy (" // trim(adjustl(char_log2)) // &
                                      "): x = " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, exp_fmt) discrepancy(j)
      call log_message( &
        "value of largest discrepancy (" // trim(adjustl(char_log2)) // &
                                      "): " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, int_fmt) counter(j)
      call log_message( &
        "amount of nodes not satisfying criterion (" // trim(adjustl(char_log2)) // &
                                          "): " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(*,*) ""
    end do
  end subroutine standard_equil_conditions


  !> Enforces the non-adiabatic force-balance equation for the equilibrium state.
  !! This is given by
  !! $$
  !! T_0 \rho_0 \frac{\left(\varepsilon v_{01}\right)'}{\varepsilon}
  !! + \rho_0 \mathscr{L}_0
  !! - B_{01}^2\left[\frac{\kappa_{\parallel,0} - \kappa_{\perp,0}}{B_0^2} T_0'\right]'
  !! - \frac{1}{\varepsilon}\left(\varepsilon \kappa_{\perp, 0} T_0'\right)'
  !! + \frac{1}{(\gamma - 1)}T_0'\rho_0 v_{01} = 0
  !! $$
  !! This subroutine essentially sets $\mathscr{L}_0$ in such a way that this equation
  !! is satisfied. If the heating is assumed to only depend on the equilibrium,
  !! and if there is no $B_{01}$, $v_{01}$ or perpendicular thermal conduction,
  !! then $\mathscr{L}_0 = 0$. If one (or more) of these effects are present, then
  !! $\mathscr{L}_0 = 0$ is no longer true.
  !!  The <tt>rc_field % heat_loss</tt> attribute is modified on exit.
  subroutine set_energy_balance( &
    rho_field, T_field, B_field, v_field, rc_field, kappa_field &
  )
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT, gamma_1
    use mod_grid, only: eps_grid, d_eps_grid_dr

    !> the type containing the density attributes
    type(density_type), intent(in)  :: rho_field
    !> the type containing the temperature attributes
    type(temperature_type), intent(in)  :: T_field
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in) :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in) :: v_field
    !> the type containing the radiative cooling attributes
    type(cooling_type), intent(inout)  :: rc_field
    !> the type containing the thermal conduction attributes
    type(conduction_type), intent(in) :: kappa_field

    real(dp)  :: rho, drho, T0, dT0, ddT0
    real(dp)  :: B01, B02, dB02, B03, dB03, B0, dB0
    real(dp)  :: v01, dv01
    real(dp)  :: kappa_perp, dkappa_perp_dr, Kp, dKp
    real(dp)  :: eps, deps
    integer   :: i

    B01 = B_field % B01
    do i = 1, gauss_gridpts
      rho = rho_field % rho0(i)
      drho = rho_field % d_rho0_dr(i)
      T0 = T_field % T0(i)
      dT0 = T_field % d_T0_dr(i)
      ddT0 = T_field % dd_T0_dr(i)
      B02 = B_field % B02(i)
      dB02 = B_field % d_B02_dr(i)
      B03 = B_field % B03(i)
      dB03 = B_field % d_B03_dr(i)
      B0 = B_field % B0(i)
      dB0 = (B02 * dB02 + B03 * dB03) / B0
      v01 = v_field % v01(i)
      dv01 = v_field % d_v01_dr(i)
      eps = eps_grid(i)
      deps = d_eps_grid_dr(i)
      kappa_perp = kappa_field % kappa_perp(i)
      dkappa_perp_dr = kappa_field % d_kappa_perp_dr(i)
      Kp = kappa_field % prefactor(i)
      dKp = kappa_field % d_prefactor_dr(i)

      ! set L0, this is equal to 0 if there is no B01, v01 or kappa_perp. The extra
      ! rho * lambda(T0) factor cancels out with the radiative cooling contribution
      rc_field % heat_loss(i) = ( &
        (T0 / eps) * (deps * v01 + eps * dv01) &
        + dT0 * v01 / gamma_1 &
        - (B01**2 / rho) * (dKp * dT0 + Kp * ddT0) &
        - ( &
          deps * kappa_perp * dT0 &
          + eps * dkappa_perp_dr * dT0 &
          + eps * kappa_perp * ddT0 &
        ) / (eps * rho) &
      )
    end do

    ! log this if it's set
    if (any(abs(rc_field % heat_loss) > dp_limit)) then
      call log_message( &
        "encountered non-zero B01, v01 or kappa_perp, energy balance has been set", &
        level="info" &
      )
    end if
    ! double check if L0=0 holds
    if ( &
      abs(B01) < dp_LIMIT .and. &
      all(abs(v_field % v01) < dp_LIMIT) .and. &
      all(abs(kappa_field % kappa_perp) < dp_LIMIT) &
    ) then
      ! if B01 = v01 = kappa_perp = 0, then L0 must be zero
      if (any(abs(rc_field % heat_loss) > dp_LIMIT)) then
        call log_message("expected L0 = 0 but got non-zero values!", level="error")
      end if
    end if
  end subroutine set_energy_balance


  !> Checks the induction equation for the equilibrium state. The two (nonzero)
  !! resulting expressions are
  !! $$ (B_{01} v_{02} - (B_{02} v_{01})' = 0, $$
  !! $$
  !! \frac{1}{\varepsilon} \bigl(\varepsilon (B_{01}v_{03} - B_{03}v_{01}) \bigr)' = 0
  !! $$
  !! and should both be fulfilled.
  !! @warning   Throws a warning if the equilibrium induction equation is not satisfied.
  subroutine induction_equil_conditions(B_field, v_field)
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)       :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)     :: v_field

    real(dp)  :: B01, B02, dB02, B03, dB03, v01, v02, v03, dv01, dv02, dv03
    real(dp)  :: eps, d_eps, r(2), discrepancy(2)
    real(dp)  :: eq_cond(gauss_gridpts, 2)
    integer   :: i, j, counter(2)
    logical   :: satisfied(2)

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    B01 = B_field % B01
    do i = 1, gauss_gridpts
      B02 = B_field % B02(i)
      B03 = B_field % B03(i)
      dB02 = B_field % d_B02_dr(i)
      dB03 = B_field % d_B03_dr(i)
      v01 = v_field % v01(i)
      v02 = v_field % v02(i)
      v03 = v_field % v03(i)
      dv01 = v_field % d_v01_dr(i)
      dv02 = v_field % d_v02_dr(i)
      dv03 = v_field % d_v03_dr(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i, 1) = B01 * dv02 - B02 * dv01 - dB02 * v01
      eq_cond(i, 2) = B01 * dv03 - dB03 * v01 - B03 * dv01 &
                      + d_eps * (B01 * v03 - B03 * v01) / eps
      do j = 1, 2
        if (abs(eq_cond(i, j)) > dp_LIMIT) then
          counter(j) = counter(j) + 1
          satisfied(j) = .false.
          if (abs(eq_cond(i, j)) > discrepancy(j)) then
            discrepancy(j) = abs(eq_cond(i, j))
            r(j) = grid_gauss(i)
          end if
        end if
      end do
    end do

    do j = 1, 2
      if (satisfied(j)) then
        cycle
      end if
      call log_message( &
        "induction equilibrium conditions not satisfied!", &
        level="warning" &
      )
      write(char_log2, int_fmt) j
      write(char_log, dp_fmt) r(j)
      call log_message( &
        "location of largest discrepancy (" // trim(adjustl(char_log2)) // &
                                      "): x = " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, exp_fmt) discrepancy(j)
      call log_message( &
        "value of largest discrepancy (" // trim(adjustl(char_log2)) // &
                                      "): " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, int_fmt) counter(j)
      call log_message( &
        "amount of nodes not satisfying criterion (" // trim(adjustl(char_log2)) // &
                                          "): " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
    end do
  end subroutine induction_equil_conditions



  !> Checks the continuity equation for the equilibrium state. This is given by
  !! $$ \frac{1}{\varepsilon} \bigl( \varepsilon \rho_0 v_{01} \bigr)' = 0. $$
  !! @warning   Throws a warning if equilibrium continuity is not satisfied.
  subroutine continuity_equil_conditions(rho_field, v_field)
    use mod_global_variables, only: gauss_gridpts, dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    !> the type containing the density attributes
    type(density_type), intent(in)      :: rho_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)  :: v_field

    real(dp)  :: rho, drho, v01, dv01
    real(dp)  :: eps, d_eps, r, discrepancy
    real(dp)  :: eq_cond(gauss_gridpts)
    integer   :: i, counter
    logical   :: satisfied

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, gauss_gridpts-1
      rho = rho_field % rho0(i)
      drho = rho_field % d_rho0_dr(i)
      v01 = v_field % v01(i)
      dv01 = v_field % d_v01_dr(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)

      eq_cond(i) = drho * v01 + rho * dv01 + rho * v01 * d_eps / eps

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
      call log_message( &
        "continuity equilibrium conditions not satisfied!", &
        level='warning' &
      )
      write(char_log, dp_fmt) r
      call log_message( &
        "location of largest discrepancy: x = " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, exp_fmt) discrepancy
      call log_message( &
        "value of largest discrepancy: " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(char_log, int_fmt) counter
      call log_message( &
        "amount of nodes not satisfying criterion: " // adjustl(trim(char_log)), &
        level='warning', &
        use_prefix=.false. &
      )
      write(*,*) ""
    end if
  end subroutine continuity_equil_conditions

end module mod_inspections
