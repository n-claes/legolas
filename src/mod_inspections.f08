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
  use mod_logging, only: log_message, str
  use mod_settings, only: settings_t
  implicit none

  private

  public :: perform_NaN_and_negative_checks
  public :: perform_sanity_checks
  public :: check_wavenumbers

contains

  !> General routine to do initial sanity checks on the various equilibrium attributes.
  !! We check the equilibrium arrays for NaN and see if all density and temperature
  !! values are positive.
  subroutine perform_NaN_and_negative_checks( &
    rho_field, T_field, B_field, v_field, grav_field &
  )
    use mod_check_values, only: is_NaN, is_negative

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

    character(50) :: name

    name = ""

    ! TODO: is there an easier way to do this?
    if (any(is_negative(rho_field % rho0))) then
      call log_message("negative density encountered!", level="error")
      return
    end if
    if (any(is_negative(T_field % T0))) then
      call log_message("negative temperature encountered!", level="error")
      return
    end if
    if (any(is_NaN(rho_field % rho0))) then
      name = "density"
    else if (any(is_NaN(T_field % T0))) then
      name = "temperature"
    else if (is_NaN(B_field % B01)) then
      name = "B01"
    else if (any(is_NaN(B_field % B02))) then
      name = "B02"
    else if (any(is_NaN(B_field % B03))) then
      name = "B03"
    else if (any(is_NaN(v_field % v01))) then
      name = "v01"
    else if (any(is_NaN(v_field % v02))) then
      name = "v02"
    else if (any(is_NaN(v_field % v03))) then
      name = "v03"
    else if (any(is_NaN(grav_field % grav))) then
      name = "gravity"
    end if

    if (name /= "") then
      call log_message("NaN encountered in " // adjustl(trim(name)), level="error")
      return
    end if
  end subroutine perform_NaN_and_negative_checks


  !> General routine to do sanity checks on the different equilibrium types.
  !! We check the wavenumbers and on-axis values, as well as standard
  !! and non-adiabatic equilibrium force balance.
  subroutine perform_sanity_checks( &
    settings, rho_field, T_field, B_field, v_field, grav_field, rc_field, kappa_field &
  )
    !> the settings object
    type(settings_t), intent(in) :: settings
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

    call check_wavenumbers(geometry=settings%grid%get_geometry())
    call check_on_axis_values(settings, B_field, v_field)
    call standard_equil_conditions( &
      settings, rho_field, T_field, B_field, v_field, grav_field &
    )
    call continuity_equil_conditions(settings, rho_field, v_field)
    call induction_equil_conditions(settings, B_field, v_field)
    ! set the energy balance based on the equilibrium conditions
    call set_energy_balance( &
      settings, rho_field, T_field, B_field, v_field, rc_field, kappa_field &
    )
  end subroutine perform_sanity_checks


  !> Sanity check on the wavenumbers.
  !! Checks if k2 is an integer in cylindrical geometry.
  !! @warning An error if thrown if the geometry is cylindrical and k2 is not
  !!          an integer.
  subroutine check_wavenumbers(geometry)
    use mod_global_variables, only: dp_LIMIT
    use mod_equilibrium_params, only: k2

    character(len=*), intent(in) :: geometry

    if (geometry == "cylindrical") then
      ! in cylindrical geometry k2 should be an integer
      if (abs(int(k2) - k2) > dp_LIMIT) then
        call log_message( &
          "cylindrical geometry but k2 is not an integer! Value: " // str(k2), &
          level="error" &
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
  subroutine check_on_axis_values(settings, B_field, v_field)
    type(settings_t), intent(in) :: settings
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)   :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in) :: v_field

    real(dp)  :: on_axis_limit

    if (settings%grid%get_geometry() == "Cartesian") then
      return
    end if

    on_axis_limit = 1.0d-3
    if (settings%grid%get_grid_start() > on_axis_limit) then
      return
    end if

    ! LCOV_EXCL_START
    if (abs(B_field % B02(1)) > on_axis_limit) then
      call log_message( &
        "B_theta non-zero on axis! Value: " // str(B_field % B02(1)), &
        level="warning" &
      )
    end if
    if (abs(B_field % d_B03_dr(1)) > on_axis_limit) then
      call log_message( &
        "dBz/dr non-zero on axis! Value: " // str(B_field % d_B03_dr(1)), &
        level="warning" &
      )
    end if
    if (abs(v_field % v02(1)) > on_axis_limit) then
      call log_message( &
        "v_theta non-zero on axis! Value: " // str(v_field % v02(1)), &
        level="warning" &
      )
    end if
    if (abs(v_field % d_v03_dr(1)) > on_axis_limit) then
      call log_message( &
        "dvz_dr non-zero on axis! Value: " // str(v_field % d_v03_dr(1)), &
        level="warning" &
      )
    end if
    ! LCOV_EXCL_STOP
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
  subroutine standard_equil_conditions( &
    settings, rho_field, T_field, B_field, v_field, grav_field &
  )
    use mod_global_variables, only: dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    type(settings_t), intent(in) :: settings
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
    real(dp)  :: eq_cond(settings%grid%get_gauss_gridpts(), 3)
    integer   :: i, j, counter(3)
    logical   :: satisfied(3)

    B01 = B_field % B01
    if ( &
      (settings%grid%get_geometry() == "cylindrical") &
      .and. (abs(B01) > dp_LIMIT) &
    ) then
      call log_message( &
        "B01 component currently not supported for cylindrical geometries!", &
        level="error" &
      )
    end if

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts()
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
      ! LCOV_EXCL_START
      call log_message( &
        "standard equilibrium conditions not satisfied!", level="warning" &
      )
      call log_message( &
        "location of largest discrepancy (" // str(j) // "): x = " // str(r(j)), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "value of largest discrepancy (" // str(j) // "): " &
        // str(discrepancy(j), fmt="e20.8"), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "amount of nodes not satisfying criterion (" // str(j) // "): " &
        // str(counter(j)), &
        level="warning", &
        use_prefix=.false. &
      )
      write(*,*) ""
      ! LCOV_EXCL_STOP
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
    settings, rho_field, T_field, B_field, v_field, rc_field, kappa_field &
  )
    use mod_global_variables, only: dp_LIMIT
    use mod_grid, only: eps_grid, d_eps_grid_dr

    !> the settings object
    type(settings_t), intent(in) :: settings
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
    do i = 1, settings%grid%get_gauss_gridpts()
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
        + dT0 * v01 / settings%physics%get_gamma_1() &
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
  subroutine induction_equil_conditions(settings, B_field, v_field)
    use mod_global_variables, only: dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    type(settings_t), intent(in) :: settings
    !> the type containing the magnetic field attributes
    type(bfield_type), intent(in)       :: B_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)     :: v_field

    real(dp)  :: B01, B02, dB02, B03, dB03, v01, v02, v03, dv01, dv02, dv03
    real(dp)  :: eps, d_eps, r(2), discrepancy(2)
    real(dp)  :: eq_cond(settings%grid%get_gauss_gridpts(), 2)
    integer   :: i, j, counter(2)
    logical   :: satisfied(2)

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    B01 = B_field % B01
    do i = 1, settings%grid%get_gauss_gridpts()
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
      ! LCOV_EXCL_START
      call log_message( &
        "induction equilibrium conditions not satisfied!", level="warning" &
      )
      call log_message( &
        "location of largest discrepancy (" // str(j) // "): x = " // str(r(j)), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "value of largest discrepancy (" // str(j) // "): " &
        // str(discrepancy(j), fmt="e20.8"), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "amount of nodes not satisfying criterion (" // str(j) // "): " &
        // str(counter(j)), &
        level="warning", &
        use_prefix=.false. &
      )
    end do
    ! LCOV_EXCL_STOP
  end subroutine induction_equil_conditions



  !> Checks the continuity equation for the equilibrium state. This is given by
  !! $$ \frac{1}{\varepsilon} \bigl( \varepsilon \rho_0 v_{01} \bigr)' = 0. $$
  !! @warning   Throws a warning if equilibrium continuity is not satisfied.
  subroutine continuity_equil_conditions(settings, rho_field, v_field)
    use mod_global_variables, only: dp_LIMIT
    use mod_grid, only: grid_gauss, eps_grid, d_eps_grid_dr

    type(settings_t), intent(in) :: settings
    !> the type containing the density attributes
    type(density_type), intent(in)      :: rho_field
    !> the type containing the velocity attributes
    type(velocity_type), intent(in)  :: v_field

    real(dp)  :: rho, drho, v01, dv01
    real(dp)  :: eps, d_eps, r, discrepancy
    real(dp)  :: eq_cond(settings%grid%get_gauss_gridpts())
    integer   :: i, counter
    logical   :: satisfied

    satisfied = .true.
    discrepancy = 0.0d0
    counter = 0
    do i = 1, settings%grid%get_gauss_gridpts()-1
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

    ! LCOV_EXCL_START
    if (.not. satisfied) then
      call log_message( &
        "continuity equilibrium conditions not satisfied!", &
        level="warning" &
      )
      call log_message( &
        "location of largest discrepancy: x = " // str(r), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "value of largest discrepancy: " // str(discrepancy, fmt="e20.8"), &
        level="warning", &
        use_prefix=.false. &
      )
      call log_message( &
        "amount of nodes not satisfying criterion: " // str(counter), &
        level="warning", &
        use_prefix=.false. &
      )
      write(*,*) ""
    end if
    ! LCOV_EXCL_STOP
  end subroutine continuity_equil_conditions

end module mod_inspections
