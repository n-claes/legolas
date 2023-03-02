! =============================================================================
!> This submodule defines an equilibrium containing magnetothermal instabilities
!! in a cylindrical geometry. The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Van der Linden, R. A. M., Goossens, M., & Hood, A. W. (1992).
!!  The relevance of the ballooning approximation for magnetic, thermal,
!!  and coalesced magnetothermal instabilities. Solar physics, 140(2), 317-342_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_T0</tt> = 1 : used to set the temperature.
!! - cooling_curve = 'rosner'
!! - parallel thermal conduction, no perpendicular conduction
!!
!! and normalisations given by
!!
!! - <tt>unit_temperature</tt> = 2.6e6 K
!! - <tt>unit_magneticfield</tt> = 10 Gauss
!! - <tt>unit_length</tt> = 1e8 cm
!!
!! and can all be changed in the parfile. @endnote
!! @note The default setup handles _CASE I_ in the original paper. For _CASE II_,
!!       you can set k2 = 10. To reproduce the thermal continuum plot,
!!       set <tt>unit_length</tt> = 2.44e9 cm in _CASE I_. @endnote
submodule(mod_equilibrium) smod_equil_magnetothermal_instabilities
  implicit none

contains

  !> Sets the equilibrium
  module procedure magnetothermal_instability_eq
    use mod_equilibrium_params, only: cte_T0

    real(dp)  :: r
    real(dp), allocatable :: p_r(:)
    integer   :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_cooling(cooling_curve="rosner")
      call settings%physics%enable_parallel_conduction()
      call settings%units%set_units_from_temperature( &
        unit_temperature=2.6d6, &
        unit_magneticfield=10.0d0, &
        unit_length=1.00d8, &
        mean_molecular_weight=1.0d0 & ! this work assumes pure proton plasma
      )

      cte_T0 = 1.0d0
      k2 = 0.0d0
      k3 = 1.0d0
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    allocate(p_r(settings%grid%get_gauss_gridpts()))
    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      B_field % B02(i) = r / (1.0d0 + r**2)
      B_field % B03(i) = 0.0d0
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      p_r(i) = 1.0d0 / ( 2.0d0 * (1.0d0 + r**2)**2 )
      T_field % T0(i) = cte_T0
      rho_field % rho0(i) = p_r(i) / cte_T0

      rho_field % d_rho0_dr(i) = -2.0d0 * r / (cte_T0 * (r**2 + 1.0d0)**3)
      B_field % d_B02_dr(i) = (1.0d0 - r**2) / (r**4 + 2.0d0 * r**2 + 1.0d0)
      ! Temperature derivative is zero due to isothermal
    end do
    deallocate(p_r)
  end procedure magnetothermal_instability_eq

end submodule smod_equil_magnetothermal_instabilities
