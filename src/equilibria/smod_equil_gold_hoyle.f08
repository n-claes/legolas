! =============================================================================
!> This submodule defines a Gold-Hoyle equilibrium in cylindrical geometry.
!! This equilibrium configuration models a filament with a uniform twist
!! such that all fieldlines perform an equal amount of turns around the cylinder axis.
!! The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from
!! _Van der Linden, R. A. M., Goossens, M. (1991).
!! "The thermal continuum in coronal loops: instability criteria and the influence of
!!  perpendicular thermal conduction." Solar physics, 134.2, 247-273._
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = 1
!! - <tt>cte_T0</tt> = 1 : used to set the temperature.
!! - <tt>cte_rho0</tt> = 1 : used to set the density
!! - <tt>alpha</tt> = 20 : used in the magnetic field components
!! - cooling_curve = 'rosner'
!! - parallel thermal conduction, no perpendicular conduction
!!
!! and can all be changed in the parfile. @endnote
!! @note To reproduce the three profiles of the original paper you can supply
!!       one of the following normalisations:
!!
!! - cold plasma: <tt>unit_magneticfield</tt> = 10 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e10 cm,
!!                corresponding to a temperature of 5.7e4 K.
!! - cool plasma: <tt>unit_magneticfield</tt> = 22.5 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e10 cm,
!!                corresponding to a temperature of 2.9e5 K (default).
!! - hot plasma:  <tt>unit_magneticfield</tt> = 67 Gauss,
!!                <tt>unit_density</tt> = 1.6726e-15 g/cm3,
!!                <tt>unit_length</tt> = 1.0e11 cm,
!!                corresponding to a temperature of 2.6e6 K. @endnote
submodule (mod_equilibrium) smod_equil_gold_hoyle
  implicit none

contains

  module procedure gold_hoyle_eq
    use mod_equilibrium_params, only: cte_T0, cte_rho0, alpha

    real(dp)  :: r
    integer   :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_cooling(cooling_curve="rosner")
      call settings%physics%enable_parallel_conduction()

      k2 = 1.0d0
      k3 = 1.0d0
      cte_rho0 = 1.0d0
      cte_T0 = 0.001d0
      alpha = 20.0d0

      ! cool plasma: B = 22.5 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 2.9e5 K
      call settings%units%set_units_from_density( &
        unit_density=1.6727d-15, &
        unit_magneticfield=22.5d0, &
        unit_length=1.0d10, &
        mean_molecular_weight=1.0d0 & ! pure proton plasma
      )
      ! hot plasma: B = 67.0 G, rho = 1.6726e-12 kg/m3, R = 1e9 m, T = 2.6e6 K
      ! call settings%units%set_units_from_density( &
      !   unit_density=1.6727d-15, &
      !   unit_magneticfield=67.0d0, &
      !   unit_length=1.0d11, &
      !   mean_molecular_weight=1.0d0 & ! pure proton plasma
      ! )
      ! cold plasma: B = 10.0 G, rho = 1.6726e-12 kg/m3, R = 1e8 m, T = 5.7e4 K
      ! call settings%units%set_units_from_density( &
      !   unit_density=1.6727d-15, &
      !   unit_magneticfield=10.0d0, &
      !   unit_length=1.0d10, &
      !   mean_molecular_weight=1.0d0 & ! pure proton plasma
      ! )
    end if ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      T_field % T0(i) = cte_T0
      B_field % B02(i) = alpha * r / (1.0d0 + alpha**2 * r**2)
      B_field % B03(i) = 1.0d0 / (1.0d0 + alpha**2 * r**2)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)

      B_field % d_B02_dr(i) = -(alpha**3 * r**2 - alpha) / ( &
        alpha**4 * r**4 + 2.0d0 * alpha**2 * r**2 + 1.0d0 &
      )
      B_field % d_B03_dr(i) = -2.0d0 * alpha**2 * r / (alpha**2 * r**2 + 1.0d0)**2
    end do
  end procedure gold_hoyle_eq

end submodule smod_equil_gold_hoyle
