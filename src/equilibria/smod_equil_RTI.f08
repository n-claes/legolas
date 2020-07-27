! =============================================================================
!> This submodule defines Rayleigh-Taylor instabilities in Cartesian geometry.
!! This equilibrium is a specific case of the flow driven instabilities.
!!
!!
!! This equilibrium is taken from section 13.2, p. 487 (case a) in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : used in the density profile.
!! - <tt>cte_p0</tt> = 1000 : used in the pressure profile, must be large so T0 > 0.
!! - <tt>delta</tt> = -5 : used in the density profile.
!! - <tt>g</tt> = 15 : gravitational constant.
!! - <tt>alpha</tt> = 0 : magnetic shear.
!! - <tt>theta</tt> = 0.35\(\pi\) : angle used in the velocity profile.
!! - <tt>p1</tt> = 0.2 : \(v_0\), used in the velocity profile.
!! - <tt>p2</tt> = 0.6 : \(v_1\), used in the velocity profile.
!! - <tt>p3</tt> = 0 : \(v_2\), used in the velocity profile.
!! - <tt>p4</tt> = -0.35\(\pi\) : \(\phi_0\), used in the magnetic field profile.
!! - <tt>tau</tt> = 0 : used in the velocity profile.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_RTI
  implicit none

contains

  module subroutine RTI_eq()
    use mod_equilibrium_params, only: g, delta, theta, p1, p2, p3, tau, &
                                      p4, alpha, cte_rho0, cte_p0

    if (use_defaults) then
      external_gravity = .true.

      k2 = 0.0d0
      k3 = 1.0d0
      cte_rho0 = 1.0d0
      cte_p0 = 1000.0d0
      delta = -5.0d0
      g = 15.0d0
      alpha = 0.0d0
      theta = 0.35d0 * dpi
      p1 = 0.2d0
      p2 = 0.6d0
      p3 = 0.0d0
      p4 = -0.35d0 * dpi
      tau = 0.0d0
    end if

    call flow_driven_instabilities_eq()
  end subroutine RTI_eq

end submodule smod_equil_RTI
