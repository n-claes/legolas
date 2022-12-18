! =============================================================================
!> This submodule defines Kelvin-Helmholtz instabilities in Cartesian geometry.
!! This equilibrium is a specific case of the flow driven instabilities.
!! No magnetic fields are considered in this case (pure HD).
!!
!!
!! This equilibrium is taken from section 13.2, p. 487 (case b) in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : used in the density profile.
!! - <tt>cte_p0</tt> = 10 : used in the pressure profile.
!! - <tt>delta</tt> = 0 : used in the density profile.
!! - <tt>g</tt> = 0 : gravitational constant.
!! - <tt>alpha</tt> = 0 : magnetic shear.
!! - <tt>theta</tt> = 0 : angle used in the velocity profile.
!! - <tt>p1</tt> = 0 : \(v_0\), used in the velocity profile.
!! - <tt>p2</tt> = 0 : \(v_1\), used in the velocity profile.
!! - <tt>p3</tt> = 1 : \(v_2\), used in the velocity profile.
!! - <tt>p4</tt> = 0 : \(\phi_0\), used in the magnetic field profile.
!! - <tt>tau</tt> = 11 : used in the velocity profile.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_KHI
  implicit none

contains

  module procedure KHI_eq
    use mod_equilibrium_params, only: g, delta, theta, p1, p2, p3, tau, &
      p4, alpha, cte_rho0, cte_p0

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      k2 = 0.0d0
      k3 = 1.0d0
      cte_rho0 = 1.0d0
      cte_p0 = 1000d0
      delta = 0.0d0
      g = 0.0d0
      alpha = 0.0d0
      theta = 0.0d0
      p1 = 0.0d0
      p2 = 0.0d0
      p3 = 1.0d0
      p4 = 0.0d0
      tau = 11.0d0
    end if ! LCOV_EXCL_STOP

    call flow_driven_instabilities_eq(settings)
    ! manually force the magnetic field to zero, it's HD here (the above set
    ! of parameters still yield a B03 component)
    B_field % B02 = 0.0d0
    B_field % B03 = 0.0d0
    B_field % B0 = 0.0d0
    B_field % d_B02_dr = 0.0d0
    B_field % d_B03_dr = 0.0d0
  end procedure KHI_eq

end submodule smod_equil_KHI
