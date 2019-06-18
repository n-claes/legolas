module mod_setup_matrix_a
  use mod_global_variables
  implicit none

  private

  ! Sets up the A-matrix for the eigenvalue problem wBX = AX
  real(dp)                 :: h_cubic(4), h_quadratic(4)
  real(dp)                 :: dh_cubic_dr(4), dh_quadratic_dr(4)
  ! Factors and positions are allocatable so they are dynamic, in case we add
  ! additional equations (self-gravity)
  complex(dp), allocatable :: factors_A(:)
  integer, allocatable     :: positions(:, :)

  public  :: construct_A
  public  :: matrix_A_clean

contains

  subroutine construct_A(matrix_A)
    use mod_grid
    use mod_equilibrium
    use mod_spline_functions

    complex(dp), intent(out)  :: matrix_A(matrix_gridpts, matrix_gridpts)
    complex(dp)               :: quadblock(dim_quadblock, dim_quadblock)
    real(dp)                  :: r_lo, r_hi, eps, d_eps_dr, curr_weight
    real(dp)                  :: r
    integer                   :: i, j, gauss_idx, k, l

    ! Initialise matrix to zero (A is complex)
    matrix_A = (0.0d0, 0.0d0)

    do i = 2, gridpts-1
      ! Set quadblock to zero
      quadblock = (0.0d0, 0.0d0)

      r_lo = grid(i)
      r_hi = grid(i + 1)

      do j = 1, n_gauss
        ! Current grid index (from 4*n_gauss points)
        gauss_idx = (i - 1)*n_gauss + j
        r = grid_gauss(gauss_idx)
        if ((geometry .eq. "cartesian") .or. (geometry .eq. "Cartesian")) then
          eps      = 1.0d0
          d_eps_dr = 0.0d0
        else if (geometry .eq. "cylindrical") then
          eps      = r
          d_eps_dr = 1.0d0
        else
          write(*,*) "Geometry not defined correctly."
          write(*,*) "Currently set on:    ", geometry
          stop
        end if

        curr_weight = gaussian_weights(j)

        call quadratic_factors(r, r_lo, r_hi, h_quadratic)
        call quadratic_factors_deriv(r, r_lo, r_hi, dh_quadratic_dr)
        call cubic_factors(r, r_lo, r_hi, h_cubic)
        call cubic_factors_deriv(r, r_lo, r_hi, dh_cubic_dr)

        call get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock)

      end do  ! end do iteration Gaussian points

      ! Gridpoint i = 1, place quadblock in upper left corner of A
      if (i == 1) then
        do k = 1, dim_quadblock
          do l = 1, dim_quadblock
            matrix_A(k, l) = quadblock(k, l)
          end do
        end do
      end if

      ! TODO: fill matrix A with quadblocks

    end do    ! end do iteration gridpoints

    deallocate(factors_A)
    deallocate(positions)

  end subroutine construct_A

  subroutine get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock)
    use mod_grid
    use mod_equilibrium
    use mod_equilibrium_derivatives
    use mod_make_subblock
    use mod_gravity, only: grav

    integer, intent(in)       :: gauss_idx
    real(dp), intent(in)      :: eps, d_eps_dr, curr_weight
    complex(dp), intent(out)  :: quadblock(dim_quadblock, dim_quadblock)

    real(dp)                  :: eps_inv
    real(dp)                  :: rho0, T0, B01, B02, B03, B2_inv
    real(dp)                  :: v01, v02, v03
    real(dp)                  :: tc_para, tc_perp, L0, eta

    real(dp)                  :: drho0, drB02, dB02_r, dB03, dT0, dB02
    real(dp)                  :: drv02, dv03
    real(dp)                  :: dtc_perp_dT, dtc_perp_drho, dtc_perp_dB2
    real(dp)                  :: L_T, L_rho
    real(dp)                  :: deta, ddB03, ddB02


    complex(dp)               :: ic

    ic      = (0.0d0, 1.0d0)
    eps_inv = 1.0d0 / eps


    !! Equilibrium quantities
    !! Default variables
    rho0    = rho0_eq(gauss_idx)
    T0      = T0_eq(gauss_idx)
    B01     = B01_eq(gauss_idx)
    B02     = B02_eq(gauss_idx)
    B03     = B03_eq(gauss_idx)
    B2_inv  = 1.0d0 / (B0_eq(gauss_idx)**2)
    !! Flow
    v01     = v01_eq(gauss_idx)
    v02     = v02_eq(gauss_idx)
    v03     = v03_eq(gauss_idx)
    !! Thermal conduction
    tc_para = tc_para_eq(gauss_idx)
    tc_perp = tc_perp_eq(gauss_idx)
    !! Radiative cooling
    L0      = heat_loss_eq(gauss_idx)
    !! Resistivity
    eta     = eta_eq(gauss_idx)

    !! Derivatives of equilibrium quantities
    !! Default derivatives
    drho0   = d_rho0_dr(gauss_idx)
    drB02   = d_rB02_dr(gauss_idx)
    dB02_r  = d_B02_r_dr(gauss_idx)
    dB03    = d_B03_dr(gauss_idx)
    dT0     = d_T0_dr(gauss_idx)
    dB02    = d_B02_dr(gauss_idx)
    !! Flow
    drv02   = d_rv02_dr(gauss_idx)
    dv03    = d_v03_dr(gauss_idx)
    !! Thermal conduction
    dtc_perp_dT   = d_tc_perp_eq_dT(gauss_idx)
    dtc_perp_drho = d_tc_perp_eq_drho(gauss_idx)
    dtc_perp_dB2  = d_tc_perp_eq_dB2(gauss_idx)
    !! Radiative cooling
    L_T     = d_L_dT(gauss_idx)
    L_rho   = d_L_drho(gauss_idx)
    !! Resistivity
    deta    = d_eta_dT(gauss_idx)
    ddB02   = dd_B02_dr(gauss_idx)
    ddB03   = dd_B03_dr(gauss_idx)


    !! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors_A(19)
    call reset_positions(19)

    ! A(1, 1)
    factors_A(1) = eps_inv * (eps_inv * v02 * k2 + v03 * k3) * d_eps_dr
    positions(1, :) = [1, 1]
    ! A(1, 3)
    factors_A(2) = eps_inv * rho0 * k2
    positions(2, :) = [1, 3]
    ! A(1, 4)
    factors_A(3) = eps_inv * rho0 * k3
    positions(3, :) = [1, 4]
    ! A(3, 1)
    factors_A(4) = eps_inv**2 * T0 * k2
    positions(4, :) = [3, 1]
    ! A(3, 3)
    factors_A(5) = rho0 * (eps_inv * k2 * v02 + k3 * v03)
    positions(5, :) = [3, 3]
    ! A(3, 5)
    factors_A(6) = eps_inv**2 * rho0 * k2
    positions(6, :) = [3, 5]
    ! A(3, 6)
    factors_A(7) = -B03 * (k3**2 + k2**2 * eps_inv**2)
    positions(7, :) = [3, 6]
    ! A(4, 1)
    factors_A(8) = eps_inv * T0 * k3
    positions(8, :) = [4, 1]
    ! A(4, 4)
    factors_A(9) = eps_inv * rho0 * (eps_inv * k2 * v02 + k3 * v03)
    positions(9, :) = [4, 4]
    ! A(4, 5)
    factors_A(10) = eps_inv * rho0 * k3
    positions(10, :) = [4, 5]
    ! A(4, 6)
    factors_A(11) = B02 * (k3**2 + k2**2 * eps_inv**2)
    positions(11, :) = [4, 6]
    ! A(5, 1)
    factors_A(12) = ic * gamma_1 * eps_inv * &
                    (d_eps_dr * eps_inv * dT0 * dtc_perp_dT - L0 - rho0*L_rho)
    positions(12, :) = [5, 1]
    ! A(5, 3)
    factors_A(13) = gamma_1 * eps_inv * rho0 * T0 * k2
    positions(13, :) = [5, 3]
    ! A(5, 4)
    factors_A(14) = gamma_1 * eps_inv * rho0 * T0 * k3
    positions(14, :) = [5, 4]
    ! A(5, 5)
    factors_A(15) = -eps_inv * ic * gamma_1 * ( &
                    (tc_para - tc_perp) * B2_inv * (k2 * eps_inv * B02 + k3*B03) &
                     + tc_perp * (d_eps_dr * eps_inv)**2 &
                     + tc_perp * (k2**2 * eps_inv**2 + k3**2) &
                     + rho0 * L_T - eps_inv * dT0 * dtc_perp_dT &
                                              )
    positions(15, :) = [5, 5]
    ! A(5, 6)
    factors_A(16) = ic * gamma_1 * d_eps_dr * eps_inv * 2 * dT0 * &
                    (eps * B02 * k3 - B03 * k2) * dtc_perp_dB2
    positions(16, :) = [5, 6]
    ! A(6, 3)
    factors_A(17) = -B03
    positions(17, :) = [6, 3]
    ! A(6, 4)
    factors_A(18) = -eps_inv * B02
    positions(18, :) = [6, 4]
    ! A(6, 6)
    factors_A(19) = k2 * eps_inv * v02 + k2 * v03
    positions(19, :) = [6, 6]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_quadratic, h_quadratic)


    ! Quadratic * Cubic
    call reset_factors_A(10)
    call reset_positions(10)

    ! A(1, 2)
    factors_A(1) = -eps_inv * drho0
    positions(1, :) = [1, 2]
    ! A(3, 2)
    factors_A(2) = -eps_inv**2 * rho0 * drv02
    positions(2, :) = [3, 2]
    ! A(3, 7)
    factors_A(3) = eps_inv**2 * drB02 * k3
    positions(3, :) = [3, 7]
    ! A(3, 8)
    factors_A(4) = -eps_inv**2 * drB02 * k2
    positions(4, :) = [3, 8]
    ! A(4, 2)
    factors_A(5) = -eps_inv * rho0 * dv03
    positions(5, :) = [4, 2]
    ! A(4, 7)
    factors_A(6) = eps_inv * dB03 * k3
    positions(6, :) = [4, 7]
    ! A(4, 8)
    factors_A(7) = -eps_inv * dB03 * k2
    positions(7, :) = [4, 8]
    ! A(5, 2)
    factors_A(8) = -eps_inv * rho0 * dT0
    positions(8, :) = [5, 2]
    ! A(5, 7)
    factors_A(9) = ic * gamma_1 * eps_inv * (tc_para - tc_perp) * B2_inv &
                   * dT0 * (k2 * k3 * B02 + k3**2 * B02)
    positions(9, :) = [5, 7]
    ! A(5, 8)
    factors_A(10) = ic * gamma_1 * eps_inv * (tc_para - tc_perp) * B2_inv &
                    * dT0 * (k2**2 * B02 + k2 * k3 * B03)
    positions(10, :) = [5, 8]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_quadratic, h_cubic)


    ! Quadratic * d(Cubic)/dr
    call reset_factors_A(10)
    call reset_positions(10)

    ! A(1, 2)
    factors_A(1) = -eps_inv * rho0
    positions(1, :) = [1, 2]
    ! A(3, 7)
    factors_A(2) = eps_inv**2 * B03 * k2
    positions(2, :) = [3, 7]
    ! A(3, 8)
    factors_A(3) = B03 * k3
    positions(3, :) = [3, 8]
    ! A(4, 7)
    factors_A(4) = -eps_inv**2 * B02 * k2
    positions(4, :) = [4, 7]
    ! A(4, 8)
    factors_A(5) = -B02 * k3
    positions(5, :) = [4, 8]
    ! A(5, 2)
    factors_A(6) = -gamma_1 * eps_inv * rho0 * T0
    positions(6, :) = [5, 2]
    ! A(5, 7)
    factors_A(7) = 2*ic*gamma_1*eps_inv * &
                  (dT0 * B03 * d_eps_dr * eps_inv * dtc_perp_dB2 - eta * ddB03)
    positions(7, :) = [5, 7]
    ! A(5, 8)    (derivative of eta-term has been rewritten)
    factors_A(8) = -2*ic*gamma_1 * &
                    (dT0 * d_eps_dr * eps_inv * B02 * dtc_perp_dB2 &
                     - eta*ddB02 + 2*(d_eps_dr * eps_inv)**2 * eta * B02)
    positions(8, :) = [5, 8]
    ! A(6, 7)
    factors_A(9) = -eps_inv * v02 + ic * eta * eps_inv**2 * k2
    positions(9, :) = [6, 7]
    ! A(6, 8)
    factors_A(10) = -v03 + ic * eta * k3
    positions(10, :) = [6, 8]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_quadratic, dh_cubic_dr)


    ! Cubic * Quadratic
    call reset_factors_A(6)
    call reset_positions(6)

    ! A(2, 1)
    factors_A(1) = -d_eps_dr * eps_inv**2 * v02**2 + eps_inv * grav
    positions(1, :) = [2, 1]
    ! A(2, 3)
    factors_A(2) = -2 * d_eps_dr * eps_inv * rho0 * v02
    positions(2, :) = [2, 3]
    ! A(2, 6)
    factors_A(3) = eps_inv * drB02 * k3 - eps * k3 * db02_r
    positions(3, :) = [2, 6]
    ! A(7, 5)
    factors_A(4) = -ic * eps_inv * deta * dB03
    positions(4, :) = [7, 5]
    ! A(8, 5)
    factors_A(5) = -ic * eps_inv**2 * deta * drB02
    positions(5, :) = [8, 5]
    ! A(8, 6)
    factors_A(6) = -ic * eta * d_eps_dr * eps_inv * k3
    positions(6, :) = [8, 6]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_cubic, h_quadratic)


    ! d(Cubic)/dr * Quadratic
    call reset_factors_A(5)
    call reset_positions(5)

    ! A(2, 1)
    factors_A(1) = -T0 * eps_inv
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors_A(2) = -rho0 * eps_inv
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors_A(3) = eps_inv * k2 * B03 - k3 * B02
    positions(3, :) = [2, 6]
    ! A(7, 6)
    factors_A(4) = -ic * eta * eps_inv * k2
    positions(4, :) = [7, 6]
    ! A(8, 6)
    factors_A(5) = ic * eta * k3
    positions(5, :) = [8, 6]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  dh_cubic_dr, h_quadratic)


    ! Cubic * Cubic
    call reset_factors_A(9)
    call reset_positions(9)

    ! A(2, 2)
    factors_A(1) = eps_inv * (eps_inv * v02 * k2 + v03 * k3) * rho0
    positions(1, :) = [2, 2]
    ! A(2, 7)
    factors_A(2) = -B03 * eps_inv * k3**2 - B02 * eps_inv**2 * k2 * k3
    positions(2, :) = [2, 7]
    ! A(2, 8)
    factors_A(3) = B03 * eps_inv * k2 * k3 + B02 * eps_inv**2 * k2**2
    positions(3, :) = [2, 8]
    ! A(7, 2)
    factors_A(4) = -eps_inv * B03
    positions(4, :) = [7, 2]
    ! A(7, 7)
    factors_A(5) = eps_inv * k3 * (v03 - ic * eta * k3)
    positions(5, :) = [7, 7]
    ! A(7, 8)
    factors_A(6) = -eps_inv * k2 * (v03 - ic * eta * k3)
    positions(6, :) = [7, 8]
    ! A(8, 2)
    factors_A(7) = eps_inv * B02
    positions(7, :) = [8, 2]
    ! A(8, 7)
    factors_A(8) = -eps_inv * k3 * (v02 - ic * eta * eps_inv * k2)
    positions(8, :) = [8, 7]
    ! A(8, 8)
    factors_A(9) = eps_inv * k2 * (v02 - ic * eta * eps_inv * k2)
    positions(9, :) = [8, 8]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_cubic, h_cubic)


    ! Cubic * d(Cubic)/dr
    call reset_factors_A(2)
    call reset_positions(2)

    ! A(2, 8)
    factors_A(1) = -2 * B02 * d_eps_dr * eps_inv
    positions(1, :) = [2, 8]
    ! A(8, 8)
    factors_A(2) = ic * eta * d_eps_dr * eps_inv
    positions(2, :) = [8, 8]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_cubic, dh_cubic_dr)


    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors_A(4)
    call reset_positions(4)

    ! A(2, 7)
    factors_A(1) = -B03 * eps_inv
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors_A(2) = B02
    positions(2, :) = [2, 8]
    ! A(7, 7)
    factors_A(3) = -ic * eta * eps_inv
    positions(3, :) = [7, 7]
    ! A(8, 8)
    factors_A(4) = -ic * eta
    positions(4, :) = [8, 8]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  dh_cubic_dr, dh_cubic_dr)


    ! d(Quadratic)/dr * Quadratic
    call reset_factors_A(3)
    call reset_positions(3)

    ! A(5, 1)
    factors_A(1) = -ic * gamma_1 * eps_inv * dT0 * dtc_perp_drho
    positions(1, :) = [5, 1]
    ! A(5, 5)
    factors_A(2) = ic * gamma_1 * eps_inv * (d_eps_dr * eps_inv * tc_perp &
                                             - dT0 * dtc_perp_dT)
    positions(2, :) = [5, 5]
    ! A(5, 6)
    factors_A(3) = -2 * ic * gamma_1 * eps_inv * &
                   (dT0 * (eps * B02 * k3 - B03 * k2) * dtc_perp_dB2  &
                    - eta * k2 * dB03 + eta * k3 * drB02 )
    positions(3, :) = [5, 6]

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  dh_quadratic_dr, h_quadratic)



    return

  end subroutine get_A_elements



  subroutine reset_factors_A(size_factors_A)
    integer, intent(in)       :: size_factors_A

    if (allocated(factors_A)) then
      deallocate(factors_A)
    end if

    allocate(factors_A(size_factors_A))
    factors_A = (0.0d0, 0.0d0)

  end subroutine reset_factors_A



  subroutine reset_positions(size_positions)
    integer, intent(in)       :: size_positions

    if (allocated(positions)) then
      deallocate(positions)
    end if

    allocate(positions(size_positions, 2))
    positions = 0

  end subroutine reset_positions


  subroutine matrix_A_clean()
    if (allocated(positions)) then
      deallocate(positions)
    end if
    if (allocated(factors_A)) then
      deallocate(factors_A)
    end if
  end subroutine matrix_A_clean

end module mod_setup_matrix_a
