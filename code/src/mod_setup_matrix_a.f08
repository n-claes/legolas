module mod_setup_matrix_a
  use mod_global_variables
  implicit none

  ! Sets up the A-matrix for the eigenvalue problem wBX = AX
  real(dp)                 :: h_cubic(4), h_quadratic(4)
  real(dp)                 :: dh_cubic_dr(4), dh_quadratic_dr(4)
  ! Factors and positions are allocatable so they are dynamic, in case we add
  ! additional equations (self-gravity)
  complex(dp), allocatable :: factors_A(:)
  integer, allocatable     :: positions(:, :)

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

    integer, intent(in)       :: gauss_idx
    real(dp), intent(in)      :: eps, d_eps_dr, curr_weight
    complex(dp), intent(out)  :: quadblock(dim_quadblock, dim_quadblock)

    real(dp)                  :: eps_inv
    real(dp)                  :: rho0, v01, v02, v03, T0, B01, B02, B03
    real(dp)                  :: B2_inv
    real(dp)                  :: drho0, drB02, dB03_r, drv02, dv03, dB03, dT0

    real(dp)                  :: L0, L_T, L_rho
    real(dp)                  :: tc_para, tc_perp, d_tc_perp_dT, &
                                 d_tc_perp_drho, d_tc_perp_dB2

    complex(dp)               :: ic

    ic      = (0.0d0, 1.0d0)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
    rho0    = rho0_eq(gauss_idx)
    v01     = v01_eq(gauss_idx)
    v02     = v02_eq(gauss_idx)
    v03     = v03_eq(gauss_idx)
    T0      = T0_eq(gauss_idx)
    B01     = B01_eq(gauss_idx)
    B02     = B02_eq(gauss_idx)
    B03     = B03_eq(gauss_idx)
    B2_inv  = 1.0d0 / (B01**2 + B02**2 + B03**2)

    ! Derivatives of equilibrium quantities
    drho0   = d_rho0_dr(gauss_idx)
    drB02   = d_rB02_dr(gauss_idx)
    dB03_r  = d_B03_r_dr(gauss_idx)
    drv02   = d_rv02_dr(gauss_idx)
    dv03    = d_v03_dr(gauss_idx)
    dB03    = d_B03_dr(gauss_idx)
    dT0     = d_T0_dr(gauss_idx)

    ! Radiative cooling quantities (all 0.0d0 if set to false)
    L0      = heat_loss_eq(gauss_idx)
    L_T     = d_L_dT(gauss_idx)
    L_rho   = d_L_drho(gauss_idx)

    ! Thermal conduction quantities (all 0.0d0 if set to false)
    tc_para = tc_para_eq(gauss_idx)
    tc_perp = tc_perp_eq(gauss_idx)
    d_tc_perp_dT   = d_tc_perp_eq_dT(gauss_idx)
    d_tc_perp_drho = d_tc_perp_eq_drho(gauss_idx)
    d_tc_perp_dB2  = d_tc_perp_eq_dB2(gauss_idx)

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
                    (d_eps_dr * eps_inv * dT0 * d_tc_perp_dT - L0 - rho0*L_rho)
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
                     + rho0 * L_T - eps_inv * dT0 * d_tc_perp_dT &
                                              )
    positions(15, :) = [5, 5]
    ! A(5, 6)
    factors_A(16) = ic * gamma_1 * d_eps_dr * eps_inv * 2 * dT0 * &
                    (eps * B02 * k3 - B03 * k2) * d_tc_perp_dB2
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
    factors_A(4) = -eps_inv**2 drB02 * k2
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
