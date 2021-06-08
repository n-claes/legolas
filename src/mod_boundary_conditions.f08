! =============================================================================
!> Imposes boundary conditions on the finite element matrices after assembly.
!! Both natural and essential conditions are handled.
module mod_boundary_conditions
  use mod_global_variables, only: dp, matrix_gridpts, dim_quadblock, dim_subblock
  implicit none

  private

  !> boolean to check presence of perpendicular thermal conduction
  logical, save   :: kappa_perp_is_zero

  public :: apply_boundary_conditions

contains


  !> Main routine to handle the boundary conditions. This subroutine first performs some checks
  !! and then handles the natural and essential boundary conditions.
  !! @note  Perpendicular thermal conduction is hard-checked, that is, not just the
  !!        global logical. There must be a value in the corresponding array that is non-zero.
  subroutine apply_boundary_conditions(matrix_A, matrix_B)
    use mod_equilibrium, only: kappa_field
    ! use mod_viscosity, only: viscosity_boundaries
    use mod_global_variables, only: dp_LIMIT
    use mod_equilibrium, only: kappa_field

    !> the A-matrix with boundary conditions imposed on exit
    complex(dp), intent(inout)  :: matrix_A(matrix_gridpts, matrix_gridpts)
    !> the B-matrix with boundary conditions imposed on exit
    real(dp), intent(inout)     :: matrix_B(matrix_gridpts, matrix_gridpts)
    complex(dp)                 :: quadblock(dim_quadblock, dim_quadblock)
    ! complex(dp)                 :: quadblock_visc(dim_quadblock, dim_quadblock)
    ! complex(dp)                 :: quadblock_Hall(dim_quadblock, dim_quadblock)
    integer                     :: idx_end_left, idx_start_right
    ! real(dp)                    :: hf

    ! check if perpendicular thermal conduction is present
    kappa_perp_is_zero = .true.
    if (any(abs(kappa_field % kappa_perp) > dp_LIMIT)) then
      kappa_perp_is_zero = .false.
    end if

    ! end of index first-gridpoint quadblock
    idx_end_left = dim_quadblock
    ! start of index last-gridpoint quadblock
    idx_start_right = matrix_gridpts - dim_quadblock + 1

    ! matrix B left-edge quadblock
    quadblock = matrix_B(1:idx_end_left, 1:idx_end_left)
    call essential_boundaries(quadblock, edge='l_edge', matrix='B')
    matrix_B(1:idx_end_left, 1:idx_end_left) = real(quadblock)
    ! matrix B right-edge quadblock
    quadblock = matrix_B(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts)
    call essential_boundaries(quadblock, edge='r_edge', matrix='B')
    matrix_B(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts) = real(quadblock)

    ! matrix A left-edge quadblock
    quadblock = matrix_A(1:idx_end_left, 1:idx_end_left)
    call essential_boundaries(quadblock, edge='l_edge', matrix='A')
    call natural_boundaries(quadblock, edge='l_edge')
    ! if (viscosity) then
    !   call viscosity_boundaries(quadblock_visc, edge='l_edge')
    !   quadblock = quadblock + quadblock_visc
    ! end if
    ! if (hall_mhd) then
    !   call hall_boundaries(quadblock_Hall, kappa_perp_is_zero, edge='l_edge')
    !   hf = hall_field % hallfactor(1)
    !   quadblock = quadblock - hf * quadblock_Hall
    ! end if
    matrix_A(1:idx_end_left, 1:idx_end_left) = quadblock
    ! matrix A right-edge quadblock
    quadblock = matrix_A(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts)
    call essential_boundaries(quadblock, edge='r_edge', matrix='A')
    call natural_boundaries(quadblock, edge='r_edge')
    ! if (hall_mhd) then
    !   call hall_boundaries(quadblock_Hall, kappa_perp_is_zero, edge='r_edge')
    !   hf = hall_field % hallfactor(gauss_gridpts)
    !   quadblock = quadblock - hf * quadblock_Hall
    ! end if
    matrix_A(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts) = quadblock
  end subroutine apply_boundary_conditions


  !> Imposes essential boundary conditions on a given matrix, that is,
  !! the ones that have to be implemented explicitly.
  !! For a wall, that means handling \(v1, a2\) and \(a3\) (and \(T1\) if there is
  !! perpendicular thermal conduction).
  !! @note    The contribution from the 0 quadratic basis function automatically
  !!          zeroes out the odd rows/columns for the quadratic variables on the
  !!          left edge. These are explicitly handled. @endnote
  !! @note    In the case of the Arnoldi solver, setting boundary conditions to
  !!          1e20 throws off the solver such that it returns completely wrong results.
  !!          Setting them to a "normal" value, 0 here, seems to fix it. We now have
  !!          an exact match with the QR/QZ solvers in this case.
  !!          Eigenvalues will be introduced at 0 instead of at 1, but this does not
  !!          affect the final spectrum in any way. @endnote
  !! @warning Throws an error if <tt>boundary_type</tt> is not known,
  !!          or if <tt>edge</tt> is not known.
  subroutine essential_boundaries(quadblock, edge, matrix)
    use mod_global_variables, only: boundary_type, solver, viscosity, geometry, &
                                    coaxial
    use mod_logging, only: log_message

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(inout)    :: quadblock(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge
    !> the matrix, either "A" or "B"
    character, intent(in)         :: matrix

    complex(dp)                   :: diagonal_factor
    integer                       :: i, j, qua_zeroes(5), wall_idx_left(4), wall_idx_right(4)
    integer                       :: noslip_idx_left(2), noslip_idx_right(2)

    if (matrix == 'B') then
      diagonal_factor = (1.0d0, 0.0d0)
    else if (matrix == 'A') then
      if (solver == "arnoldi") then
        diagonal_factor = (0.0d0, 0.0d0)
      else
        diagonal_factor = (1.0d20, 0.0d0)
      end if
    else
      call log_message("essential boundaries: invalid matrix argument", level='error')
    end if

    ! explicitly handle zeroed out rows/cols due to 0 term quadratic basis function
    qua_zeroes = [1, 5, 7, 9, 11]
    if (edge == 'l_edge') then
      do i = 1, size(qua_zeroes)
        j = qua_zeroes(i)
        quadblock(j, j) = diagonal_factor
      end do
    end if

    ! wall/regularity conditions: handling of v1, a2 and a3 (and T if conduction).
    ! v1, a2 and a3 are cubic elements, so omit non-zero basis functions (odd rows/columns)
    ! T is a quadratic element, so omit even row/columns
    wall_idx_left = [3, 13, 15, 10]
    wall_idx_right = [19, 29, 31, 26]
    ! For a no-slip condition for viscous fluids v2 and v3 should equal the wall's
    ! tangential velocities, here zero. (Quadratic elements, so omit even row/columns)
    noslip_idx_left = [6, 8]
    noslip_idx_right = [22, 24]

    select case(boundary_type)
    case('wall')
      if (edge == 'l_edge') then
        ! left regularity/wall conditions
        do i = 1, size(wall_idx_left)
          j = wall_idx_left(i)
          if (j == 10 .and. kappa_perp_is_zero) then
            cycle
          end if
          quadblock(j, :) = (0.0d0, 0.0d0)
          quadblock(:, j) = (0.0d0, 0.0d0)
          quadblock(j, j) = diagonal_factor
        end do
        ! No-slip condition does not apply at left edge for a cylindrical
        ! geometry unless two coaxial walls are used
        if (viscosity .and. (geometry == 'Cartesian' .or. coaxial)) then
          do i = 1, size(noslip_idx_left)
            j = noslip_idx_left(i)
            quadblock(j, :) = (0.0d0, 0.0d0)
            quadblock(:, j) = (0.0d0, 0.0d0)
            quadblock(j, j) = diagonal_factor
          end do
        end if
      else if (edge == 'r_edge') then
        do i = 1, size(wall_idx_right)
          j = wall_idx_right(i)
          if ((j == 26) .and. kappa_perp_is_zero) then
            cycle
          end if
          quadblock(j, :) = (0.0d0, 0.0d0)
          quadblock(:, j) = (0.0d0, 0.0d0)
          quadblock(j, j) = diagonal_factor
        end do
        if (viscosity) then
          do i = 1, size(noslip_idx_right)
            j = noslip_idx_right(i)
            quadblock(j, :) = (0.0d0, 0.0d0)
            quadblock(:, j) = (0.0d0, 0.0d0)
            quadblock(j, j) = diagonal_factor
          end do
        end if
      else
        call log_message("essential boundaries: invalid edge argument", level='error')
      end if

    case default
      call log_message( "essential boundaries: invalid boundary_type", level='error')
    end select
  end subroutine essential_boundaries


  !> Imposes the natural boundary conditions, that is, the ones originating
  !! from the weak formulation. This is only applicable for the A-matrix.
  !! @note    For now only the terms of a solid wall boundary are implemented, hence
  !!          we only have additional resistive and conductive terms in the energy equation.
  !!          If perpendicular thermal conduction is included we have \(T1 = 0\), such that
  !!          then there are no surface terms to be added and we simply return.
  !!          If free boundary conditions are considered (e.g. vacuum-wall) additional terms
  !!          have to be added. @endnote
  !! @warning Throws an error if <tt>edge</tt> is unknown.
  subroutine natural_boundaries(quadblock, edge)
    use mod_global_variables, only: boundary_type, gauss_gridpts, gamma_1, ic, dim_subblock
    use mod_logging, only: log_message
    use mod_equilibrium, only: B_field, eta_field
    use mod_grid, only: eps_grid, d_eps_grid_dr
    use mod_equilibrium_params, only: k2, k3

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(inout)    :: quadblock(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge

    complex(dp), allocatable  :: surface_terms(:)
    real(dp)                  :: eps, d_eps_dr, eta, B02, dB02, B03, dB03, drB02
    integer, allocatable      :: positions(:, :)
    integer                   :: idx, i

    if (.not. boundary_type == 'wall') then
      call log_message('natural boundaries: only wall is implemented!', level='error')
    end if

    ! For a wall with perpendicular thermal conduction we also have the condition T1 = 0.
    ! Hence, in that case there are no surface terms to be added, and we simply return.
    ! So for now we only have resistive terms in the energy equation.
    if (.not. kappa_perp_is_zero) then
      return
    end if

    if (edge == 'l_edge') then
      idx = 1
    else if (edge == 'r_edge') then
      idx = gauss_gridpts
    else
      call log_message('natural boundaries: wrong edge supplied' // edge, level='error')
    end if

    ! retrieve variables at current edge
    eps = eps_grid(idx)
    d_eps_dr = d_eps_grid_dr(idx)
    B02 = B_field % B02(idx)
    dB02 = B_field % d_B02_dr(idx)
    B03 = B_field % B03(idx)
    dB03 = B_field % d_B03_dr(idx)
    drB02 = d_eps_dr * B02 + eps * dB02
    eta = eta_field % eta(idx)

    allocate(positions(3, 2))
    allocate(surface_terms(3))

    ! surface term for element (5, 6)
    surface_terms(1) = 2.0d0 * ic * gamma_1 * (1.0d0 / eps) * eta * (k3 * drB02 - k2 * dB03)
    positions(1, :) = [5, 6]
    ! surface term for element (5, 7)
    surface_terms(2) = 2.0d0 * ic * gamma_1 * (1.0d0 / eps) * eta * dB03
    positions(2, :) = [5, 7]
    ! surface term for element (5, 8)
    surface_terms(3) = -2.0d0 * ic * gamma_1 * (1.0d0 / eps) * eta * drB02
    positions(3, :) = [5, 8]

    ! l_edge: add to bottom-right of 2x2 block, for top-left subblock only
    ! r_edge: add to bottom-right of 2x2 block, for bottom-right subblock only
    if (edge == 'l_edge') then
      positions = 2 * positions
    else if (edge == 'r_edge') then
      positions = 2 * positions + dim_subblock
    end if

    do i = 1, size(surface_terms)
      quadblock(positions(i, 1), positions(i, 2)) = quadblock(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do

    deallocate(positions)
    deallocate(surface_terms)
  end subroutine natural_boundaries

end module mod_boundary_conditions
