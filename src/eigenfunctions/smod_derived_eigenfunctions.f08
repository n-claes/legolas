! =============================================================================
!> This module composes perturbed quantities of interest based on the
!! eigenfunctions and eigenvectors. As a result, all quantities are defined on
!! the same grid as the eigenfunctions, ef_grid.
!! The computed quantities are
!!    - div v1
!!    - all 3 components of curl v1
!!    - all 3 components of B1
!!    - div B1
!!    - all 3 components of curl B1
!!    - 3 components of v1, parallel/perpendicular to B0
!!    - 2 components of curl v1, parallel and perpendicular to B0
!! It is assumed that B0 has no B01 component, such that the 1-direction is
!! already perpendicular to B0, and the parallel/perpendicular directions are
!! taken in the 23-plane. The right-handed triad is e1, B0/|B0|, (e1 x B0) / |B0|.

submodule(mod_eigenfunctions) smod_derived_eigenfunctions
  use mod_global_variables, only: ic
  use mod_equilibrium_params, only: k2, k3
  implicit none

  !> Determines if parallel/perpendicular quantities can be computed
  logical, save :: can_calculate_pp_quantities
  !> array containing the interpolated B02
  real(dp), allocatable :: B02_on_ef_grid(:)
  !> array containing the interpolated B03
  real(dp), allocatable :: B03_on_ef_grid(:)
  !> array containing the interpolated rho0
  real(dp), allocatable :: rho0_on_ef_grid(:)
  !> array containing the interpolated T0
  real(dp), allocatable :: T0_on_ef_grid(:)

contains

  !> Determines if parallel/perpendicular quantities can be calculated and
  !! sets the corresponding module flag. This flag will be false if there is a
  !! non-zero B01 component present, or no magnetic field.
  subroutine check_if_para_perp_quantities_can_be_calculated()
    use mod_check_values, only: is_zero
    use mod_equilibrium, only: B_field
    use mod_logging, only: log_message

    can_calculate_pp_quantities = .true.
    if (.not. is_zero(B_field % B01)) then
      call log_message( &
        "parallel/perpendicular derived quantities currently not supported &
        &for non-zero B01 components", level="warning" &
      )
      can_calculate_pp_quantities = .false.
    end if
    if (all(is_zero(B_field % B02)) .and. all(is_zero(B_field % B03))) then
      can_calculate_pp_quantities = .false.
    end if
  end subroutine check_if_para_perp_quantities_can_be_calculated


  !> Interpolates the B02, B03, rho0 and T0 equilibrium components, used for
  !! calculating the derived eigenfunction quantities.
  subroutine interpolate_equilibrium_on_eigenfunction_grid()
    use mod_global_variables, only: gridpts
    use mod_interpolation, only: interpolate_table, lookup_table_value
    use mod_grid, only: grid_gauss
    use mod_equilibrium, only: B_field, rho_field, T_field

    call interpolate_and_set_values(B_field % B02, B02_on_ef_grid)
    call interpolate_and_set_values(B_field % B03, B03_on_ef_grid)
    call interpolate_and_set_values(rho_field % rho0, rho0_on_ef_grid)
    call interpolate_and_set_values(T_field % T0, T0_on_ef_grid)

    contains
      !> interpolates an equilibrium array and sets it on the eigenfunction grid.
      subroutine interpolate_and_set_values(basic_array, expected_array)
        use mod_interpolation, only: interpolate_table, lookup_table_value
        use mod_grid, only: grid_gauss

        !> basic array (equilibrium)
        real(dp), intent(in)  :: basic_array
        !> expected array (set on eigenfunction grid)
        real(dp), allocatable, intent(out) :: expected_array(:)
        !> number of points used for interpolation
        integer :: ip_pts
        !> interpolated grid
        real(dp), allocatable  :: interpolated_grid(:)
        !> interpolated equilibrium values
        real(dp), allocatable  :: interpolated_values(:)

        ip_pts = 25 * (gridpts - 1) + 1
        allocate(interpolated_grid(ip_pts))
        allocate(interpolated_values(ip_pts))
        call interpolate_table(ip_pts, grid_gauss, basic_array, interpolated_values)

        allocate(expected_array(size(ef_grid)))
        do i = 1, size(ef_grid)
          expected_array(i) = lookup_table_value( &
            ef_grid(i), interpolated_grid, interpolated_values &
          )
        end do
        deallocate(interpolated_grid)
        deallocate(interpolated_values)
      end subroutine interpolate_and_set_values
  end subroutine interpolate_equilibrium_on_eigenfunction_grid


  !> Initialised the derived eigenfunction array, sets the corresponding names and
  !! vector indices, allocates the (subset of) derived eigenfunctions
  module procedure initialise_derived_eigenfunctions
    use mod_equilibrium, only: B_field
    use mod_check_values, only: is_equal

    integer :: i

    if (can_calculate_pp_quantities) then
      derived_ef_names = [ &
        character(len=str_len_arr) :: "S", "div v", "(curl v)1", "(curl v)2", &
        "(curl v)3", "B1", "B2", "B3", "div B", &
        "(curl B)1", "(curl B)2", "(curl B)3", &
        "B_para", "B_perp", "(curl B)_para", &
        "(curl B)_perp", "v1", "v_para", "v_perp", &
        "(curl v)_para", "(curl v)_perp" &
      ]
    else
      derived_ef_names = [ &
        character(len=str_len_arr) :: "S", "div v", "(curl v)1", "(curl v)2", &
        "(curl v)3", "B1", "B2", "B3", "div B", &
        "(curl B)1", "(curl B)2", "(curl B)3" &
      ]
    end if
    allocate(derived_eigenfunctions(size(derived_ef_names)))

    do i = 1, size(derived_eigenfunctions)
      derived_eigenfunctions(i) % state_vector_index = i
      derived_eigenfunctions(i) % name = derived_ef_names(i)
      allocate(derived_eigenfunctions(i)%quantities(size(ef_grid), nb_eigenfuncs))
    end do

    call interpolate_equilibrium_on_eigenfunction_grid()
    derived_efs_initialised = .true.
  end procedure initialise_derived_eigenfunctions


  !> Calculates the derived eigenfunction quantities corresponding to the requested
  !! eigenvalues and sets them as attributes for the corresponding types.
  module procedure calculate_derived_eigenfunctions
    !> quantities values
    complex(dp) :: single(size(ef_grid))
    complex(dp) :: double(size(ef_grid), 2)
    complex(dp) :: triple(size(ef_grid), 3)
    integer     :: i, j

    call set_entropy(derived_eigenfunctions(1))
    call set_velocity_divergence(derived_eigenfunctions(2), right_eigenvectors)
    call set_velocity_curl_1_component(derived_eigenfunctions(3))
    call set_velocity_curl_2_component(derived_eigenfunctions(4), right_eigenvectors)
    call set_velocity_curl_3_component(derived_eigenfunctions(5), right_eigenvectors)

    if (can_calculate_pp_quantities) then
      call set_velocity_curl_parallel_component(derived_eigenfunctions(20))
      call set_velocity_curl_perpendicular_component(derived_eigenfunctions(21))
    end if

    !! Magnetic field B1
    do i = 1, size(vr, dim=2)
      call get_magnetic_field(i, vr(:, i), triple, double)
      derived_eigenfunctions(6) % quantities(:, i) = triple(:, 1)
      derived_eigenfunctions(7) % quantities(:, i) = triple(:, 2)
      derived_eigenfunctions(8) % quantities(:, i) = triple(:, 3)
      !! In parallel/perpendicular components w.r.t B0 (no B01, so set becomes B1, Bpara, Bperp)
      if (.not. B_zero) then
        derived_eigenfunctions(13) % quantities(:, i) = double(:, 1)
        derived_eigenfunctions(14) % quantities(:, i) = double(:, 2)
      end if
    end do

    !! Divergence of B1
    do i = 1, size(vr, dim=2)
      call get_div_magnetic(i, vr(:, i), single)
      derived_eigenfunctions(9) % quantities(:, i) = single
    end do

    !! Curl of B1
    do i = 1, size(vr, dim=2)
      call get_curl_magnetic(i, vr(:, i), triple, double)
      derived_eigenfunctions(10) % quantities(:, i) = triple(:, 1)
      derived_eigenfunctions(11) % quantities(:, i) = triple(:, 2)
      derived_eigenfunctions(12) % quantities(:, i) = triple(:, 3)
      !! In parallel/perpendicular components w.r.t B0 (no B01, so set becomes
      !! (curlB1)1, (curlB1)para, (curlB1)perp)
      if (.not. B_zero) then
        derived_eigenfunctions(15) % quantities(:, i) = double(:, 1)
        derived_eigenfunctions(16) % quantities(:, i) = double(:, 2)
      end if
    end do

    !! v1 components in parallel/perpendicular w.r.t. B0 (no B01, so set is
    !! v1, v_para, v_perp)
    if (.not. B_zero) then
      do i = 1, size(vr, dim=2)
        call get_v_paraperp(i, triple)
        derived_eigenfunctions(17) % quantities(:, i) = triple(:, 1)
        derived_eigenfunctions(18) % quantities(:, i) = triple(:, 2)
        derived_eigenfunctions(19) % quantities(:, i) = triple(:, 3)
      end do
    end if
  end procedure calculate_derived_eigenfunctions


  !> Calculates the entropy S1.
  subroutine set_entropy(derived_ef)
    !> derived eigenfunction type at entropy position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> density eigenfunction
    complex(dp)  :: rho_ef(size(ef_grid))
    !> temperature eigenfunction
    complex(dp)  :: T_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      rho_ef = base_eigenfunctions(findloc(ef_names, "rho"))%quantities(:, i)
      T_ef = base_eigenfunctions(findloc(ef_names, "T"))%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        T_ef / rho0_on_ef_grid ** (2.0d0 / 3.0d0) &
        - (2.0d0/3.0d0) * rho_ef * T0_on_ef_grid / rho0_on_ef_grid ** (5.0d0/3.0d0) &
      )
    end do
  end subroutine set_entropy


  !> Calculates the divergence of the perturbed velocity v1
  subroutine set_velocity_divergence(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at div(v) position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> v2 eigenfunction
    complex(dp) :: v2_ef(size(ef_grid))
    !> v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    !> derivative of v1 eigenfunction
    complex(dp) :: dv1_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2"))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3"))%quantities(:, i)
      dv1_ef = get_assembled_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "v1")), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        ic * (-dv1_ef / ef_eps + k2 * v2_ef / eps + k3 * v3_ef) &
      )
    end do
  end subroutine set_velocity_divergence


  !> Calculates the 1-component of the curl of the perturbed velocity v1
  subroutine set_velocity_curl_1_component(derived_ef)
    !> derived eigenfunction type at curl(v)1 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    ! v1 eigenfunction
    complex(dp) :: v1_ef(size(ef_grid))
    ! v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      v1_ef = base_eigenfunctions(findloc(ef_names, "v1"))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3"))%quantities(:, i)
      derived_ef%quantities(:, i) = ic * (k2 * v3_ef / ef_eps - k3 * v2_ef)
    end do
  end subroutine set_velocity_curl_1_component


  !> Calculates the 2-component of the curl of the perturbed velocity v1
  subroutine set_velocity_curl_2_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at curl(v)2 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    ! v1 eigenfunction
    complex(dp) :: v1_ef(size(ef_grid))
    ! v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    !> dv3 eigenfunction
    complex(dp) :: dv3_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      v1_ef = base_eigenfunctions(findloc(ef_names, "v1"))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3"))%quantities(:, i)
      dv3_ef = get_assembled_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "v3")), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        ic * k3 * v1_ef - dv3_ef / ef_eps - ef_deps * v3_ef / ef_eps &
      )
    end do
  end subroutine set_velocity_curl_2_component


  !> Calculates the 3-component of the curl of the perturbed velocity v1
  subroutine set_velocity_curl_3_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at curl(v)3 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    ! v1 eigenfunction
    complex(dp) :: v1_ef(size(ef_grid))
    ! v2 eigenfunction
    complex(dp) :: v2_ef(size(ef_grid))
    !> dv2 eigenfunction
    complex(dp) :: dv2_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      v1_ef = base_eigenfunctions(findloc(ef_names, "v1"))%quantities(:, i)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2"))%quantities(:, i)
      dv2_ef = get_assembled_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "v2")), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        dv2_ef(m) + (ef_deps * v2_ef - ic * k2 * v1_ef) / ef_eps &
      )
    end do
  end subroutine set_velocity_curl_3_component


  !> Sets the velocity curl component parallel to the magnetic field.
  subroutine set_velocity_curl_parallel_component(derived_ef)
    !> derived eigenfunction type at curl(v)_parallel position in the array
    type(ef_type), intent(inout)  :: derived_ef
    complex(dp) :: curlv2(size(ef_grid))
    complex(dp) :: curlv3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlv2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)2") &
      )%quantities(:, i)
      curlv3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)3") &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * curlv2 + B03_on_ef_grid * curlv3 &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_velocity_curl_parallel_component


  !> Sets the velocity curl component perpendicular to the magnetic field
  subroutine set_velocity_curl_perpendicular_component(derived_ef)
    !> derived eigenfunction type at curl(v)_perpendicular position in the array
    type(ef_type), intent(inout)  :: derived_ef
    complex(dp) :: curlv2(size(ef_grid))
    complex(dp) :: curlv3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlv2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)2") &
      )%quantities(:, i)
      curlv3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)3") &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * curlv3 - B03_on_ef_grid * curlv2 &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_velocity_curl_perpendicular_component


  !> Calculates the perturbed magnetic field B1
  subroutine get_magnetic_field(index, eigenvector, field, field2)
    use mod_global_variables, only: gauss_gridpts
    use mod_equilibrium, only: B_field

    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: field(ef_gridpts, 3), field2(ef_gridpts, 2)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da2(ef_gridpts), da3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps, x, B02, B03, B0
    complex(dp) :: B2, B3

    a1 = base_eigenfunctions(6) % eigenfunctions(:, index)
    a2 = base_eigenfunctions(7) % eigenfunctions(:, index)
    a3 = base_eigenfunctions(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      B2 = ic * k3 * a1(m) - da3(m)
      B3 = (da2(m) - ic * k2 * a1(m)) / eps

      field(m, 1) = ic * (k2 * a3(m) / eps - k3 * a2(m))
      field(m, 2) = B2
      field(m, 3) = B3

      !! In parallel/perpendicular components
      if (.not. B_zero) then
        x = ef_grid(m)
        if (m == 1) then
          B02 = B_field % B02(1)
          B03 = B_field % B03(1)
        else if (m == ef_gridpts) then
          B02 = B_field % B02(gauss_gridpts)
          B03 = B_field % B03(gauss_gridpts)
        else
          B02 = lookup_table_value(x, grid_ip(:, 1), eq_ip(:, 1))
          B03 = lookup_table_value(x, grid_ip(:, 2), eq_ip(:, 2))
        end if
        B0 = sqrt(B02**2 + B03**2)

        field2(m, 1) = (B02 * B2 + B03 * B3) / B0
        field2(m, 2) = (B02 * B3 - B03 * B2) / B0
      end if
    end do
  end subroutine get_magnetic_field


  !> Calculates the divergence of the perturbed magnetic field B1
  subroutine get_div_magnetic(index, eigenvector, divergence)
    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: divergence(ef_gridpts)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da2(ef_gridpts), da3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps

    a1 = base_eigenfunctions(6) % eigenfunctions(:, index)
    a2 = base_eigenfunctions(7) % eigenfunctions(:, index)
    a3 = base_eigenfunctions(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      divergence(m) = ic * (k2 * da3(m) - k3 * da2(m)) / eps &
                      + ic * k2 * (ic * k3 * a1(m) - da3(m)) / eps &
                      + ic * k3 * (da2(m) - ic * k2 * a1(m)) / eps
    end do
  end subroutine get_div_magnetic


  !> Calculates the curl of the perturbed magnetic field B1
  subroutine get_curl_magnetic(index, eigenvector, curl, curl2)
    use mod_global_variables, only: gauss_gridpts
    use mod_equilibrium, only: B_field

    integer, intent(in)       :: index
    complex(dp), intent(in)   :: eigenvector(matrix_gridpts)
    complex(dp), intent(out)  :: curl(ef_gridpts, 3), curl2(ef_gridpts, 2)

    !! Derivatives are transformed variables, regular variables are 'actual' values
    complex(dp) :: a1(ef_gridpts), a2(ef_gridpts), a3(ef_gridpts)
    complex(dp) :: da1(ef_gridpts), da2(ef_gridpts), da3(ef_gridpts)
    complex(dp) :: dda2(ef_gridpts), dda3(ef_gridpts)
    integer     :: m
    real(dp)    :: eps, deps, x, B02, B03, B0
    complex(dp) :: curlB2, curlB3

    a1 = base_eigenfunctions(6) % eigenfunctions(:, index)
    a2 = base_eigenfunctions(7) % eigenfunctions(:, index)
    a3 = base_eigenfunctions(8) % eigenfunctions(:, index)
    call get_eigenfunction_deriv(6, eigenvector, da1, 1)
    call get_eigenfunction_deriv(7, eigenvector, da2, 1)
    call get_eigenfunction_deriv(8, eigenvector, da3, 1)
    call get_eigenfunction_deriv(7, eigenvector, dda2, 2)
    call get_eigenfunction_deriv(8, eigenvector, dda3, 2)

    do m = 1, ef_gridpts
      call set_eps(m, eps, deps)
      curlB2 = -k3 * (k2 * a3(m) / eps - k3 * a2(m)) &
                   - (dda2(m) - k2 * da1(m)) / eps &
                   + deps * (da2(m) - ic * k2 * a1(m)) / eps**2
      curlB3 = k3 * da1(m) - dda3(m) + deps * (ic * k3 * a1(m) - da3(m)) / eps &
                   + k2 * (k2 * a3(m) / eps - k3 * a2(m)) / eps

      curl(m, 1) = ic * k2 * (da2(m) - ic * k2 * a1(m)) / eps**2 &
                   - ic * k3 * (ic * k3 * a1(m) - da3(m))
      curl(m, 2) = curlB2
      curl(m, 3) = curlB3

     !! In parallel/perpendicular components
     if (.not. B_zero) then
       x = ef_grid(m)
       if (m == 1) then
         B02 = B_field % B02(1)
         B03 = B_field % B03(1)
       else if (m == ef_gridpts) then
         B02 = B_field % B02(gauss_gridpts)
         B03 = B_field % B03(gauss_gridpts)
       else
         B02 = lookup_table_value(x, grid_ip(:, 1), eq_ip(:, 1))
         B03 = lookup_table_value(x, grid_ip(:, 2), eq_ip(:, 2))
       end if
       B0 = sqrt(B02**2 + B03**2)

       curl2(m, 1) = (B02 * curlB2 + B03 * curlB3) / B0
       curl2(m, 2) = (B02 * curlB3 - B03 * curlB2) / B0
     end if
    end do
  end subroutine get_curl_magnetic


  subroutine get_v_paraperp(index, field)
    use mod_global_variables, only: gauss_gridpts
    use mod_equilibrium, only: B_field

    integer, intent(in)       :: index
    complex(dp), intent(out)  :: field(ef_gridpts, 3)

    complex(dp) :: v2(ef_gridpts), v3(ef_gridpts)
    integer     :: m
    real(dp)    :: x, B02, B03, B0

    field(:, 1) = base_eigenfunctions(2) % eigenfunctions(:, index)
    v2 = base_eigenfunctions(3) % eigenfunctions(:, index)
    v3 = base_eigenfunctions(4) % eigenfunctions(:, index)

    do m = 1, ef_gridpts
      !! In parallel/perpendicular components
      x = ef_grid(m)
      if (m == 1) then
        B02 = B_field % B02(1)
        B03 = B_field % B03(1)
      else if (m == ef_gridpts) then
        B02 = B_field % B02(gauss_gridpts)
        B03 = B_field % B03(gauss_gridpts)
      else
        B02 = lookup_table_value(x, grid_ip(:, 1), eq_ip(:, 1))
        B03 = lookup_table_value(x, grid_ip(:, 2), eq_ip(:, 2))
      end if
      B0 = sqrt(B02**2 + B03**2)

      field(m, 2) = (B02 * v2(m) + B03 * v3(m)) / B0
      field(m, 3) = (B02 * v3(m) - B03 * v2(m)) / B0
    end do
  end subroutine get_v_paraperp


  !> Calculates the derivative of a single eigenfunction based on an eigenvector.
  !! The eigenvector supplied is the one corresponding to one particular
  !! eigenvalue.
  subroutine get_eigenfunction_deriv(ef_idx, eigenvector, Y, deriv)
    use mod_global_variables, only: dim_subblock, gridpts
    use mod_grid, only: grid
    use mod_logging, only: log_message

    !> the index of the variable in the state vector
    integer, intent(in)           :: ef_idx
    !> the eigenvector for this particular eigenvalue
    complex(dp), intent(in)       :: eigenvector(matrix_gridpts)
    !> the assembled eigenfunction
    complex(dp), intent(out)      :: Y(ef_gridpts)

    integer                       :: idx1, idx2, grid_idx, i
    real(dp)                      :: r, r_lo, r_hi
    real(dp)                      :: spline(4)
    !> the deriv-th spline derivative, defaults to 0
    integer, optional             :: deriv

    if (.not. present(deriv)) then
      deriv = 0
    end if

    ! initialise eigenfunction to zero
    Y = (0.0d0, 0.0d0)
    ! map ef_idx to actual subblock index. So 1 -> 1, 2 -> 3, 3 -> 5 etc.
    idx1 = 2 * ef_idx - 1
    idx2 = idx1 + 1

    ! Contribution from first gridpoint, left edge
    r_lo = grid(1)
    r_hi = grid(2)
    r    = grid(1)
    grid_idx = 1
    ef_grid(grid_idx) = r
    call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
    Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                              + eigenvector(idx2) * spline(4) &
                              + eigenvector(idx1 + dim_subblock) * spline(1) &
                              + eigenvector(idx2 + dim_subblock) * spline(3)

    ! Contribution from other gridpoints
    do i = 1, gridpts-1
      ! Contribution from centre
      r_lo = grid(i)
      r_hi = grid(i + 1)
      r    = 0.5 * (r_lo + r_hi)
      ! save gridpoint
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r

      call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Contribution from end point
      r = r_hi
      grid_idx = grid_idx + 1
      ef_grid(grid_idx) = r
      call get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
      Y(grid_idx) = Y(grid_idx) + eigenvector(idx1) * spline(2) &
                                + eigenvector(idx2) * spline(4) &
                                + eigenvector(idx1 + dim_subblock) * spline(1) &
                                + eigenvector(idx2 + dim_subblock) * spline(3)

      ! Increment indices of eigenvector elements
      idx1 = idx1 + dim_subblock
      idx2 = idx2 + dim_subblock
    end do
  end subroutine get_eigenfunction_deriv


  !> Retrieves the spline or one of its derivatives for a particular variable by
  !! calling the correct finite element basis function depending on the
  !! variable passed.
  subroutine get_spline_deriv(ef_idx, r, r_lo, r_hi, spline, deriv)
    use mod_spline_functions, only: cubic_factors, quadratic_factors, &
                                    cubic_factors_deriv, quadratic_factors_deriv, &
                                    cubic_factors_deriv2

    !> the index of the variable in the state vector
    integer, intent(in)     :: ef_idx
    !> the current position in the grid
    real(dp), intent(in)    :: r
    !> left edge of the current grid interval
    real(dp), intent(in)    :: r_lo
    !> right edge of the current grid interval
    real(dp), intent(in)    :: r_hi
    !> the finite element basis functions for this grid position
    real(dp), intent(out)   :: spline(4)
    !> the deriv-th spline derivative, defaults to 0 (no derivative)
    integer, optional       :: deriv

    if (.not. present(deriv)) then
      deriv = 0
    end if

    if (ef_idx == 2 .or. ef_idx == 7 .or. ef_idx == 8) then
      if (deriv == 2) then
        call cubic_factors_deriv2(r, r_lo, r_hi, spline)
      else if (deriv == 1) then
        call cubic_factors_deriv(r, r_lo, r_hi, spline)
      else
        call cubic_factors(r, r_lo, r_hi, spline)
      end if
    else
      if (deriv == 1) then
        call quadratic_factors_deriv(r, r_lo, r_hi, spline)
      else
        call quadratic_factors(r, r_lo, r_hi, spline)
      end if
    end if
  end subroutine get_spline_deriv


  subroutine set_eps(index, eps, deps)
    use mod_global_variables, only: geometry

    integer, intent(in)   :: index
    real(dp), intent(out) :: eps, deps

    if (geometry == 'cylindrical') then
      eps  = ef_grid(index)
      deps = 1.0d0
    else
      eps  = 1.0d0
      deps = 1.0d0
    end if
  end subroutine set_eps


  !> Cleaning routine, deallocated the derived quantities
  module procedure clean_derived_eigenfunctions()
    deallocate(B02_on_ef_grid)
    deallocate(B03_on_ef_grid)
    deallocate(rho0_on_ef_grid)
    deallocate(T0_on_ef_grid)
  end procedure clean_derived_eigenfunctions

end submodule smod_derived_eigenfunctions
