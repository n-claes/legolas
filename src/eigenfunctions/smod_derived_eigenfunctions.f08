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
  subroutine check_if_pp_quantities_can_be_calculated()
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
  end subroutine check_if_pp_quantities_can_be_calculated


  !> Sets the B02, B03, rho0 and T0 equilibrium components on the eigenfunction grid
  !! instead of on the Gaussian grid, these are used for calculating the derived
  !! eigenfunction quantities. No interpolation needed since the Gaussian grid is
  !! always finer than the eigenfunction grid, so we can get away by simply
  !! looking up the nearest values for every point in the eigenfunction grid.
  subroutine set_equilibrium_on_eigenfunction_grid()
    use mod_interpolation, only: interpolate_table, lookup_table_value
    use mod_equilibrium, only: B_field, rho_field, T_field

    call lookup_and_set_values(B_field % B02, B02_on_ef_grid)
    call lookup_and_set_values(B_field % B03, B03_on_ef_grid)
    call lookup_and_set_values(rho_field % rho0, rho0_on_ef_grid)
    call lookup_and_set_values(T_field % T0, T0_on_ef_grid)

    contains
      !> Sets an equilibrium array on the eigenfunction grid
      subroutine lookup_and_set_values(equil_array, equil_array_on_ef_grid)
        use mod_interpolation, only: interpolate_table, lookup_table_value
        use mod_grid, only: grid_gauss

        !> basic equilibrium array
        real(dp), intent(in)  :: equil_array(:)
        !> equilibrium array set on eigenfunction grid
        real(dp), allocatable, intent(out) :: equil_array_on_ef_grid(:)
        integer :: i

        allocate(equil_array_on_ef_grid(size(ef_grid)))
        do i = 1, size(ef_grid)
          ! allow outside since edges of eigenfunction grid are beyond Gaussian grid
          equil_array_on_ef_grid(i) = lookup_table_value( &
            ef_grid(i), grid_gauss, equil_array, allow_outside=.true. &
          )
        end do
      end subroutine lookup_and_set_values
  end subroutine set_equilibrium_on_eigenfunction_grid


  !> Initialised the derived eigenfunction array, sets the corresponding names and
  !! vector indices, allocates the (subset of) derived eigenfunctions
  module procedure initialise_derived_eigenfunctions
    use mod_check_values, only: is_equal

    integer :: i

    call check_if_pp_quantities_can_be_calculated()
    if (can_calculate_pp_quantities) then
      derived_ef_names = [ &
        character(len=str_len_arr) :: "S", "div v", "(curl v)1", "(curl v)2", &
        "(curl v)3", "B1", "B2", "B3", "div B", &
        "(curl B)1", "(curl B)2", "(curl B)3", &
        "B_para", "B_perp", "(curl B)_para", &
        "(curl B)_perp", "v_para", "v_perp", &
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

    call set_equilibrium_on_eigenfunction_grid()
    derived_efs_initialised = .true.
  end procedure initialise_derived_eigenfunctions


  !> Calculates the derived eigenfunction quantities corresponding to the requested
  !! eigenvalues and sets them as attributes for the corresponding types.
  module procedure calculate_derived_eigenfunctions
    integer :: i, eigenvalue_idx
    ! the density and temperature eigenfunctions
    type(ef_type) :: rho_efs, T_efs
    ! the v1, v2 and v3 eigenfunctions
    type(ef_type) :: v1_efs, v2_efs, v3_efs

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)

      rho_efs = retrieve_eigenfunction(name="rho")
      v1_efs = retrieve_eigenfunction(name="v1")
      v2_efs = retrieve_eigenfunction(name="v2")
      v3_efs = retrieve_eigenfunction(name="v3")
      T_efs = retrieve_eigenfunction(name="T")

      call set_entropy( &
        loc=1, &
        ef_index=i, &
        rho_ef=rho_efs%quantities(:, i), &
        T_ef=T_efs%quantities(:, i) &
      )
      call set_velocity_divergence( &
        loc=2, &
        ef_index=i, &
        v2_ef=v2_efs%quantities(:, i), &
        v3_ef=v3_efs%quantities(:, i), &
        rvec=right_eigenvectors(:, eigenvalue_idx) &
      )
    end do

    call set_velocity_curl_1_component(derived_eigenfunctions(3))
    call set_velocity_curl_2_component(derived_eigenfunctions(4), right_eigenvectors)
    call set_velocity_curl_3_component(derived_eigenfunctions(5), right_eigenvectors)
    call set_magnetic_field_ef_1_component(derived_eigenfunctions(6))
    call set_magnetic_field_ef_2_component( &
      derived_eigenfunctions(7), right_eigenvectors &
    )
    call set_magnetic_field_ef_3_component( &
      derived_eigenfunctions(8), right_eigenvectors &
    )
    call set_magnetic_field_divergence(derived_eigenfunctions(9), right_eigenvectors)
    call set_magnetic_field_curl_1_component( &
      derived_eigenfunctions(10), right_eigenvectors &
    )
    call set_magnetic_field_curl_2_component( &
      derived_eigenfunctions(11), right_eigenvectors &
    )
    call set_magnetic_field_curl_3_component( &
      derived_eigenfunctions(12), right_eigenvectors &
    )

    if (can_calculate_pp_quantities) then
      call set_magnetic_field_parallel_component(derived_eigenfunctions(13))
      call set_magnetic_field_perpendicular_component(derived_eigenfunctions(14))
      call set_magnetic_field_parallel_curl_component(derived_eigenfunctions(15))
      call set_magnetic_field_perpendicular_curl_component(derived_eigenfunctions(16))
      call set_velocity_parallel_component(derived_eigenfunctions(17))
      call set_velocity_perpendicular_component(derived_eigenfunctions(18))
      call set_velocity_curl_parallel_component(derived_eigenfunctions(19))
      call set_velocity_curl_perpendicular_component(derived_eigenfunctions(20))
    end if
  end procedure calculate_derived_eigenfunctions


  !> Calculates the entropy S1 and places it at location "loc" in the main array.
  subroutine set_entropy(loc, ef_index, rho_ef, T_ef)
    !> position index in the main array
    integer, intent(in) :: loc
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index
    !> density eigenfunction
    complex(dp), intent(in)  :: rho_ef(:)
    !> temperature eigenfunction
    complex(dp), intent(in)  :: T_ef(:)

    derived_eigenfunctions(loc)%quantities(:, ef_index) = ( &
      T_ef / rho0_on_ef_grid ** (2.0d0 / 3.0d0) &
      - (2.0d0/3.0d0) * (rho_ef * T0_on_ef_grid / rho0_on_ef_grid ** (5.0d0/3.0d0)) &
    )
  end subroutine set_entropy


  !> Calculates the divergence of the perturbed velocity v1
  subroutine set_velocity_divergence(loc, ef_index, v2_ef, v3_ef, rvec)
    !> position index in the main array
    integer, intent(in) :: loc
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index
    !> v2 eigenfunction
    complex(dp), intent(in) :: v2_ef(:)
    !> v3 eigenfunction
    complex(dp), intent(in) :: v3_ef(:)
    !> right eigenvector
    complex(dp), intent(in) :: rvec(:)
    !> derivative of v1 eigenfunction
    complex(dp) :: dv1_ef(size(ef_grid))

    dv1_ef = assemble_eigenfunction( &
      base_ef=retrieve_eigenfunction(name="v1"), eigenvector=rvec, derivative_order=1 &
    )
    derived_eigenfunctions(loc)%quantities(:, ef_index) = ( &
      ic * (-dv1_ef / ef_eps + k2 * v2_ef / ef_eps + k3 * v3_ef) &
    )
  end subroutine set_velocity_divergence


  !> Calculates the 1-component of the curl of the perturbed velocity v1
  subroutine set_velocity_curl_1_component(derived_ef)
    !> derived eigenfunction type at curl(v)1 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    ! v2 eigenfunction
    complex(dp) :: v2_ef(size(ef_grid))
    ! v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2", dim=1))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3", dim=1))%quantities(:, i)
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
      v1_ef = base_eigenfunctions(findloc(ef_names, "v1", dim=1))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3", dim=1))%quantities(:, i)
      dv3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "v3", dim=1)), &
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
      v1_ef = base_eigenfunctions(findloc(ef_names, "v1", dim=1))%quantities(:, i)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2", dim=1))%quantities(:, i)
      dv2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "v2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        dv2_ef + (ef_deps * v2_ef - ic * k2 * v1_ef) / ef_eps &
      )
    end do
  end subroutine set_velocity_curl_3_component


  !> Sets the velocity curl component parallel to the magnetic field.
  subroutine set_velocity_curl_parallel_component(derived_ef)
    !> derived eigenfunction type at curl(v)_parallel position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> curl v2 eigenfunction
    complex(dp) :: curlv2(size(ef_grid))
    !> curl v3 eigenfunction
    complex(dp) :: curlv3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlv2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)2", dim=1) &
      )%quantities(:, i)
      curlv3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)3", dim=1) &
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
    !> curl v2 eigenfunction
    complex(dp) :: curlv2(size(ef_grid))
    !> curl v3 eigenfunction
    complex(dp) :: curlv3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlv2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)2", dim=1) &
      )%quantities(:, i)
      curlv3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl v)3", dim=1) &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * curlv3 - B03_on_ef_grid * curlv2 &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_velocity_curl_perpendicular_component


  !> Calculates the 1-component of the perturbed magnetic field
  subroutine set_magnetic_field_ef_1_component(derived_ef)
    !> derived eigenfunction type at B1-1 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> a2 eigenfunction
    complex(dp) :: a2_ef(size(ef_grid))
    !> a3 eigenfunction
    complex(dp) :: a3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      a2_ef = base_eigenfunctions(findloc(ef_names, "a2", dim=1))%quantities(:, i)
      a3_ef = base_eigenfunctions(findloc(ef_names, "a3", dim=1))%quantities(:, i)
      derived_ef%quantities(:, i) = ic * (k2 * a3_ef / ef_eps - k3 * a2_ef)
    end do
  end subroutine set_magnetic_field_ef_1_component


  !> Calculates the 2-component of the perturbed magnetic field
  subroutine set_magnetic_field_ef_2_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at B1-2 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> da3 eigenfunction
    complex(dp) :: da3_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      da3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a3", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ic * k3 * a1_ef - da3_ef
    end do
  end subroutine set_magnetic_field_ef_2_component


  !> Calculates the 3-component of the perturbed magnetic field
  subroutine set_magnetic_field_ef_3_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at B1-3 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> da3 eigenfunction
    complex(dp) :: da2_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      da2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = (da2_ef - ic * k2 * a1_ef) / ef_eps
    end do
  end subroutine set_magnetic_field_ef_3_component


  !> Sets the magnetic field component parallel to B0.
  subroutine set_magnetic_field_parallel_component(derived_ef)
    !> derived eigenfunction type at B_parallel position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> B2 eigenfunction
    complex(dp) :: B2_ef(size(ef_grid))
    !> B3 eigenfunction
    complex(dp) :: B3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      B2_ef = derived_eigenfunctions( &
        findloc(derived_ef_names, "B2", dim=1) &
      )%quantities(:, i)
      B3_ef = derived_eigenfunctions( &
        findloc(derived_ef_names, "B3", dim=1) &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * B2_ef + B03_on_ef_grid * B3_ef &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_magnetic_field_parallel_component


  !> Sets the magnetic field component parallel to the B0.
  subroutine set_magnetic_field_perpendicular_component(derived_ef)
    !> derived eigenfunction type at B_perpendicular position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> B2 eigenfunction
    complex(dp) :: B2_ef(size(ef_grid))
    !> B3 eigenfunction
    complex(dp) :: B3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      B2_ef = derived_eigenfunctions( &
        findloc(derived_ef_names, "B2", dim=1) &
      )%quantities(:, i)
      B3_ef = derived_eigenfunctions( &
        findloc(derived_ef_names, "B3", dim=1) &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * B3_ef - B03_on_ef_grid * B2_ef &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_magnetic_field_perpendicular_component


  !> Calculates the divergence of the perturbed magnetic field
  subroutine set_magnetic_field_divergence(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at div(B) position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> a2 eigenfunction
    complex(dp) :: a2_ef(size(ef_grid))
    !> a3 eigenfunction
    complex(dp) :: a3_ef(size(ef_grid))
    !> da2 eigenfunction
    complex(dp) :: da2_ef(size(ef_grid))
    !> da3 eigenfunction
    complex(dp) :: da3_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      a2_ef = base_eigenfunctions(findloc(ef_names, "a2", dim=1))%quantities(:, i)
      a3_ef = base_eigenfunctions(findloc(ef_names, "a3", dim=1))%quantities(:, i)
      da2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      da3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a3", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        ic * (k2 * da3_ef - k3 * da2_ef) / ef_eps &
        + ic * k2 * (ic * k3 * a1_ef - da3_ef) / ef_eps &
        + ic * k3 * (da2_ef - ic * k2 * a1_ef) / ef_eps &
      )
    end do
  end subroutine set_magnetic_field_divergence


  !> Calculates the curl 1-component of the perturbed magnetic field
  subroutine set_magnetic_field_curl_1_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at curl(B)1 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> da2 eigenfunction
    complex(dp) :: da2_ef(size(ef_grid))
    !> da3 eigenfunction
    complex(dp) :: da3_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      da2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      da3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a3", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      derived_ef%quantities(:, i) = ( &
        ic * k2 * (da2_ef - ic * k2 * a1_ef) / ef_eps ** 2 &
        - ic * k3 * (ic * k3 * a1_ef - da3_ef) &
      )
    end do
  end subroutine set_magnetic_field_curl_1_component


  !> Calculates the curl 2-component of the perturbed magnetic field
  subroutine set_magnetic_field_curl_2_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at curl(B)2 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> da1 eigenfunction
    complex(dp) :: da1_ef(size(ef_grid))
    !> a2 eigenfunction
    complex(dp) :: a2_ef(size(ef_grid))
    !> a3 eigenfunction
    complex(dp) :: a3_ef(size(ef_grid))
    !> da2 eigenfunction
    complex(dp) :: da2_ef(size(ef_grid))
    !> dda2 eigenfunction
    complex(dp) :: dda2_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      a2_ef = base_eigenfunctions(findloc(ef_names, "a2", dim=1))%quantities(:, i)
      a3_ef = base_eigenfunctions(findloc(ef_names, "a3", dim=1))%quantities(:, i)
      da1_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a1", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      da2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      dda2_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a2", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=2 &
      )
      derived_ef%quantities(:, i) = ( &
        -k3 * (k2 * a3_ef / ef_eps - k3 * a2_ef) &
        - (dda2_ef - k2 * da1_ef) / ef_eps &
        + ef_deps * (da2_ef - ic * k2 * a1_ef) / ef_eps ** 2 &
      )
    end do
  end subroutine set_magnetic_field_curl_2_component


  !> Calculates the curl 3-component of the perturbed magnetic field
  subroutine set_magnetic_field_curl_3_component(derived_ef, right_eigenvectors)
    !> derived eigenfunction type at curl(B)3 position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> right eigenvectors
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    !> a1 eigenfunction
    complex(dp) :: a1_ef(size(ef_grid))
    !> da1 eigenfunction
    complex(dp) :: da1_ef(size(ef_grid))
    !> a2 eigenfunction
    complex(dp) :: a2_ef(size(ef_grid))
    !> a3 eigenfunction
    complex(dp) :: a3_ef(size(ef_grid))
    !> da3 eigenfunction
    complex(dp) :: da3_ef(size(ef_grid))
    !> dda3 eigenfunction
    complex(dp) :: dda3_ef(size(ef_grid))
    integer :: i, eigenvalue_idx

    do i = 1, size(ef_written_idxs)
      eigenvalue_idx = ef_written_idxs(i)
      a1_ef = base_eigenfunctions(findloc(ef_names, "a1", dim=1))%quantities(:, i)
      a2_ef = base_eigenfunctions(findloc(ef_names, "a2", dim=1))%quantities(:, i)
      a3_ef = base_eigenfunctions(findloc(ef_names, "a3", dim=1))%quantities(:, i)
      da1_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a1", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      da3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a3", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=1 &
      )
      dda3_ef = assemble_eigenfunction( &
        base_ef=base_eigenfunctions(findloc(ef_names, "a3", dim=1)), &
        eigenvector=right_eigenvectors(:, eigenvalue_idx), &
        derivative_order=2 &
      )
      derived_ef%quantities(:, i) = ( &
        k3 * da1_ef &
        - dda3_ef &
        + ef_deps * (ic * k3 * a1_ef - da3_ef) / ef_eps &
        + k2 * (k2 * a3_ef / ef_eps - k3 * a2_ef) / ef_eps &
      )
    end do
  end subroutine set_magnetic_field_curl_3_component


  !> Sets the magnetic field curl component parallel to B0.
  subroutine set_magnetic_field_parallel_curl_component(derived_ef)
    !> derived eigenfunction type at B_parallel position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> B2 eigenfunction
    complex(dp) :: curlB2(size(ef_grid))
    !> B3 eigenfunction
    complex(dp) :: curlB3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlB2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl B)2", dim=1) &
      )%quantities(:, i)
      curlB3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl B)3", dim=1) &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * curlB2 + B03_on_ef_grid * curlB3 &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_magnetic_field_parallel_curl_component


  !> Sets the magnetic field curl component perpendicular to B0.
  subroutine set_magnetic_field_perpendicular_curl_component(derived_ef)
    !> derived eigenfunction type at (curl B) perpendicular position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> B2 eigenfunction
    complex(dp) :: curlB2(size(ef_grid))
    !> B3 eigenfunction
    complex(dp) :: curlB3(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      curlB2 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl B)2", dim=1) &
      )%quantities(:, i)
      curlB3 = derived_eigenfunctions( &
        findloc(derived_ef_names, "(curl B)3", dim=1) &
      )%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * curlB3 - B03_on_ef_grid * curlB2 &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_magnetic_field_perpendicular_curl_component


  !> Sets the velocity component parallel to the magnetic field.
  subroutine set_velocity_parallel_component(derived_ef)
    !> derived eigenfunction type at v_parallel position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> v2 eigenfunction
    complex(dp) :: v2_ef(size(ef_grid))
    !> v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2", dim=1))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3", dim=1))%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * v2_ef + B03_on_ef_grid * v3_ef &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_velocity_parallel_component


  !> Sets the velocity component perpendicular to the magnetic field.
  subroutine set_velocity_perpendicular_component(derived_ef)
    !> derived eigenfunction type at v_perpendicular position in the array
    type(ef_type), intent(inout)  :: derived_ef
    !> v2 eigenfunction
    complex(dp) :: v2_ef(size(ef_grid))
    !> v3 eigenfunction
    complex(dp) :: v3_ef(size(ef_grid))
    integer :: i

    do i = 1, size(ef_written_idxs)
      v2_ef = base_eigenfunctions(findloc(ef_names, "v2", dim=1))%quantities(:, i)
      v3_ef = base_eigenfunctions(findloc(ef_names, "v3", dim=1))%quantities(:, i)
      derived_ef%quantities(:, i) = ( &
        B02_on_ef_grid * v3_ef - B03_on_ef_grid * v2_ef &
      ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    end do
  end subroutine set_velocity_perpendicular_component


  !> Cleaning routine, deallocated the derived quantities
  module procedure clean_derived_eigenfunctions
    deallocate(B02_on_ef_grid)
    deallocate(B03_on_ef_grid)
    deallocate(rho0_on_ef_grid)
    deallocate(T0_on_ef_grid)
  end procedure clean_derived_eigenfunctions

end submodule smod_derived_eigenfunctions
