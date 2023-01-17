module mod_derived_efs
  use mod_global_variables, only: dp, str_len_arr, ic
  use mod_equilibrium_params, only: k2, k3
  use mod_settings, only: settings_t
  use mod_base_efs, only: base_ef_t
  use mod_ef_assembly, only: assemble_eigenfunction, retransform_eigenfunction, &
    get_ef_eps, get_ef_deps
  use mod_derived_ef_names
  use mod_logging, only: logger
  use mod_get_indices, only: get_index
  implicit none

  private

  type, extends(base_ef_t), public :: derived_ef_t
  contains

    procedure,  public :: assemble

  end type derived_ef_t

  interface
    function derived_ef_func(settings, eigenvector, ef_grid) result(ef)
      use mod_global_variables, only: dp
      use mod_settings, only: settings_t
      type(settings_t), intent(in) :: settings
      complex(dp), intent(in) :: eigenvector(:)
      real(dp), intent(in) :: ef_grid(:)
      complex(dp) :: ef(size(ef_grid))
    end function derived_ef_func
  end interface

  real(dp), allocatable :: rho0_on_ef_grid(:)
  real(dp), allocatable :: T0_on_ef_grid(:)
  real(dp), allocatable :: B02_on_ef_grid(:)
  real(dp), allocatable :: B03_on_ef_grid(:)

  public :: deallocate_derived_ef_module_variables

contains

  subroutine assemble(this, settings, idxs_to_assemble, right_eigenvectors, ef_grid)
    class(derived_ef_t), intent(inout) :: this
    type(settings_t), intent(in) :: settings
    integer, intent(in) :: idxs_to_assemble(:)
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    real(dp), intent(in) :: ef_grid(:)
    integer :: i, idx
    procedure(derived_ef_func), pointer :: get_derived_ef

    call set_equilibrium_arrays_on_ef_grid(ef_grid)

    do i = 1, size(idxs_to_assemble)
      idx = idxs_to_assemble(i)
      select case(this%name)
      case(S_name)
        get_derived_ef => get_entropy
      case(div_v_name)
        get_derived_ef => get_div_v
      case(curl_v_1_name)
        get_derived_ef => get_curl_v_1
      case(curl_v_2_name)
        get_derived_ef => get_curl_v_2
      case(curl_v_3_name)
        get_derived_ef => get_curl_v_3
      case default
        call logger%error( &
          "derived ef assembly -- unknown eigenfunction name: "// trim(this%name) &
        )
        nullify(get_derived_ef)
        return
      end select

      this%quantities(:, i) = get_derived_ef( &
        settings=settings, &
        eigenvector=right_eigenvectors(:, idx), &
        ef_grid=ef_grid &
      )
      nullify(get_derived_ef)
    end do
  end subroutine assemble


  function get_entropy(settings, eigenvector, ef_grid) result(entropy)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: entropy(size(ef_grid))
    complex(dp) :: rho(size(ef_grid)), T(size(ef_grid))

    rho = get_base_eigenfunction("rho", settings, ef_grid, eigenvector)
    T = get_base_eigenfunction("T", settings, ef_grid, eigenvector)
    entropy = ( &
      T / rho0_on_ef_grid ** (2.0_dp / 3.0_dp) &
      - (2.0_dp / 3.0_dp) * ( &
        rho * T0_on_ef_grid / rho0_on_ef_grid ** (5.0_dp / 3.0_dp) &
      ) &
    )
  end function get_entropy


  function get_div_v(settings, eigenvector, ef_grid) result(div_v)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: div_v(size(ef_grid))
    complex(dp) :: dv1(size(ef_grid)), v2(size(ef_grid)), v3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    dv1 = get_base_eigenfunction("v1", settings, ef_grid, eigenvector, diff_order=1)
    v2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector)
    v3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector)
    ef_eps = get_ef_eps(settings, ef_grid)
    div_v = ic * (-dv1 / ef_eps + k2 * v2 / ef_eps + k3 * v3)
  end function get_div_v


  function get_curl_v_1(settings, eigenvector, ef_grid) result(curl_v_1)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_v_1(size(ef_grid))
    complex(dp) :: v2(size(ef_grid)), v3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    v2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector)
    v3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector)
    ef_eps = get_ef_eps(settings, ef_grid)
    curl_v_1 = ic * (k2 * v3 / ef_eps - k3 * v2)
  end function get_curl_v_1


  function get_curl_v_2(settings, eigenvector, ef_grid) result(curl_v_2)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_v_2(size(ef_grid))
    complex(dp) :: v1(size(ef_grid)), v3(size(ef_grid)), dv3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid)), ef_deps

    v1 = get_base_eigenfunction("v1", settings, ef_grid, eigenvector)
    v3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector)
    dv3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector, diff_order=1)
    ef_eps = get_ef_eps(settings, ef_grid)
    ef_deps = get_ef_deps(settings)
    curl_v_2 = ic * k3 * v1 - dv3 / ef_eps - ef_deps * v3 / ef_eps
  end function get_curl_v_2


  function get_curl_v_3(settings, eigenvector, ef_grid) result(curl_v_3)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_v_3(size(ef_grid))
    complex(dp) :: v1(size(ef_grid)), v2(size(ef_grid)), dv2(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid)), ef_deps

    v1 = get_base_eigenfunction("v1", settings, ef_grid, eigenvector)
    v2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector)
    dv2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector, diff_order=1)
    ef_eps = get_ef_eps(settings, ef_grid)
    ef_deps = get_ef_deps(settings)
    curl_v_3 = dv2 + (ef_deps * v2 - ic * k2 * v1) / ef_eps
  end function get_curl_v_3


  function get_base_eigenfunction( &
    name, settings, ef_grid, eigenvector, diff_order &
  ) result(base_ef)
    character(len=*), intent(in) :: name
    type(settings_t), intent(in) :: settings
    real(dp), intent(in) :: ef_grid(:)
    complex(dp), intent(in) :: eigenvector(:)
    integer, intent(in), optional :: diff_order

    complex(dp) :: base_ef(size(ef_grid))
    integer :: derivative_order

    derivative_order = 0
    if (present(diff_order)) derivative_order = diff_order

    base_ef = assemble_eigenfunction( &
      settings=settings, &
      ef_name=name, &
      ef_grid=ef_grid, &
      state_vector_index=get_index(name, settings%get_state_vector()), &
      eigenvector=eigenvector, &
      derivative_order=derivative_order &
    )

    if (derivative_order /= 0) return
    ! only retransform for "true" base eigenfunctions, not their derivatives
    base_ef = retransform_eigenfunction( &
      ef_name=name, ef_eps=get_ef_eps(settings, ef_grid), eigenfunction=base_ef &
    )
  end function get_base_eigenfunction


  subroutine set_equilibrium_arrays_on_ef_grid(ef_grid)
    use mod_equilibrium, only: rho_field, B_field, T_field

    real(dp), intent(in) :: ef_grid(:)

    if (.not. allocated(rho0_on_ef_grid)) then
      rho0_on_ef_grid = interpolate_array_on_ef_grid(rho_field%rho0, ef_grid)
      T0_on_ef_grid = interpolate_array_on_ef_grid(T_field%T0, ef_grid)
      B02_on_ef_grid = interpolate_array_on_ef_grid(B_field%B02, ef_grid)
      B03_on_ef_grid = interpolate_array_on_ef_grid(B_field%B03, ef_grid)
    end if
  end subroutine set_equilibrium_arrays_on_ef_grid


  function interpolate_array_on_ef_grid(array, ef_grid) result(array_on_ef_grid)
    use mod_interpolation, only: lookup_table_value
    use mod_grid, only: grid_gauss

    real(dp), intent(in) :: array(:)
    real(dp), intent(in) :: ef_grid(:)
    real(dp), allocatable :: array_on_ef_grid(:)
    integer :: i

    allocate(array_on_ef_grid(size(ef_grid)))
    do i = 1, size(ef_grid)
      array_on_ef_grid(i) = lookup_table_value( &
        ef_grid(i), grid_gauss, array, allow_outside=.true. &
      )
    end do
  end function interpolate_array_on_ef_grid


  subroutine deallocate_derived_ef_module_variables()
    if (allocated(rho0_on_ef_grid)) deallocate(rho0_on_ef_grid)
    if (allocated(T0_on_ef_grid)) deallocate(T0_on_ef_grid)
    if (allocated(B02_on_ef_grid)) deallocate(B02_on_ef_grid)
    if (allocated(B03_on_ef_grid)) deallocate(B03_on_ef_grid)
  end subroutine deallocate_derived_ef_module_variables


end module mod_derived_efs
