module mod_derived_efs
  use mod_global_variables, only: dp, str_len_arr, ic
  use mod_equilibrium_params, only: k2, k3
  use mod_settings, only: settings_t
  use mod_background, only: background_t
  use mod_ef_assembly, only: assemble_eigenfunction, retransform_eigenfunction, &
    get_ef_eps, get_ef_deps
  use mod_derived_ef_names
  use mod_logging, only: logger
  use mod_get_indices, only: get_index
  implicit none

  private

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

  type, public :: derived_ef_t
    character(str_len_arr) :: name
    complex(dp), allocatable :: quantities(:, :)
    procedure(derived_ef_func), pointer, private, nopass :: get_derived_ef

  contains

    procedure, public :: initialise
    procedure, public :: assemble
    procedure, public :: delete

    procedure, private :: set_function_pointer

  end type derived_ef_t

  real(dp), allocatable :: rho0_on_ef_grid(:)
  real(dp), allocatable :: T0_on_ef_grid(:)
  real(dp), allocatable :: B02_on_ef_grid(:)
  real(dp), allocatable :: B03_on_ef_grid(:)

  public :: deallocate_derived_ef_module_variables

contains

  subroutine initialise(this, name, ef_grid_size, nb_efs)
    class(derived_ef_t), intent(inout) :: this
    character(str_len_arr), intent(in) :: name
    integer, intent(in) :: ef_grid_size
    integer, intent(in) :: nb_efs

    this%name = name
    allocate(this%quantities(ef_grid_size, nb_efs))
    call this%set_function_pointer()
  end subroutine initialise


  subroutine assemble( &
    this, settings, background, idxs_to_assemble, right_eigenvectors, ef_grid &
  )
    class(derived_ef_t), intent(inout) :: this
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    integer, intent(in) :: idxs_to_assemble(:)
    complex(dp), intent(in) :: right_eigenvectors(:, :)
    real(dp), intent(in) :: ef_grid(:)
    integer :: i, idx

    call set_equilibrium_arrays_on_ef_grid(background, ef_grid)

    do i = 1, size(idxs_to_assemble)
      idx = idxs_to_assemble(i)
      this%quantities(:, i) = this%get_derived_ef( &
        settings=settings, &
        eigenvector=right_eigenvectors(:, idx), &
        ef_grid=ef_grid &
      )
    end do
  end subroutine assemble


  pure subroutine delete(this)
    class(derived_ef_t), intent(inout) :: this
    if (allocated(this%quantities)) deallocate(this%quantities)
    nullify(this%get_derived_ef)
  end subroutine delete


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


  subroutine set_function_pointer(this)
    class(derived_ef_t), intent(inout) :: this

    select case(this%name)
    case(S_name)
      this%get_derived_ef => get_entropy
    case(div_v_name)
      this%get_derived_ef => get_div_v
    case(curl_v_1_name)
      this%get_derived_ef => get_curl_v_1
    case(curl_v_2_name)
      this%get_derived_ef => get_curl_v_2
    case(curl_v_3_name)
      this%get_derived_ef => get_curl_v_3
    case(B1_name)
      this%get_derived_ef => get_B1
    case(B2_name)
      this%get_derived_ef => get_B2
    case(B3_name)
      this%get_derived_ef => get_B3
    case(div_B_name)
      this%get_derived_ef => get_div_B
    case(curl_B_1_name)
      this%get_derived_ef => get_curl_B_1
    case(curl_B_2_name)
      this%get_derived_ef => get_curl_B_2
    case(curl_B_3_name)
      this%get_derived_ef => get_curl_B_3
    case(B_para_name)
      this%get_derived_ef => get_B_para
    case(B_perp_name)
      this%get_derived_ef => get_B_perp
    case(curl_B_para_name)
      this%get_derived_ef => get_curl_B_para
    case(curl_B_perp_name)
      this%get_derived_ef => get_curl_B_perp
    case(v_para_name)
      this%get_derived_ef => get_v_para
    case(v_perp_name)
      this%get_derived_ef => get_v_perp
    case(curl_v_para_name)
      this%get_derived_ef => get_curl_v_para
    case(curl_v_perp_name)
      this%get_derived_ef => get_curl_v_perp
    case default
      call logger%error( &
        "derived ef assembly -- unknown eigenfunction name: "// trim(this%name) &
      )
      nullify(this%get_derived_ef)
      return
    end select
  end subroutine set_function_pointer


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


  function get_B1(settings, eigenvector, ef_grid) result(B1)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: B1(size(ef_grid))
    complex(dp) :: a2(size(ef_grid)), a3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    a2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector)
    a3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector)
    ef_eps = get_ef_eps(settings, ef_grid)
    B1 = ic * (k2 * a3 / ef_eps - k3 * a2)
  end function get_B1


  function get_B2(settings, eigenvector, ef_grid) result(B2)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: B2(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da3(size(ef_grid))

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector, diff_order=1)
    B2 = ic * k3 * a1 - da3
  end function get_B2


  function get_B3(settings, eigenvector, ef_grid) result(B3)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: B3(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da2(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector, diff_order=1)
    ef_eps = get_ef_eps(settings, ef_grid)
    B3 = (da2 - ic * k2 * a1) / ef_eps
  end function get_B3


  function get_div_B(settings, eigenvector, ef_grid) result(div_B)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: div_B(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da2(size(ef_grid)), da3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector, diff_order=1)
    da3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector, diff_order=1)
    ef_eps = get_ef_eps(settings, ef_grid)
    div_B = ( &
      ic * (k2 * da3 - k3 * da2) &
      + ic * k2 * (ic * k3 * a1 - da3) &
      + ic * k3 * (da2 - ic * k2 * a1) &
    ) / ef_eps
  end function get_div_B


  function get_curl_B_1(settings, eigenvector, ef_grid) result(curl_B_1)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_B_1(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da2(size(ef_grid)), da3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid))

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector, diff_order=1)
    da3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector, diff_order=1)
    ef_eps = get_ef_eps(settings, ef_grid)
    curl_B_1 = ( &
      ic * k2 * (da2 - ic * k2 * a1) / ef_eps**2 - ic * k3 * (ic * k3 * a1 - da3) &
    )
  end function get_curl_B_1


  function get_curl_B_2(settings, eigenvector, ef_grid) result(curl_B_2)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_B_2(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da1(size(ef_grid))
    complex(dp) :: a2(size(ef_grid)), da2(size(ef_grid)), dda2(size(ef_grid))
    complex(dp) :: a3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid)), ef_deps

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector, diff_order=1)
    a2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector)
    da2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector, diff_order=1)
    dda2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector, diff_order=2)
    a3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector)
    ef_eps = get_ef_eps(settings, ef_grid)
    ef_deps = get_ef_deps(settings)
    curl_B_2 = ( &
      -k3 * (k2 * a3 / ef_eps - k3 * a2) &
      - (dda2 - k2 * da1) / ef_eps &
      + ef_deps * (da2 - ic * k2 * a1) / ef_eps ** 2 &
    )
  end function get_curl_B_2


  function get_curl_B_3(settings, eigenvector, ef_grid) result(curl_B_3)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_B_3(size(ef_grid))
    complex(dp) :: a1(size(ef_grid)), da1(size(ef_grid))
    complex(dp) :: a2(size(ef_grid))
    complex(dp) :: a3(size(ef_grid)), da3(size(ef_grid)), dda3(size(ef_grid))
    real(dp) :: ef_eps(size(ef_grid)), ef_deps

    a1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector)
    da1 = get_base_eigenfunction("a1", settings, ef_grid, eigenvector, diff_order=1)
    a2 = get_base_eigenfunction("a2", settings, ef_grid, eigenvector)
    a3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector)
    da3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector, diff_order=1)
    dda3 = get_base_eigenfunction("a3", settings, ef_grid, eigenvector, diff_order=2)
    ef_eps = get_ef_eps(settings, ef_grid)
    ef_deps = get_ef_deps(settings)
    curl_B_3 = ( &
      k3 * da1 &
      - dda3 &
      + ef_deps * (ic * k3 * a1 - da3) / ef_eps &
      + k2 * (k2 * a3 / ef_eps - k3 * a2) / ef_eps &
    )
  end function get_curl_B_3


  function get_B_para(settings, eigenvector, ef_grid) result(B_para)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: B_para(size(ef_grid))
    complex(dp) :: B2(size(ef_grid)), B3(size(ef_grid))

    B2 = get_B2(settings, eigenvector, ef_grid)
    B3 = get_B3(settings, eigenvector, ef_grid)
    B_para = ( &
      (B02_on_ef_grid * B2 + B03_on_ef_grid * B3) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_B_para


  function get_B_perp(settings, eigenvector, ef_grid) result(B_perp)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: B_perp(size(ef_grid))
    complex(dp) :: B2(size(ef_grid)), B3(size(ef_grid))

    B2 = get_B2(settings, eigenvector, ef_grid)
    B3 = get_B3(settings, eigenvector, ef_grid)
    B_perp = ( &
      (B02_on_ef_grid * B3 - B03_on_ef_grid * B2) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_B_perp


  function get_curl_B_para(settings, eigenvector, ef_grid) result(curl_B_para)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_B_para(size(ef_grid))
    complex(dp) :: curl_B_2(size(ef_grid)), curl_B_3(size(ef_grid))

    curl_B_2 = get_curl_B_2(settings, eigenvector, ef_grid)
    curl_B_3 = get_curl_B_3(settings, eigenvector, ef_grid)
    curl_B_para = ( &
      (B02_on_ef_grid * curl_B_2 + B03_on_ef_grid * curl_B_3) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_curl_B_para


  function get_curl_B_perp(settings, eigenvector, ef_grid) result(curl_B_perp)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_B_perp(size(ef_grid))
    complex(dp) :: curl_B_2(size(ef_grid)), curl_B_3(size(ef_grid))

    curl_B_2 = get_curl_B_2(settings, eigenvector, ef_grid)
    curl_B_3 = get_curl_B_3(settings, eigenvector, ef_grid)
    curl_B_perp = ( &
      (B02_on_ef_grid * curl_B_3 - B03_on_ef_grid * curl_B_2) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_curl_B_perp


  function get_v_para(settings, eigenvector, ef_grid) result(v_para)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: v_para(size(ef_grid))
    complex(dp) :: v2(size(ef_grid)), v3(size(ef_grid))

    v2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector)
    v3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector)
    v_para = ( &
      (B02_on_ef_grid * v2 + B03_on_ef_grid * v3) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_v_para


  function get_v_perp(settings, eigenvector, ef_grid) result(v_perp)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: v_perp(size(ef_grid))
    complex(dp) :: v2(size(ef_grid)), v3(size(ef_grid))

    v2 = get_base_eigenfunction("v2", settings, ef_grid, eigenvector)
    v3 = get_base_eigenfunction("v3", settings, ef_grid, eigenvector)
    v_perp = ( &
      (B02_on_ef_grid * v3 - B03_on_ef_grid * v2) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_v_perp


  function get_curl_v_para(settings, eigenvector, ef_grid) result(curl_v_para)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_v_para(size(ef_grid))
    complex(dp) :: curl_v_2(size(ef_grid)), curl_v_3(size(ef_grid))

    curl_v_2 = get_curl_v_2(settings, eigenvector, ef_grid)
    curl_v_3 = get_curl_v_3(settings, eigenvector, ef_grid)
    curl_v_para = ( &
      (B02_on_ef_grid * curl_v_2 + B03_on_ef_grid * curl_v_3) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_curl_v_para


  function get_curl_v_perp(settings, eigenvector, ef_grid) result(curl_v_perp)
    type(settings_t), intent(in) :: settings
    complex(dp), intent(in) :: eigenvector(:)
    real(dp), intent(in) :: ef_grid(:)
    complex(dp) :: curl_v_perp(size(ef_grid))
    complex(dp) :: curl_v_2(size(ef_grid)), curl_v_3(size(ef_grid))

    curl_v_2 = get_curl_v_2(settings, eigenvector, ef_grid)
    curl_v_3 = get_curl_v_3(settings, eigenvector, ef_grid)
    curl_v_perp = ( &
      (B02_on_ef_grid * curl_v_3 - B03_on_ef_grid * curl_v_2) &
      / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2) &
    )
  end function get_curl_v_perp


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


  subroutine set_equilibrium_arrays_on_ef_grid(background, ef_grid)
    use mod_grid, only: grid_gauss
    use mod_function_utils, only: from_function

    type(background_t), intent(in) :: background
    real(dp), intent(in) :: ef_grid(:)

    if (.not. allocated(rho0_on_ef_grid)) then
      rho0_on_ef_grid = interpolate_array_on_ef_grid( &
        from_function(background%density%rho0, grid_gauss), ef_grid &
      )
      T0_on_ef_grid = interpolate_array_on_ef_grid( &
        from_function(background%temperature%T0, grid_gauss), ef_grid &
      )
      B02_on_ef_grid = interpolate_array_on_ef_grid( &
        from_function(background%magnetic%B02, grid_gauss), ef_grid &
      )
      B03_on_ef_grid = interpolate_array_on_ef_grid( &
        from_function(background%magnetic%B03, grid_gauss), ef_grid &
      )
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
