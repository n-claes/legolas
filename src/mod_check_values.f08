module mod_check_values
  use mod_global_variables, only: dp, dp_LIMIT
  implicit none

  private

  interface check_small_values
    module procedure small_values_real_array
    module procedure small_values_complex_array
    module procedure small_values_real_matrix
    module procedure small_values_complex_matrix
  end interface check_small_values

  public :: check_small_values
  public :: check_negative_array
  public :: check_equilibrium_conditions

contains

  subroutine small_values_real_array(array)
    real(dp), intent(inout) :: array(:)
    integer                 :: i

    do i = 1, size(array)
      if (abs(array(i)) < dp_LIMIT) then
        array(i) = 0.0d0
      end if
    end do
  end subroutine small_values_real_array

  subroutine small_values_complex_array(array)
    complex(dp), intent(inout) :: array(:)
    integer                    :: i
    real(dp)                   :: a_real, a_imag

    do i = 1, size(array)
      a_real = real(array(i))
      a_imag = aimag(array(i))

      if (abs(a_real) < dp_LIMIT) then
        a_real = 0.0d0
      end if
      if (abs(a_imag) < dp_LIMIT) then
        a_imag = 0.0d0
      end if

      array(i) = cmplx(a_real, a_imag, kind=dp)
    end do
  end subroutine small_values_complex_array

  subroutine small_values_real_matrix(matrix)
    real(dp), intent(inout)    :: matrix(:, :)
    integer                    :: i, j

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        if (abs(matrix(i, j)) < dp_LIMIT) then
          matrix(i, j) = 0.0d0
        end if
      end do
    end do
  end subroutine small_values_real_matrix

  subroutine small_values_complex_matrix(matrix)
    complex(dp), intent(inout) :: matrix(:, :)
    integer                    :: i, j
    real(dp)                   :: a_real, a_imag

    do j = 1, size(matrix(1, :))
      do i = 1, size(matrix(:, 1))
        a_real = real(matrix(i, j))
        a_imag = aimag(matrix(i, j))

        if (abs(a_real) < dp_LIMIT) then
          a_real = 0.0d0
        end if
        if (abs(a_imag) < dp_LIMIT) then
          a_imag = 0.0d0
        end if

        matrix(i, j) = cmplx(a_real, a_imag, kind=dp)
      end do
    end do
  end subroutine small_values_complex_matrix

  subroutine check_negative_array(array, variable_name)
    real(dp), intent(in)          :: array(:)
    character(len=*), intent(in)  :: variable_name
    integer                       :: i

    do i = 1, size(array)
      if (array(i) < 0.0d0) then
        write(*, *) "WARNING: ", trim(variable_name), " is negative somewhere!"
        stop
      end if
    end do
  end subroutine check_negative_array
  
  
  subroutine check_equilibrium_conditions(rho_field, T_field, B_field, v_field, grav_field)
    use mod_types, only: density_type, temperature_type, bfield_type, velocity_type, gravity_type
    use mod_global_variables, only: geometry, dp_LIMIT, gauss_gridpts
    use mod_grid, only: eps_grid, d_eps_grid_dr
    
    type(density_type), intent(in)      :: rho_field
    type(temperature_type), intent(in)  :: T_field
    type(bfield_type), intent(in)       :: B_field 
    type(velocity_type), intent(in)     :: v_field 
    type(gravity_type), intent(in)      :: grav_field
    
    real(dp)    :: rho, drho, B02, dB02, B03, dB03, T0, dT0, grav, v02, v03
    real(dp)    :: eps, d_eps, eq_limit
    real(dp)    :: eq_cond(gauss_gridpts)
    integer     :: i 
    
    do i = 1, gauss_gridpts
      rho = rho_field % rho0(i)
      drho = rho_field % d_rho0_dr(i)
      B02 = B_field % B02(i)
      B03 = B_field % B03(i) 
      dB02 = B_field % d_B02_dr(i)
      dB03 = B_field % d_B03_dr(i)
      T0 = T_field % T0(i) 
      dT0 = T_field % d_T0_dr(i)
      grav = grav_field % grav(i)
      v02 = v_field % v02(i)
      v03 = v_field % v03(i)
      eps = eps_grid(i)
      d_eps = d_eps_grid_dr(i)
    
      eq_cond(i) = drho * T0 + rho * dT0 + B02 * dB02 + B03 * dB03 + rho * grav &
                   - (d_eps/eps) * (rho * v02**2 - B02**2)

      if (eq_cond(i) > dp_LIMIT) then 
        write(*, *) "WARNING: equilibrium conditions not met!"
        stop 
      end if 
    end do 
    
    eq_limit = 1.0d0
    
    if (geometry == 'cylindrical') then 
      if (abs(B_field % B02(1)) > eq_limit) then
        write(*, *) "WARNING: B_theta(0) is non-zero!"
        stop
      else if (abs(B_field % d_B03_dr(i)) > eq_limit) then
        write(*, *) "WARNING: dB_z/dr(0) is non-zero!"
        stop
      else if (abs(v_field % v02(i)) > eq_limit) then
        write(*, *) "WARNING: v_theta(0) is non-zero!"
        stop
      else if (abs(v_field % d_v03_dr(1)) > eq_limit) then
        write(*, *) "WARNING: d_v03_dr(0) is non-zero!"
        stop
      end if
    end if
    
  end subroutine check_equilibrium_conditions
    


end module mod_check_values
