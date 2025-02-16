module mod_iv_initial_conditions
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_logging, only: logger
  use mod_iv_globals, only: iv_prof_fcn_ptr_t, profile_fcn
  implicit none

contains

  subroutine get_f_lists(f_list, df_list)
    type(iv_prof_fcn_ptr_t), intent(out) :: f_list(:), df_list(:)
    integer :: num_profiles

    ! FIXME: set this automatically
    num_profiles = 2  ! Change this manually for now
    if ((size(f_list) /= num_profiles) .or. size(df_list) /= num_profiles) then
      call logger%error("Number of function profiles does not match provided list size")
    end if

    ! Add the profiles here
    f_list(1)%ptr => rho
    f_list(2)%ptr => zeros
    ! f_list(3)%ptr => zeros

    df_list(1)%ptr => drho
    df_list(2)%ptr => zeros
    ! df_list(3)%ptr => zeros
  end subroutine get_f_lists

  ! -----------------------------------------------------------------
  ! Initial profiles
  ! -----------------------------------------------------------------
  function rho(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))

    ! res = 0.0d0
    ! where (x > 0.4 .and. x < 0.6)
    !   res = sin(5 * dpi * (x - 0.4))
    ! end where

    res = 0.1 * gaussian(x, 0.5d0, 0.05d0)
  end function rho

  function drho(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))

    ! res = 0.0d0
    ! where (x > 0.4 .and. x < 0.6)
    !   res = 5 * dpi * cos(5 * dpi * (x - 0.4))
    ! end where
    
    res = 0.1 * dgaussian(x, 0.5d0, 0.05d0)
  end function drho


  function T(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 0.0d0
  end function T

  function dT(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 0.0d0
  end function dT


  function v1(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 0.1* gaussian(x, 0.5d0, 0.05d0)
  end function v1

  function dv1(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 0.1 * dgaussian(x, 0.5d0, 0.05d0)
  end function dv1

  function zeros(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 0.0d0
  end function zeros

  function ones(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = 1.0d0
  end function ones


  real(dp) elemental function gaussian(x, mean, std_dev)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: mean
    real(dp), intent(in) :: std_dev

    gaussian = (1.0d0 / (std_dev * sqrt(2.0d0*dpi))) * exp(-0.5d0 * ((x - mean)/std_dev)**2)

  end function gaussian

  real(dp) elemental function dgaussian(x, mean, std_dev)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: mean
    real(dp), intent(in) :: std_dev

    dgaussian = -((x - mean) / (std_dev**2)) * gaussian(x, mean, std_dev)

  end function dgaussian

end module mod_iv_initial_conditions