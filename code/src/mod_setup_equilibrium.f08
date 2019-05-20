module mod_setup_equilibrium
  implicit none

  real, allocatable         ::  rho_0(:)
  real, allocatable         ::  v_0(:, :, :)
  real, allocatable         ::  T_0(:)
  real, allocatable         ::  B_0(:, :, :)


contains

  subroutine initialise_equilibrium()
    use mod_global_variables

    integer                 ::  i

    allocate(rho_0(gridpts))
    allocate(v_0(gridpts, gridpts, gridpts))
    allocate(T_0(gridpts))
    allocate(B_0(gridpts, gridpts, gridpts))

    ! Temporary initialisations
    do i = 1, gridpts
      rho_0(i)     = 1.0d0
      v_0(i, i, i) = 1.0d0
      T_0(i)       = 1.0d0
      B_0(i, i, i) = 1.0d0
    end do

  end subroutine initialise_equilibrium


end module mod_setup_equilibrium
