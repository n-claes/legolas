module mod_suite_utils
  use mod_global_variables, only: dp
  implicit none

  real(dp), parameter :: TOL = 1.0d-12

contains

  subroutine set_name(name)
    character(len=*), intent(in)  :: name

    write(*, *)
    write(*, "('" // " [Test]: " // name // "')")
  end subroutine set_name

  subroutine reset_globals()
    use mod_global_variables, only: initialise_globals, logging_level
    use mod_equilibrium_params, only: init_equilibrium_params, k2, k3

    call initialise_globals()
    call init_equilibrium_params()
    logging_level = 1 ! also print warnings
    k2 = 1.0d0
    k3 = 2.5d0
  end subroutine reset_globals


  subroutine reset_fields(init_fields)
    use mod_equilibrium, only: rho_field, equilibrium_clean, initialise_equilibrium

    logical, intent(in) :: init_fields

    if (allocated(rho_field % rho0)) then
      call equilibrium_clean()
    end if
    if (init_fields) then
      call initialise_equilibrium()
    end if
  end subroutine reset_fields


  subroutine clean_up()
    use mod_global_variables, only: radiative_cooling
    use mod_grid, only: grid, grid_clean
    use mod_radiative_cooling, only: radiative_cooling_clean

    if (allocated(grid)) then
      call grid_clean()
    end if
    call reset_fields(init_fields=.false.)
    if (radiative_cooling) then
      call radiative_cooling_clean()
    end if
  end subroutine clean_up


  subroutine create_test_grid(pts, geom, start, end)
    use mod_global_variables, only: x_start, x_end, geometry, set_gridpts
    use mod_grid, only: initialise_grid

    integer, intent(in)             :: pts
    character(len=*), intent(in)    :: geom
    real(dp), intent(in), optional  :: start, end

    geometry = geom
    if (present(start)) then
      x_start = start
    else
      x_start = 0.0d0
    end if
    if (present(end)) then
      x_end = end
    else
      x_end = 1.0d0
    end if
    call set_gridpts(pts)
    call initialise_grid()
  end subroutine create_test_grid


  subroutine set_default_units()
    use mod_global_variables, only: cgs_units
    use mod_units, only: set_normalisations

    cgs_units = .true.
    call set_normalisations( &
      new_unit_temperature=1.0d6, &
      new_unit_magneticfield=5.0d0, &
      new_unit_length=1.0d10, &
      new_mean_molecular_weight=1.0d0 &
    )
  end subroutine set_default_units


  function linspace(x0, x1, xvals) result(xarray)
    real(dp), intent(in)  :: x0, x1
    integer, intent(in)   :: xvals
    real(dp)  :: dx, xarray(xvals)
    integer   :: i

    dx = (x1 - x0) / (xvals - 1)
    do i = 1, xvals
      xarray(i) = x0 + (i - 1) * dx
    end do
  end function linspace


  subroutine sort_complex_array(array)
    complex(dp), intent(inout)  :: array(:)
    complex(dp) :: temp
    integer     :: i, minidx

    ! sort array using selection sort, based on real part
    do i = 1, size(array) - 1
      minidx = minloc(real(array(i:)), 1) + i - 1
      if (real(array(i)) > real(array(minidx))) then
        temp = array(i)
        array(i) = array(minidx)
        array(minidx) = temp
      end if
    end do
  end subroutine sort_complex_array


  subroutine create_banded_array(subdiags, superdiags, mat)
    integer, intent(in) :: subdiags, superdiags
    complex(dp), intent(out) :: mat(8, 8)
    integer :: i

    if (superdiags > 2 .or. subdiags > 2) then
      write(*, *) "can only create banded array up to 2 sub/super diagonals"
      stop
    end if

    mat = (0.0d0, 0.0d0)
    ! diagonal
    do i = 1, 8
      mat(i, i) = (1.0d0, 2.0d0) * i
    end do
    ! diagonals 1
    do i = 1, 7
      if (superdiags >= 1) mat(i, i + 1) = -2.0d0 * i + (3.0d0, 5.0d0)
      if (subdiags >= 1) mat(i + 1, i) =  cmplx(1.5d0, 0.5d0, kind=dp) * i
    end do
    ! dieagonals 2
    do i = 1, 6
      if (superdiags >= 2) mat(i, i + 2) = cmplx(8.0d0 - 2.0d0 * i, 1.5d0 * i, kind=dp)
      if (subdiags >= 2) mat(i + 2, i) = cmplx(1.0d0, -2.5d0 * i, kind=dp)
    end do
  end subroutine create_banded_array


end module mod_suite_utils
