! =============================================================================
!> This submodule calculates quantities derived from the base eigenfunctions that are
!! parallel or perpendicular to the background magnetic field:
!!  - parallel and perpendicular components of \(\mathbf{B}_1\)
!!  - parallel and perpendicular components of \(\nabla \times \mathbf{B}_1\)
!!  - parallel and perpendicular components of \(\mathbf{v}_1\)
!!  - parallel and perpendicular components of \(\nabla \times \mathbf{v}_1\)
submodule(mod_eigenfunctions:smod_derived_efs) smod_derived_pp_efs
  implicit none

contains

  !> Sets the various parallel/perpendicular quantities.
  module procedure set_pp_quantities
    call set_magnetic_field_pp(loc=locB, ef_index=ef_index, state_vector=state_vector)
    call set_magnetic_field_curl_pp( &
      loc=locCurlB, ef_index=ef_index, state_vector=state_vector &
    )
    call set_velocity_pp(loc=locV, ef_index=ef_index)
    call set_velocity_curl_pp( &
      loc=locCurlV, ef_index=ef_index, state_vector=state_vector &
    )
  end procedure set_pp_quantities


  !> Sets the parallel and perpendicular components of the perturbed magnetic field
  !! with respect to the background magnetic field.
  subroutine set_magnetic_field_pp(loc, ef_index, state_vector)
    !> position indices in the main array
    integer, intent(in) :: loc(2)
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index
    !> state vector
    character(len=*), intent(in) :: state_vector(:)
    !> perturbed B2 eigenfunction
    complex(dp) :: B2_ef(size(ef_grid))
    !> perturbed B3 eigenfunction
    complex(dp) :: B3_ef(size(ef_grid))

    B2_ef = retrieve_eigenfunction_from_index("B2", state_vector, ef_index=ef_index)
    B3_ef = retrieve_eigenfunction_from_index("B3", state_vector, ef_index=ef_index)
    ! parallel component
    derived_eigenfunctions(loc(1))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * B2_ef + B03_on_ef_grid * B3_ef &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
    ! perpendicular component
    derived_eigenfunctions(loc(2))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * B3_ef - B03_on_ef_grid * B2_ef &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
  end subroutine set_magnetic_field_pp


  !> Sets the parallel and perpendicular components of the perturbed
  !! magnetic field curl with respect to the background magnetic field.
  subroutine set_magnetic_field_curl_pp(loc, ef_index, state_vector)
    !> position indices in the main array
    integer, intent(in) :: loc(2)
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index
    !> state vector
    character(len=*), intent(in) :: state_vector(:)
    !> curl of B2
    complex(dp):: curlB2(size(ef_grid))
    !> curl of B3
    complex(dp):: curlB3(size(ef_grid))

    curlB2 = retrieve_eigenfunction_from_index( &
      "(curl B)2", state_vector, ef_index=ef_index &
    )
    curlB3 = retrieve_eigenfunction_from_index( &
      "(curl B)3", state_vector, ef_index=ef_index &
    )
    ! parallel component
    derived_eigenfunctions(loc(1))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * curlB2 + B03_on_ef_grid * curlB3 &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
    ! perpendicular component
    derived_eigenfunctions(loc(2))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * curlB3 - B03_on_ef_grid * curlB2 &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
  end subroutine set_magnetic_field_curl_pp


  !> Sets the parallel and perpendicular components of the perturbed velocity with
  !! respect to the background magnetic field.
  subroutine set_velocity_pp(loc, ef_index)
    !> position indices in the main array
    integer, intent(in) :: loc(2)
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index

    ! parallel component
    derived_eigenfunctions(loc(1))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * v2_ef + B03_on_ef_grid * v3_ef &
    ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
    ! perpendicular component
    derived_eigenfunctions(loc(2))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * v3_ef - B03_on_ef_grid * v2_ef &
    ) / sqrt(B02_on_ef_grid**2 + B03_on_ef_grid**2)
  end subroutine set_velocity_pp


  !> Sets the parallel and perpendicular components of the velocity curl with
  !! respect to the background magnetic field.
  subroutine set_velocity_curl_pp(loc, ef_index, state_vector)
    !> position indices in the main array
    integer, intent(in) :: loc(2)
    !> index of the eigenfunction in the "quantities" array attribute
    integer, intent(in) :: ef_index
    !> state vector
    character(len=*), intent(in) :: state_vector(:)
    !> curl of v2
    complex(dp):: curlv2(size(ef_grid))
    !> curl of v3
    complex(dp):: curlv3(size(ef_grid))

    curlv2 = retrieve_eigenfunction_from_index( &
      "(curl v)2", state_vector, ef_index=ef_index &
    )
  curlv3 = retrieve_eigenfunction_from_index( &
      "(curl v)3", state_vector, ef_index=ef_index &
    )
    ! parallel component
    derived_eigenfunctions(loc(1))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * curlv2 + B03_on_ef_grid * curlv3 &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
    ! perpendicular component
    derived_eigenfunctions(loc(2))%quantities(:, ef_index) = ( &
      B02_on_ef_grid * curlv3 - B03_on_ef_grid * curlv2 &
    ) / sqrt(B02_on_ef_grid ** 2 + B03_on_ef_grid ** 2)
  end subroutine set_velocity_curl_pp

end submodule smod_derived_pp_efs
