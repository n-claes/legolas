submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_selfgravity
  use mod_global_variables, only: selfgravity
  implicit none

contains

  module procedure add_natural_selfgravity_Bterms
    real(dp)  :: eps

    if (.not. selfgravity) then
      return
    end if

    eps = eps_grid(grid_idx)

    ! Cubic * Cubic
    call reset_factor_positions(new_size=1)
    ! G(9, 9)
    factors(1) = eps
    positions(1, :) = [9, 9]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_cubic)
  end procedure add_natural_selfgravity_Bterms

end submodule smod_natural_bounds_selfgravity
