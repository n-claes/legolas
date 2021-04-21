submodule (mod_matrix_manager) smod_conduction_matrix
  implicit none

contains

  module procedure add_conduction_matrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: dT0

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! temperature variables
    dT0 = T_field % d_T0_dr(gauss_idx)


    return


  end procedure add_conduction_matrix_terms

end submodule smod_conduction_matrix