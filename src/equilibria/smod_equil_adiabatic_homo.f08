submodule (mod_equilibrium) smod_equil_adiabatic_homo
  implicit none
  
contains
  
  module subroutine adiabatic_homo_eq()
    
    geometry = 'Cartesian'
    call initialise_grid()
    
    k2 = dpi 
    k3 = dpi 
    
    rho_field % rho0 = 1.0d0 
    T_field % T0 = 1.0d0
    B_field % B02 = 1.0d0 
    B_field % B03 = 1.0d0 
    B_field % B0 = sqrt((B_field % B02)**2 + (B_field % B03)**2)
    
  end subroutine adiabatic_homo_eq
    
end submodule smod_equil_adiabatic_homo
    
