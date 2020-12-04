set (sources
        mod_version.f08
        mod_global_variables.f08
        dataIO/mod_painting.f08
        mod_equilibrium_params.f08
        dataIO/mod_exceptions.f08
        dataIO/mod_logging.f08
        mod_physical_constants.f08
        mod_units.f08
        mod_grid.f08
        mod_types.f08
        mod_check_values.f08
        solvers/mod_matrix_operations.f08
        mod_spline_functions.f08
        mod_interpolation.f08
        mod_integration.f08
        dataIO/mod_input.f08
        physics/mod_resistivity.f08
        physics/mod_cooling_curves.f08
        physics/mod_atmosphere_curves.f08
        physics/mod_solar_atmosphere.f08
        physics/mod_radiative_cooling.f08
        physics/mod_thermal_conduction.f08
        mod_inspections.f08
        mod_equilibrium.f08
        equilibria/smod_equil_adiabatic_homo.f08
        equilibria/smod_equil_constant_current.f08
        equilibria/smod_equil_coronal_flux_tube.f08
        equilibria/smod_equil_discrete_alfven.f08
        equilibria/smod_equil_flow_driven_instabilities.f08
        equilibria/smod_equil_gold_hoyle.f08
        equilibria/smod_equil_gravito_acoustic.f08
        equilibria/smod_equil_gravito_mhd.f08
        equilibria/smod_equil_interchange_modes.f08
        equilibria/smod_equil_internal_kink_instability.f08
        equilibria/smod_equil_isothermal_atmosphere.f08
        equilibria/smod_equil_kelvin_helmholtz_cd.f08
        equilibria/smod_equil_KHI.f08
        equilibria/smod_equil_magnetothermal_instabilities.f08
        equilibria/smod_equil_MRI_accretion.f08
        equilibria/smod_equil_photospheric_flux_tube.f08
        equilibria/smod_equil_resistive_homo.f08
        equilibria/smod_equil_resistive_tearing.f08
        equilibria/smod_equil_resistive_tearing_flow.f08
        equilibria/smod_equil_resonant_absorption.f08
        equilibria/smod_equil_rotating_plasma_cylinder.f08
        equilibria/smod_equil_RTI.f08
        equilibria/smod_equil_RTI_KHI.f08
        equilibria/smod_equil_RTI_theta_pinch.f08
        equilibria/smod_equil_suydam_cluster.f08
        equilibria/smod_equil_harris_sheet.f08
        ${USR_SMOD_LOC}/smod_user_defined.f08
        mod_make_subblock.f08
        mod_boundary_conditions.f08
        mod_matrix_creation.f08
        solvers/mod_solvers.f08
        solvers/mod_arpack_type.f08
        solvers/smod_qr_invert.f08
        solvers/smod_qz_direct.f08
        solvers/smod_arnoldi.f08
        mod_eigenfunctions.f08
        dataIO/mod_output.f08
        mod_hallmhd.f08
)

set_source_files_properties(${sources} PROPERTIES LANGUAGE Fortran)
set_source_files_properties(main.f08 PROPERTIES LANGUAGE Fortran)

set(LIBSTOLINK
    ${ARPACK_LIBRARIES}
    ${PARPACK_LIBRARIES}
#    ${MUMPS_LIBRARIES}
#    ${SCALAPACK_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

add_library(${LEGOLASLIB} ${sources})
add_executable(${LEGOLASEXEC} main.f08)
target_link_libraries(${LEGOLASEXEC} ${LEGOLASLIB} ${LIBSTOLINK})
