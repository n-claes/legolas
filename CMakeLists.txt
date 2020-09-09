cmake_minimum_required(VERSION 3.2)
project(Legolas)
set(VERSION 1.0)
enable_language(Fortran)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_FLAGS "-fcheck=all -fbounds-check -Wall \
                         -Wextra -Wconversion -pedantic -fbacktrace")

if (DEFINED ENV{LEGOLASDIR})
    set(SRC $ENV{LEGOLASDIR}/src)
    set(setup_tools $ENV{LEGOLASDIR}/setup_tools)
    set(BUILD $ENV{LEGOLASDIR}/build)
    file(MAKE_DIRECTORY ${BUILD})
else()
    message(WARNING
            "\nEnvironment variable LEGOLASDIR is not defined! \nTrying to use current source directory."
            )
    set(SRC ${CMAKE_SOURCE_DIR}/src)
    set(setup_tools "${CMAKE_SOURCE_DIR}/setup_tools")
    set(BUILD ${CMAKE_CURRENT_BINARY_DIR})
endif()
set(BIN ${BUILD}/bin)
set(OUTPUT ${CMAKE_SOURCE_DIR}/output)
# place .mod files in bin directory
set (CMAKE_Fortran_MODULE_DIRECTORY ${BIN})
# place executable in top directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
# place library in lib directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD}/lib)
# create directories
file(MAKE_DIRECTORY ${OUTPUT})

include(${setup_tools}/findBLAS.cmake)
include(${setup_tools}/findLAPACK.cmake)

# name of library
set(LEGOLASLIB lego)
# name of executable
set(LEGOLASEXEC legolas)
# location of default user-defined submodule
set(USR_SMOD_LOC equilibria)
add_subdirectory(${SRC} ${BIN})
