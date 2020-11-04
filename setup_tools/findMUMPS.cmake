# This file is a modified version of the original one in the MUMPS repository:
# https://github.com/scivision/mumps/blob/master/cmake/Modules/FindMUMPS.cmake

function(unable_to_find)
    set(MUMPS_FOUND False PARENT_SCOPE)
    message(STATUS "Unable to find MUMPS!")
endfunction()


# check if mumps compiles
function(check_if_mumps_links)
    set(CMAKE_REQUIRED_INCLUDES
        ${MUMPS_INCLUDE_DIR}
    )
    set(CMAKE_REQUIRED_LIBRARIES
        ${MUMPS_LIBRARY}
        ${SCALAPACK_LIBRARIES}
        ${LAPACK_LIBRARIES}
    )
    include(CheckFortranSourceCompiles)
    foreach(comp IN LISTS _ariths)
        check_fortran_source_compiles("
            program test_mumps
                implicit none (type, external)
                include '${comp}mumps_struc.h'
                external :: ${comp}mumps
                type(${comp}mumps_struc) :: mumps_par
            end program"
            MUMPS_${comp}_links
            SRC_EXT
            f90
        )
        if(NOT ${MUMPS_${comp}_links})
            set(MUMPS_links False PARENT_SCOPE)
        endif()
    endforeach()
endfunction()

message(STATUS "Finding MUMPS")

# MUMPS needs scalapack, already checked before this script is called
if(NOT ${SCALAPACK_FOUND})
    message(STATUS "MUMPS needs SCALAPACK which was not found!")
    unable_to_find()
    return()
endif()

# check if MUMPS libraries are already defined
if(DEFINED ENV{MUMPS_LIBRARIES} AND DEFINED ENV{MUMPS_INCLUDE_DIR})
    message(STATUS "MUMPS variables have already been set")
    message(STATUS "MUMPS_LIBRARIES: $ENV{MUMPS_LIBRARIES}")
    message(STATUS "MUMPS_INCLUDE_DIR: $ENV{MUMPS_INCLUDE_DIR}")
    set(MUMPS_FOUND True)
    return()
endif()

set(MUMPS_FOUND False)

if(NOT DEFINED ENV{MUMPS_ROOT})
    message(
        STATUS "Environment variable $MUMPS_ROOT not set, searching default locations"
    )
endif()

# try to find MUMPS include dir
find_path(MUMPS_INCLUDE_DIR
    NAMES mumps_compat.h
    HINTS $ENV{MUMPS_ROOT}
    PATH_SUFFIXES include
)
if(MUMPS_INCLUDE_DIR)
    # if-statement is False if -NOTFOUND appended
    message(STATUS "Looking for MUMPS include directory - found")
    message(STATUS "Found in ${MUMPS_INCLUDE_DIR}")
else()
    message(STATUS "Looking for MUMPS include directory - not found")
    unable_to_find()
    return()
endif()

# try to find MUMPS common library
find_library(MUMPS_COMMON
    NAMES mumps_common
    HINTS $ENV{MUMPS_ROOT}
    PATH_SUFFIXES lib build/src
)
if(MUMPS_COMMON)
    message(STATUS "Looking for libmumps_common - found")
    message(STATUS "Found in ${MUMPS_COMMON}")
else()
    message(STATUS "Looking for libmumps_common - not found")
    unable_to_find()
    return()
endif()

# try to find PORD library
find_library(PORD
    NAMES pord
    HINTS $ENV{MUMPS_ROOT}
    PATH_SUFFIXES lib build/PORD
)
if(PORD)
    message(STATUS "Looking for libport - found")
    message(STATUS "Found in ${PORD}")
else()
    message(STATUS "Looking for libport - not found")
    unable_to_find()
    return()
endif()

# try to find MUMPS precision libraries
set(_ariths)
set(libsfound False)
foreach(comp s d c z)
    find_library(MUMPS_${comp}_lib
        NAMES ${comp}mumps
        HINTS $ENV{MUMPS_ROOT}
        PATH_SUFFIXES lib build/src
    )
    if(MUMPS_${comp}_lib)
        message(STATUS "Looking for lib${comp}mumps - found")
        set(libsfound True)
        list(APPEND _ariths ${comp})
    else()
        message(STATUS "Looking for lib${comp}mumps - not found")
    endif()

    set(MUMPS_${comp}_FOUND True)
    list(APPEND MUMPS_LIBRARY ${MUMPS_${comp}_lib})
endforeach()
# if no libs are found, abort
if(NOT ${libsfound})
    message(STATUS "Could not find any MUMPS libraries!")
    unable_to_find()
    return()
endif()

set(MUMPS_LIBRARY ${MUMPS_LIBRARY} ${MUMPS_COMMON} ${PORD})
set(MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR})

# check if we are able to link and compile Fortran programs containing MUMPS
set(MUMPS_links True)
check_if_mumps_links()
if(NOT MUMPS_links)
    message(STATUS "MUMPS libraries found but linking failed!")
    unable_to_find()
    return()
endif()

# find MUMPS package and set corresponding variables
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
    REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR MUMPS_links
    HANDLE_COMPONENTS
)

set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
if(NOT TARGET MUMPS::MUMPS)
    add_library(MUMPS::MUMPS INTERFACE IMPORTED)
    set_target_properties(MUMPS::MUMPS PROPERTIES
        INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
    )
endif()
