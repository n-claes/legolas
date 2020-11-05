
function(unable_to_find)
    set(ARPACK_FOUND False PARENT_SCOPE)
    set(PARPACK_FOUND False PARENT_SCOPE)
    message(STATUS "Unable to find ARPACK!")
endfunction()

message("")
message(STATUS "Finding ARPACK")

set(ARPACK_FOUND False)
set(PARPACK_FOUND False)

if(DEFINED ENV{ARPACK_LIBRARIES})
    message(STATUS "ARPACK found: $ENV{ARPACK_LIBRARIES}")
    set(ARPACK_FOUND True)
    set(ARPACK_LIBRARIES $ENV{ARPACK_LIBRARIES})
    if(DEFINED ENV{PARPACK_LIBRARIES})
        message(STATUS "PARPACK found: $ENV{PARPACK_LIBRARIES}")
        set(PARPACK_FOUND True)
        set(PARPACK_LIBRARIES $ENV{PARPACK_LIBRARIES})
    endif()
    return()
endif()

if (NOT DEFINED ENV{ARPACK_ROOT})
    message(
        STATUS "Environment variable $ARPACK_ROOT not set, searching default locations"
    )
endif()

# try to find ARPACK include dir
find_path(ARPACK_INCLUDE_DIR
    NAMES arpackdef.h
    HINTS $ENV{ARPACK_ROOT}
    PATH_SUFFIXES include/arpack installed/include/arpack
)
if(ARPACK_INCLUDE_DIR)
    message(STATUS "Looking for ARPACK include directory - found")
    message(STATUS "Found in ${ARPACK_INCLUDE_DIR}")
else()
    message(STATUS "Looking for ARPACK include directory - not found")
    unable_to_find()
    return()
endif()

# try to find ARPACK library
find_library(ARPACK_LIBRARY
    NAMES arpack
    HINTS $ENV{ARPACK_ROOT}
    PATH_SUFFIXES lib installed/lib
)
if(ARPACK_LIBRARY)
    message(STATUS "Looking for libarpack - found")
else()
    message(STATUS "Looking for libarpack - not found")
    unable_to_find()
    return()
endif()

# try to find PARPACK library
find_library(PARPACK_LIBRARY
    NAMES parpack
    HINTS $ENV{ARPACK_ROOT}
    PATH_SUFFIXES lib installed/lib
)
if(PARPACK_LIBRARY)
    message(STATUS "Looking for libparpack - found")
else()
    # don't abort if not found
    message(STATUS "Looking for libparpack - not found")
endif()

# find (P)ARPACK package and set corresponding variables
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK
    REQUIRED_VARS ARPACK_LIBRARY ARPACK_INCLUDE_DIR
    HANDLE_COMPONENTS
)
find_package_handle_standard_args(PARPACK
    REQUIRED_VARS PARPACK_LIBRARY
    HANDLE_COMPONENTS
)

set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
set(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
set(ARPACK_INCLUDE_DIRS ${ARPACK_INCLUDE_DIR})
if(NOT TARGET ARPACK::ARPACK)
    add_library(ARPACK::ARPACK INTERFACE IMPORTED)
    set_target_properties(ARPACK::ARPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES ${ARPACK_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIR}
    )
endif()
if(PARPACK_FOUND AND NOT TARGET PARPACK::PARPACK)
    add_library(PARPACK::PARPACK INTERFACE IMPORTED)
    set_target_properties(PARPACK::PARPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES ${PARPACK_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIR}
    )
endif()
