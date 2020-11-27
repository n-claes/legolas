message("")
message(STATUS "Finding SCALAPACK")

if(DEFINED ENV{SCALAPACK_LIBRARIES})
    message(STATUS "SCALAPACK found: $ENV{SCALAPACK_LIBRARIES}")
    set(SCALAPACK_FOUND TRUE)
    set(SCALAPACK_LIBRARIES $ENV{SCALAPACK_LIBRARIES})
    return()
endif()

find_library(
        SCALAPACK_LIBRARY
        NAMES scalapack scalapack-mpi scalapack-mpich scalapack-openmpi
        PATHS ${SCALAPACK_DIR} /usr/lib /usr/local /usr/lib64 /usr/local/lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK DEFAULT_MSG SCALAPACK_LIBRARY)
if (NOT ${SCALAPACK_FOUND})
    message(STATUS "Try setting $SCALAPACK_LIBRARIES or $SCALAPACK_DIR")
endif()
set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
