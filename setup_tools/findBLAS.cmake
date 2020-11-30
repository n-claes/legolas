message("")
message(STATUS "Finding BLAS")

if(DEFINED ENV{BLAS_LIBRARIES})
    message(STATUS "BLAS found: $ENV{BLAS_LIBRARIES}")
    set(BLAS_FOUND TRUE)
    set(BLAS_LIBRARIES $ENV{BLAS_LIBRARIES})
    return()
endif()

find_package(BLAS REQUIRED)
if (NOT ${BLAS_FOUND})
    message(STATUS "Try setting $BLAS_LIBRARIES")
endif()
