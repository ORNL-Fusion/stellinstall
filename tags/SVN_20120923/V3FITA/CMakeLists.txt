cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3fit_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xv3fit ${v3fit_sources})
add_dependencies (xv3fit vmec stell)
target_link_libraries (xv3fit vmec stell)

if ((NOT ${NETCDF_C} STREQUAL "") AND (NOT ${NETCDF_F} STREQUAL ""))
    target_link_libraries (xv3fit ${NETCDF_C} ${NETCDF_F})
elseif (NOT ${NETCDF_C} STREQUAL "")
    target_link_libraries (xv3fit ${NETCDF_C})
endif ()

if (BLAS AND LAPACK)
    target_link_libraries (xv3fit ${BLAS} ${LAPACK})
endif ()
