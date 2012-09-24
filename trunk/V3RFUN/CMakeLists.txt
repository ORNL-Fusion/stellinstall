cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3rfun_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xv3rfun ${v3rfun_sources})
add_dependencies (xv3rfun stell)
target_link_libraries (xv3rfun stell)

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
    target_link_libraries (xv3rfun ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
    target_link_libraries (xv3rfun ${NETCDF_F})
endif ()

if (LAPACK_FOUND)
    target_link_libraries (xv3rfun ${LAPACK_LIBRARIES})
endif ()
