cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3rfun_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xv3rfun ${v3rfun_sources})
add_dependencies (xv3rfun stell)
target_link_libraries (xv3rfun stell)

if ((NOT ${NETCDF_C} STREQUAL "") AND (NOT ${NETCDF_F} STREQUAL ""))
   target_link_libraries (xv3rfun ${NETCDF_C} ${NETCDF_F})
endif ()
