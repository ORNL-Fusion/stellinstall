cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (mgrid_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (mgrid ${mgrid_sources})
add_dependencies (mgrid stell)
target_link_libraries (mgrid stell)

if (NETCDF_C AND NETCDF_F)
   target_link_libraries (mgrid ${NETCDF_C} ${NETCDF_F})
endif ()
