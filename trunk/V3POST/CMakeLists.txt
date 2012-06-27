cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3post_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xv3post ${v3post_sources})
add_dependencies (xv3post stell)
target_link_libraries (xv3post stell)

if (NETCDF_C AND NETCDF_F)
   target_link_libraries (xv3post ${NETCDF_C} ${NETCDF_F})
endif ()
