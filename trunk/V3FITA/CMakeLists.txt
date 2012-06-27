cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3fit_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xv3fit ${v3fit_sources})
add_dependencies (xv3fit vmec stell)
target_link_libraries (xv3fit vmec stell)

if (NETCDF_C AND NETCDF_F)
    target_link_libraries (xv3fit ${NETCDF_C} ${NETCDF_F})
endif ()
