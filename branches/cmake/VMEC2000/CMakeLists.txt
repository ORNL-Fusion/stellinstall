cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (vmec_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

add_library (vmec STATIC ${vmec_sources})
add_dependencies (vmec stell)

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/General/vmec_main.f)
add_dependencies (xvmec vmec)
target_link_libraries (xvmec vmec stell)

if (NETCDF_C AND NETCDF_F)
    target_link_libraries (xvmec ${NETCDF_C} ${NETCDF_F})
endif ()
