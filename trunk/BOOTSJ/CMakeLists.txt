cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (xbootsj_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xbootsj ${xbootsj_sources})
add_dependencies (xbootsj stell)
target_link_libraries (xbootsj stell)
set_target_properties (xbootsj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xbootsj)
set_target_properties (xbootsj PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)

# HACK: Work around build failures on Linux. The set_target_properties INCLUDE_DIRECTORIES command is getting ignored.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set_target_properties (xbootsj PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif()

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
    target_link_libraries (xbootsj ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
    target_link_libraries (xbootsj ${NETCDF_F})
endif ()

