cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (xbooz_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xbooz_xform ${xbooz_sources})
add_dependencies (xbooz_xform stell)
target_link_libraries (xbooz_xform stell)
set_target_properties (xbooz_xform PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xbooz_xform)
set_target_properties (xbooz_xform PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
    target_link_libraries (xbooz_xform ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
    target_link_libraries (xbooz_xform ${NETCDF_F})
endif ()

if (BLAS_FOUND)
    target_link_libraries (xbooz_xform ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
    target_link_libraries (xbooz_xform ${LAPACK_LIBRARIES})
endif ()
