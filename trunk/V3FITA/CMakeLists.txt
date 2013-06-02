cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3fit_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xv3fit ${v3fit_sources})
add_dependencies (xv3fit vmec stell)
target_link_libraries (xv3fit vmec stell)
set_target_properties (xv3fit PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xv3fit)
set_target_properties (xv3fit PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
	target_link_libraries (xv3fit ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
	target_link_libraries (xv3fit ${NETCDF_F})
endif ()

if (BLAS_FOUND)
	target_link_libraries (xv3fit ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xv3fit ${LAPACK_LIBRARIES})
endif ()
