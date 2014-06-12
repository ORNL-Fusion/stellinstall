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
set_target_properties (xbooz_xform PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)

# HACK: Work around build failures on Linux. The set_target_properties INCLUDE_DIRECTORIES command is getting ignored.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties (xbooz_xform PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif()

target_link_libraries (xbooz_xform ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xbooz_xform ${MPI_Fortran_LIBRARIES})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xbooz_xform ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xbooz_xform ${LAPACK_LIBRARIES})
endif ()
