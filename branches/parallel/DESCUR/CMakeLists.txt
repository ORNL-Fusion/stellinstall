cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (xdescur_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xdescur ${xdescur_sources})
add_dependencies (xdescur stell)
target_link_libraries (xdescur stell)
set_target_properties (xdescur PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xdescur)
set_target_properties (xdescur PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)

# HACK: Work around build failures on Linux. The set_target_properties INCLUDE_DIRECTORIES command is getting ignored.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties (xdescur PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif()

target_link_libraries (xdescur ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xdescur ${MPI_Fortran_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xdescur ${LAPACK_LIBRARIES})
endif ()
