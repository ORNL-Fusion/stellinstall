cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (v3rfun_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xv3rfun ${v3rfun_sources})
add_dependencies (xv3rfun stell)
target_link_libraries (xv3rfun stell)
set_target_properties (xv3rfun PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xv3rfun)
set_target_properties (xv3rfun PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)

# HACK: Work around build failures on Linux. The set_target_properties INCLUDE_DIRECTORIES command is getting ignored.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties (xv3rfun PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif()

target_link_libraries (xv3rfun ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xv3rfun ${MPI_Fortran_LIBRARIES})
endif ()
if (LAPACK_FOUND)
    target_link_libraries (xv3rfun ${LAPACK_LIBRARIES})
endif ()
