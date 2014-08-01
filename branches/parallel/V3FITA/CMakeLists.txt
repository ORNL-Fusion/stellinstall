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
if (MPI_Fortran_FOUND)
	set_target_properties (xv3fit PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec;${MPI_Fortran_INCLUDE_PATH}")
else ()
	set_target_properties (xv3fit PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
endif ()

# HACK: Work around build failures on Linux. The set_target_properties INCLUDE_DIRECTORIES command is getting ignored.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	if (MPI_Fortran_FOUND)
		set_target_properties (xv3fit PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec -I${MPI_Fortran_INCLUDE_PATH}")
	else ()
                set_target_properties (xv3fit PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	endif ()
endif ()

target_link_libraries (xv3fit ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xv3fit ${MPI_Fortran_LIBRARIES})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xv3fit ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xv3fit ${LAPACK_LIBRARIES})
endif ()
