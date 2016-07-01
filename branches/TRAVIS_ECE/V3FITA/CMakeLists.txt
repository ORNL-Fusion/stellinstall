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
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xv3fit PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec;${MPI_Fortran_INCLUDE_PATH}")
	else ()
		set (mpi_include "")
		FOREACH (path ${MPI_Fortran_INLUCDE_PATH})
			set (mpi_include "${mpi_include} -I${path}")
		ENDFOREACH (path)
		set_target_properties (xv3fit PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec ${mpi_inlcude}")		
	endif ()
else ()
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xv3fit PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
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
if (SCALAPACK_LIB)
        target_link_libraries (xv3fit ${SCALAPACK_LIB})
endif ()
