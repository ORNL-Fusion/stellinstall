# Create an empty variable to hold all the source files.
set (vmec_interactive_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xvmec_interactive ${vmec_interactive_sources})
add_dependencies (xvmec_interactive vmec)
target_link_libraries (xvmec_interactive vmec bmw stell)
set_target_properties (xvmec_interactive PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xvmec_interactive)

if (MPI_Fortran_FOUND)
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xvmec_interactive PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec;${MPI_Fortran_INCLUDE_PATH}")
	else ()
		set (mpi_include "")
		FOREACH (path ${MPI_Fortran_INLUCDE_PATH})
			set (mpi_include "${mpi_include} -I${path}")
		ENDFOREACH (path)
		set_target_properties (xvmec_interactive PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec ${mpi_inlcude}")
	endif ()
else ()
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xvmec_interactive PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	else ()
		set_target_properties (xvmec_interactive PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	endif ()
endif ()

target_link_libraries (xvmec_interactive ${NETCDF_LIBRARIES})

if (SCALAPACK_LIB)
    target_link_libraries (xvmec_interactive ${SCALAPACK_LIB})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xvmec_interactive ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xvmec_interactive ${LAPACK_LIBRARIES})
endif ()
if (MPI_Fortran_FOUND)
    target_link_libraries (xvmec_interactive ${MPI_Fortran_LIBRARIES})
endif ()
