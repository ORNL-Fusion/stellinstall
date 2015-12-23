# Create an empty variable to hold all the source files.
set (poincare_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xpoincare ${poincare_sources})
add_dependencies (xpoincare stell)
target_link_libraries (xpoincare stell)
set_target_properties (xpoincare PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xpoincare)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xpoincare PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
	set_target_properties (xpoincare PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()

#target_link_libraries (xpoincare ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xpoincare ${MPI_Fortran_LIBRARIES})
endif ()
#if (LAPACK_FOUND)
#	target_link_libraries (xpoincare ${LAPACK_LIBRARIES})
#endif ()
