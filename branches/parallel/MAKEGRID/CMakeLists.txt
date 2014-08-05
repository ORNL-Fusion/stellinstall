# Create an empty variable to hold all the source files.
set (mgrid_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (mgrid ${mgrid_sources})
add_dependencies (mgrid stell)
target_link_libraries (mgrid stell)
set_target_properties (mgrid PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/mgrid)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (mgrid PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
	set_target_properties (mgrid PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()

target_link_libraries (mgrid ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (mgrid ${MPI_Fortran_LIBRARIES})
endif ()
