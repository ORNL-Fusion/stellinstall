# Create an empty variable to hold all the source files.
set (xbootsj_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xbootsj ${xbootsj_sources})
add_dependencies (xbootsj stell)
target_link_libraries (xbootsj stell)
set_target_properties (xbootsj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xbootsj)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xbootsj PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
	set_target_properties (xbootsj PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif()

target_link_libraries (xbootsj ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xbootsj ${MPI_Fortran_LIBRARIES})
endif ()
