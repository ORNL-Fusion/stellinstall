# Create an empty variable to hold all the source files.
set (wout_converter_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xwout_converter ${wout_converter_sources})
add_dependencies (xwout_converter stell)
target_link_libraries (xwout_converter stell)
set_target_properties (xwout_converter PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/wout_converter)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xwout_converter PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
else ()
	set_target_properties (xwout_converter PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()

target_link_libraries (xwout_converter ${NETCDF_LIBRARIES})

if (BLAS_FOUND)
	target_link_libraries (xwout_converter ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xwout_converter ${LAPACK_LIBRARIES})
endif ()
