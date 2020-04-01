# Create an empty variable to hold all the source files.
set (surface_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xsurface ${surface_sources})
add_dependencies (xsurface bmw stell)
target_link_libraries (xsurface bmw stell)
set_target_properties (xsurface PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xsurface)

if (CMAKE_VERSION VERSION_GREATER "2.8.8")
    set_target_properties (xsurface PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/bmw")
else ()
    set_target_properties (xsurface PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/bmw")
endif ()

target_link_libraries (xsurface ${NETCDF_LIBRARIES})
if (MPI_Fortran_FOUND)
    target_link_libraries (xsurface ${MPI_Fortran_LIBRARIES})
endif ()
if (LAPACK_FOUND)
    target_link_libraries (xv3fit ${LAPACK_LIBRARIES})
endif ()
