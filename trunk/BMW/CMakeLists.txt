# Create an empty variable to hold all the source files.
set (bmw_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

add_library (bmw STATIC ${bmw_sources})
add_dependencies (bmw stell)
set_target_properties (bmw PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/bmw)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
    set_target_properties (bmw PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
    set_target_properties (bmw PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()

# Define an executable and link all libraries.
add_executable (xbmw ${CMAKE_CURRENT_SOURCE_DIR}/Sources/bmw.f)
add_dependencies (xbmw bmw)
target_link_libraries (xbmw bmw stell)

if (CMAKE_VERSION VERSION_GREATER "2.8.8")
    set_target_properties (xbmw PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/bmw")
else ()
    set_target_properties (xbmw PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/bmw")
endif ()

target_link_libraries (xbmw ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
    target_link_libraries (xbmw ${MPI_Fortran_LIBRARIES})
endif ()
