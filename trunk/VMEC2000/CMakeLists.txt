# Create an empty variable to hold all the source files.
set (vmec_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

add_library (vmec STATIC ${vmec_sources})
add_dependencies (vmec stell)
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (vmec PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
	set_target_properties (vmec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()
	

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/TimeStep/vmec.f)
add_dependencies (xvmec vmec)
target_link_libraries (xvmec vmec stell)

if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xvmec PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
else ()
	set_target_properties (xvmec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
endif ()

target_link_libraries (xvmec ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xvmec ${MPI_Fortran_LIBRARIES})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xvmec ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xvmec ${LAPACK_LIBRARIES})
endif ()

# build anmec
add_library (animec STATIC ${vmec_sources})
add_dependencies (animec stell)
set_target_properties (animec PROPERTIES COMPILE_DEFINITIONS "_ANIMEC")
set_target_properties (animec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/animec)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (animec PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
else ()
	set_target_properties (animec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell")
endif ()


# Define an executable and link all libraries.
add_executable (xanimec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/TimeStep/vmec.f)
add_dependencies (xanimec animec)
target_link_libraries (xanimec animec stell)


if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xanimec PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/animec")
else ()
	set_target_properties (xanimec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/animec")
endif ()

target_link_libraries (xanimec ${NETCDF_LIBRARIES})

if (MPI_Fortran_FOUND)
        target_link_libraries (xanimec ${MPI_Fortran_LIBRARIES})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xanimec ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xanimec ${LAPACK_LIBRARIES})
endif ()
