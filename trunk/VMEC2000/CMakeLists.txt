cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (vmec_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

add_library (vmec STATIC ${vmec_sources})
add_dependencies (vmec stell)

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/TimeStep/vmec.f)
add_dependencies (xvmec vmec)
target_link_libraries (xvmec vmec stell)

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
    target_link_libraries (xvmec ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
    target_link_libraries (xvmec ${NETCDF_F})
endif ()

if (BLAS_FOUND)
    target_link_libraries (xvmec ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
    target_link_libraries (xvmec ${LAPACK_LIBRARIES})
endif ()

# build anmec
add_library (animec STATIC ${vmec_sources})
# Create hack to make sure the vmec modules get built before animec.
add_dependencies (animec stell vmec)
set_target_properties(animec PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS} -D_ANIMEC")

# Define an executable and link all libraries.
add_executable (xanimec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/TimeStep/vmec.f)
add_dependencies (xanimec animec)
target_link_libraries (xanimec animec stell)

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
    target_link_libraries (xanimec ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
    target_link_libraries (xanimec ${NETCDF_F})
endif ()

if (BLAS_FOUND)
    target_link_libraries (xanimec ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
    target_link_libraries (xanimec ${LAPACK_LIBRARIES})
endif ()
