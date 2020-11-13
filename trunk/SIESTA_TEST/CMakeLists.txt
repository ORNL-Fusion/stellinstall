# Create an empty variable to hold all the source files.
set (siesta_test_sources "")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define an executable and link all libraries.
add_executable (xsiesta_test ${siesta_test_sources})
target_link_libraries (xsiesta_test ${NETCDF_LIBRARIES})
