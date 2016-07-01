# Create an empty variable to hold all the source files.
set (stell_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_library (stell STATIC ${stell_sources})
set_target_properties (stell PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/stell)
