cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (stell_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_library (stell STATIC ${stell_sources})

