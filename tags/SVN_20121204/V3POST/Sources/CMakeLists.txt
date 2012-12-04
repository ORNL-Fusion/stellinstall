cmake_minimum_required (VERSION 2.8)

# Add subdirectory for all the sources.
add_subdirectory (Modules)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND v3post_sources
	${CMAKE_CURRENT_SOURCE_DIR}/spline.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3post.f
)
set (v3post_sources "${v3post_sources}" PARENT_SCOPE)