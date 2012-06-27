cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/ajax_mod.f90
	${CMAKE_CURRENT_SOURCE_DIR}/linear1_mod.f90	${CMAKE_CURRENT_SOURCE_DIR}/spec_kind_mod.f90	${CMAKE_CURRENT_SOURCE_DIR}/spline1_mod.f90	${CMAKE_CURRENT_SOURCE_DIR}/track_mod.f90	${CMAKE_CURRENT_SOURCE_DIR}/write_mod.f90
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)