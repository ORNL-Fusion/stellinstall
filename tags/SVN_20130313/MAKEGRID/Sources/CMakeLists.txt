cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND mgrid_sources
	${CMAKE_CURRENT_SOURCE_DIR}/makegrid_global_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/makegrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/write_mgrid.f
)
set (mgrid_sources "${mgrid_sources}" PARENT_SCOPE)