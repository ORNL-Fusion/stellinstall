cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/li_gbfa.f	${CMAKE_CURRENT_SOURCE_DIR}/li_gbsl.f	${CMAKE_CURRENT_SOURCE_DIR}/li_gefa.f	${CMAKE_CURRENT_SOURCE_DIR}/li_gesl.f	${CMAKE_CURRENT_SOURCE_DIR}/liprec.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)