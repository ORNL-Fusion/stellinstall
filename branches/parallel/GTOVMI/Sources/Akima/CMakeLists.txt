cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/akima.f
	${CMAKE_CURRENT_SOURCE_DIR}/doakima.f
	${CMAKE_CURRENT_SOURCE_DIR}/ssort.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
