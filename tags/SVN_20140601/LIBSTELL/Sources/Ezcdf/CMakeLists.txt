cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_attrib.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenGet.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenPut.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_inqvar.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_opncls.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf.f90
	${CMAKE_CURRENT_SOURCE_DIR}/handle_err.f90
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
