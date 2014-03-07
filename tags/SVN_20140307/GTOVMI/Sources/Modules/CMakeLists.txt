cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/gfile.f
	${CMAKE_CURRENT_SOURCE_DIR}/mapg_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/mapout.f
	${CMAKE_CURRENT_SOURCE_DIR}/polynomial.f
	${CMAKE_CURRENT_SOURCE_DIR}/precision.f
	${CMAKE_CURRENT_SOURCE_DIR}/vname.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
