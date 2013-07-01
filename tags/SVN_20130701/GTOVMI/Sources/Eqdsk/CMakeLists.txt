cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/afile.f
	${CMAKE_CURRENT_SOURCE_DIR}/getbpsq.f
	${CMAKE_CURRENT_SOURCE_DIR}/getd3dcur.f
	${CMAKE_CURRENT_SOURCE_DIR}/gfile_allocate.f
	${CMAKE_CURRENT_SOURCE_DIR}/rdeqdsk.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
