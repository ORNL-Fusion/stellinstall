cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/cfode.f
	${CMAKE_CURRENT_SOURCE_DIR}/ewset.f	${CMAKE_CURRENT_SOURCE_DIR}/intdy.f	${CMAKE_CURRENT_SOURCE_DIR}/lsode.f	${CMAKE_CURRENT_SOURCE_DIR}/prepj.f	${CMAKE_CURRENT_SOURCE_DIR}/solsy.f	${CMAKE_CURRENT_SOURCE_DIR}/stode.f	${CMAKE_CURRENT_SOURCE_DIR}/vnorm.f	${CMAKE_CURRENT_SOURCE_DIR}/xerrwv.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)