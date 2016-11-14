# Append all source files to variable bmw_sources. As new files are added this must be updated.
list (APPEND bmw_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bmw_commandline_parser.f
	${CMAKE_CURRENT_SOURCE_DIR}/bmw_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/m_grid.f
	${CMAKE_CURRENT_SOURCE_DIR}/primed_grid.f
	${CMAKE_CURRENT_SOURCE_DIR}/unprimed_grid.f
)
set (bmw_sources "${bmw_sources}" PARENT_SCOPE)
