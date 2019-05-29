# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND surface_sources
	${CMAKE_CURRENT_SOURCE_DIR}/surface.f
	${CMAKE_CURRENT_SOURCE_DIR}/surface_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/surface_commandline_parser.f
)
set (surface_sources "${surface_sources}" PARENT_SCOPE)
