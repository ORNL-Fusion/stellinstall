# Append all source files to variable vmec_interactive_sources. As new files are added this must be updated.
list (APPEND vmec_interactive_sources
	${CMAKE_CURRENT_SOURCE_DIR}/v3fita.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_equilibrium.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)
set (vmec_interactive_sources "${vmec_interactive_sources}" PARENT_SCOPE)
