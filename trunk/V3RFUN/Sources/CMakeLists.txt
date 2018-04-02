# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND v3rfun_sources
	${CMAKE_CURRENT_SOURCE_DIR}/diagnostic_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3rfun.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3rfun_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3rfun_context.f
)
set (v3rfun_sources "${v3rfun_sources}" PARENT_SCOPE)
