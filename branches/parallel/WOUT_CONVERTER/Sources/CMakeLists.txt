# Append all source files to variable wout_converter_sources. As new files are added this must be updated.
list (APPEND wout_converter_sources
	${CMAKE_CURRENT_SOURCE_DIR}/wout_converter.f
)
set (wout_converter_sources "${wout_converter_sources}" PARENT_SCOPE)
