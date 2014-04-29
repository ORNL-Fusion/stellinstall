cmake_minimum_required (VERSION 2.8)

# This program only has tyhe driver. If it ever has more than one souce file uncomment the APPEND.

list (APPEND wout_converter_sources
	${CMAKE_CURRENT_SOURCE_DIR}/wout_converter.f
)
set (wout_converter_sources "${wout_converter_sources}" PARENT_SCOPE)
