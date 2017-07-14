# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/copy_file.c
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
