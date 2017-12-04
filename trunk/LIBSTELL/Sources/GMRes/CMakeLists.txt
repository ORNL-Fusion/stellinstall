# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/dPackgmres.f
	${CMAKE_CURRENT_SOURCE_DIR}/gmres_lib.f
	${CMAKE_CURRENT_SOURCE_DIR}/gmresr.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
