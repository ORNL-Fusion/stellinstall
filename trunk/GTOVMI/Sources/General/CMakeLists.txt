# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/deriv.f
	${CMAKE_CURRENT_SOURCE_DIR}/descur_sub.f
	${CMAKE_CURRENT_SOURCE_DIR}/gtovmi2.f
	${CMAKE_CURRENT_SOURCE_DIR}/initpsi.f
	${CMAKE_CURRENT_SOURCE_DIR}/shiftx.f
	${CMAKE_CURRENT_SOURCE_DIR}/sortstuff.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
