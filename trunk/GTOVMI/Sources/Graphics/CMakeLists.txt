# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/plotbox.f
	${CMAKE_CURRENT_SOURCE_DIR}/surfaces_plot.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
