# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND gtovmi_sources
	${CMAKE_CURRENT_SOURCE_DIR}/cntourp.f
	${CMAKE_CURRENT_SOURCE_DIR}/dvdpsi.f
	${CMAKE_CURRENT_SOURCE_DIR}/equalarc.f
	${CMAKE_CURRENT_SOURCE_DIR}/errray1.f
#	${CMAKE_CURRENT_SOURCE_DIR}/furpl.f
	${CMAKE_CURRENT_SOURCE_DIR}/furplm.f
	${CMAKE_CURRENT_SOURCE_DIR}/garc.f
	${CMAKE_CURRENT_SOURCE_DIR}/mapout_nc.f
	${CMAKE_CURRENT_SOURCE_DIR}/mapperb.f
	${CMAKE_CURRENT_SOURCE_DIR}/pack.f
	${CMAKE_CURRENT_SOURCE_DIR}/rayatt.f
	${CMAKE_CURRENT_SOURCE_DIR}/rszs_allocate.f
)
set (gtovmi_sources "${gtovmi_sources}" PARENT_SCOPE)
