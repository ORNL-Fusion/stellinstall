# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND xbooz_sources
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_boozer.f
	${CMAKE_CURRENT_SOURCE_DIR}/boozer_coords.f
	${CMAKE_CURRENT_SOURCE_DIR}/boozer.f
	${CMAKE_CURRENT_SOURCE_DIR}/boozer_xform.f
	${CMAKE_CURRENT_SOURCE_DIR}/booz_jac.f
	${CMAKE_CURRENT_SOURCE_DIR}/booz_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/booz_persistent.f
	${CMAKE_CURRENT_SOURCE_DIR}/foranl.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_boozer.f
	${CMAKE_CURRENT_SOURCE_DIR}/harfun.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_wout_booz.f
	${CMAKE_CURRENT_SOURCE_DIR}/setup_booz.f
	${CMAKE_CURRENT_SOURCE_DIR}/transpmn.f
	${CMAKE_CURRENT_SOURCE_DIR}/trigfunc.f
	${CMAKE_CURRENT_SOURCE_DIR}/vcoords.f
	${CMAKE_CURRENT_SOURCE_DIR}/write_boozmn.f
)
set (xbooz_sources "${xbooz_sources}" PARENT_SCOPE)
