cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/init_bg_currents.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_modular_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_modular_currents.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_modular_wsurf.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_saddle_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_saddle_currents.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_saddle_wsurf.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_tf_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_vf_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_vf_currents.f
	${CMAKE_CURRENT_SOURCE_DIR}/initialize_coilsin.f
	${CMAKE_CURRENT_SOURCE_DIR}/load_bg_currents.f
	${CMAKE_CURRENT_SOURCE_DIR}/load_modular_coils.f	${CMAKE_CURRENT_SOURCE_DIR}/load_modular_currents.f	${CMAKE_CURRENT_SOURCE_DIR}/load_modular_wsurf.f	${CMAKE_CURRENT_SOURCE_DIR}/load_saddle_coils.f	${CMAKE_CURRENT_SOURCE_DIR}/load_saddle_currents.f	${CMAKE_CURRENT_SOURCE_DIR}/load_saddle_wsurf.f	${CMAKE_CURRENT_SOURCE_DIR}/load_tf_coils.f	${CMAKE_CURRENT_SOURCE_DIR}/load_vf_coils.f	${CMAKE_CURRENT_SOURCE_DIR}/load_vf_currents.f	${CMAKE_CURRENT_SOURCE_DIR}/write_coilsin.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
