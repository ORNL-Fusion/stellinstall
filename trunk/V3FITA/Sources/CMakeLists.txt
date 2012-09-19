cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND v3fit_sources
	${CMAKE_CURRENT_SOURCE_DIR}/b_transform_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/b_transform_vmec.f
	${CMAKE_CURRENT_SOURCE_DIR}/bfield_grid_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/bivariate.f
	${CMAKE_CURRENT_SOURCE_DIR}/coosig_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/density_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/edge_limit_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/edge_limit_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/eq_interface.f
	${CMAKE_CURRENT_SOURCE_DIR}/eq_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/extcurz_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/faraday_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/geometric_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/gsq_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/ip_global.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipch_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipch_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipsl_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/limiter_iso_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/mddc_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/mmaw.f
	${CMAKE_CURRENT_SOURCE_DIR}/model_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/plasma_edge.f
	${CMAKE_CURRENT_SOURCE_DIR}/pprofile_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/ptr_links.f
	${CMAKE_CURRENT_SOURCE_DIR}/recon_param_model.f
	${CMAKE_CURRENT_SOURCE_DIR}/recon_param_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/result_output_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_cdf.f
        ${CMAKE_CURRENT_SOURCE_DIR}/signal_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/sxrch_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/sxrch_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/thscte_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/thscte_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3f_global.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fita.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3read_wout.f
)
set (v3fit_sources "${v3fit_sources}" PARENT_SCOPE)
