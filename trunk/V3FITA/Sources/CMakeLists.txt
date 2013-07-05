cmake_minimum_required (VERSION 2.8)

# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND v3fit_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bivariate.f
	${CMAKE_CURRENT_SOURCE_DIR}/boundary.f
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/coosig_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/diagnostic.f
	${CMAKE_CURRENT_SOURCE_DIR}/edge_limit_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/edge_limit_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/eq_interface.f
	${CMAKE_CURRENT_SOURCE_DIR}/eq_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/extcurz.f
	${CMAKE_CURRENT_SOURCE_DIR}/extcurz_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/geometric.f
	${CMAKE_CURRENT_SOURCE_DIR}/geometric_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/gsq_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/intpol.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipch_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipch_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/limiter.f
	${CMAKE_CURRENT_SOURCE_DIR}/limiter_iso_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/magnetic.f
	${CMAKE_CURRENT_SOURCE_DIR}/mddc_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/mmaw.f
	${CMAKE_CURRENT_SOURCE_DIR}/model_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/pprofile_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/prior.f
	${CMAKE_CURRENT_SOURCE_DIR}/prior_range.f
	${CMAKE_CURRENT_SOURCE_DIR}/prior_gaussian.f
	${CMAKE_CURRENT_SOURCE_DIR}/recon_param_model.f
	${CMAKE_CURRENT_SOURCE_DIR}/recon_param_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/result_output.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_cdf.f
        ${CMAKE_CURRENT_SOURCE_DIR}/signal_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/sxrem.f
	${CMAKE_CURRENT_SOURCE_DIR}/sxrch_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/sxrch_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/thomson.f
	${CMAKE_CURRENT_SOURCE_DIR}/thscte_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/thscte_mc.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3f_global.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fita.f
	${CMAKE_CURRENT_SOURCE_DIR}/data_parameters.f
	${CMAKE_CURRENT_SOURCE_DIR}/parameter.f
	${CMAKE_CURRENT_SOURCE_DIR}/model.f
	${CMAKE_CURRENT_SOURCE_DIR}/equilibrium.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_equilibrium.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fit_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fit_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/commandline_parser.f
	${CMAKE_CURRENT_SOURCE_DIR}/reconstruction.f
)
set (v3fit_sources "${v3fit_sources}" PARENT_SCOPE)
