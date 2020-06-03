# Append all source files to variable vmec_interactive_sources. As new files are added this must be updated.
list (APPEND vmec_interactive_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bivariate.f
	${CMAKE_CURRENT_SOURCE_DIR}/extcurz.f
	${CMAKE_CURRENT_SOURCE_DIR}/intpol.f
	${CMAKE_CURRENT_SOURCE_DIR}/ipch_dot.f
	${CMAKE_CURRENT_SOURCE_DIR}/pprofile_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/prior_gaussian.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal.f
	${CMAKE_CURRENT_SOURCE_DIR}/signal_dot.f
#	${CMAKE_CURRENT_SOURCE_DIR}/feedback.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fita.f
	${CMAKE_CURRENT_SOURCE_DIR}/data_parameters.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fit_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/model.f
	${CMAKE_CURRENT_SOURCE_DIR}/model_state.f
	${CMAKE_CURRENT_SOURCE_DIR}/equilibrium.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_equilibrium.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fit_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3fit_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_context.f
	${CMAKE_CURRENT_SOURCE_DIR}/commandline_parser.f
	${CMAKE_CURRENT_SOURCE_DIR}/reconstruction.f
	${CMAKE_CURRENT_SOURCE_DIR}/guassian_process.f
	${CMAKE_CURRENT_SOURCE_DIR}/combination.f
)
set (vmec_interactive_sources "${vmec_interactive_sources}" PARENT_SCOPE)
