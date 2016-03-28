# Append all source files to variable siesta_sources. As new files are added this must be updated.
list (APPEND poincare_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bspline.f90
	${CMAKE_CURRENT_SOURCE_DIR}/enforce_periodicity.f
	${CMAKE_CURRENT_SOURCE_DIR}/jacobian.f
	${CMAKE_CURRENT_SOURCE_DIR}/lsode_quantities.f
	${CMAKE_CURRENT_SOURCE_DIR}/poincare.f
	${CMAKE_CURRENT_SOURCE_DIR}/reset_coords.f
	${CMAKE_CURRENT_SOURCE_DIR}/rt_hand_side_dphi.f
	${CMAKE_CURRENT_SOURCE_DIR}/rt_hand_side_dtheta.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_fields.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline3d_fit_coefs.f
)
set (poincare_sources "${poincare_sources}" PARENT_SCOPE)
