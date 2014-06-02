cmake_minimum_required (VERSION 2.8)

# Append all source files to variable xdescur_sources. As new files are added this must be updated.
list (APPEND xdescur_sources
	${CMAKE_CURRENT_SOURCE_DIR}/amplitud.f	${CMAKE_CURRENT_SOURCE_DIR}/descur.f	${CMAKE_CURRENT_SOURCE_DIR}/evolve.f	${CMAKE_CURRENT_SOURCE_DIR}/fftrans.f	${CMAKE_CURRENT_SOURCE_DIR}/fixaray.f	${CMAKE_CURRENT_SOURCE_DIR}/funct.f	${CMAKE_CURRENT_SOURCE_DIR}/getangle.f	${CMAKE_CURRENT_SOURCE_DIR}/getrz.f	${CMAKE_CURRENT_SOURCE_DIR}/inguess.f	${CMAKE_CURRENT_SOURCE_DIR}/order.f	${CMAKE_CURRENT_SOURCE_DIR}/printit.f	${CMAKE_CURRENT_SOURCE_DIR}/read_wout_curve.f	${CMAKE_CURRENT_SOURCE_DIR}/restart.f	${CMAKE_CURRENT_SOURCE_DIR}/scrunch.f	${CMAKE_CURRENT_SOURCE_DIR}/vname0.f	${CMAKE_CURRENT_SOURCE_DIR}/vname1.f
)
set (xdescur_sources "${xdescur_sources}" PARENT_SCOPE)