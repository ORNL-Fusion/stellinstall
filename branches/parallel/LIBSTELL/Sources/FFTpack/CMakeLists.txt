# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/cfft99.f
	${CMAKE_CURRENT_SOURCE_DIR}/cftfax.f90
	${CMAKE_CURRENT_SOURCE_DIR}/cftrig.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fact.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fax.f
	${CMAKE_CURRENT_SOURCE_DIR}/fft99.f
	${CMAKE_CURRENT_SOURCE_DIR}/fft99a.f
	${CMAKE_CURRENT_SOURCE_DIR}/fft99b.f
	${CMAKE_CURRENT_SOURCE_DIR}/fft991.f
	${CMAKE_CURRENT_SOURCE_DIR}/fftfax.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fftrig.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fmn_to_uv.f
	${CMAKE_CURRENT_SOURCE_DIR}/trigfact.f
	${CMAKE_CURRENT_SOURCE_DIR}/vpassm.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
