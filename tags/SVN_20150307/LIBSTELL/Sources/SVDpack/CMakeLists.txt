# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND stell_sources
	${CMAKE_CURRENT_SOURCE_DIR}/polysvd.f
	${CMAKE_CURRENT_SOURCE_DIR}/pythag.f
	${CMAKE_CURRENT_SOURCE_DIR}/SORTSVD.f
	${CMAKE_CURRENT_SOURCE_DIR}/SVBKSB.f
	${CMAKE_CURRENT_SOURCE_DIR}/svd_solve.f
	${CMAKE_CURRENT_SOURCE_DIR}/svdcmp.f
	${CMAKE_CURRENT_SOURCE_DIR}/svdinv.f
	${CMAKE_CURRENT_SOURCE_DIR}/svdinv2.f
	${CMAKE_CURRENT_SOURCE_DIR}/svdsolveb.f
	${CMAKE_CURRENT_SOURCE_DIR}/svdsolveuv.f
)
set (stell_sources "${stell_sources}" PARENT_SCOPE)
