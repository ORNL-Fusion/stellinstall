cmake_minimum_required (VERSION 2.8)

# Append all source files to variable xbootsj_sources. As new files are added this must be updated.
list (APPEND xbootsj_sources
	${CMAKE_CURRENT_SOURCE_DIR}/al31.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_angles.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_radial.f
	${CMAKE_CURRENT_SOURCE_DIR}/bongrid.f
	${CMAKE_CURRENT_SOURCE_DIR}/bootsj.f
	${CMAKE_CURRENT_SOURCE_DIR}/caprsh2.f
	${CMAKE_CURRENT_SOURCE_DIR}/datain.f
	${CMAKE_CURRENT_SOURCE_DIR}/deallocate_all.f
	${CMAKE_CURRENT_SOURCE_DIR}/denmf.f
	${CMAKE_CURRENT_SOURCE_DIR}/do_fft.f
	${CMAKE_CURRENT_SOURCE_DIR}/driver.f
	${CMAKE_CURRENT_SOURCE_DIR}/fraction.f
	${CMAKE_CURRENT_SOURCE_DIR}/grad.f
	${CMAKE_CURRENT_SOURCE_DIR}/othersums.f
	${CMAKE_CURRENT_SOURCE_DIR}/output.f
	${CMAKE_CURRENT_SOURCE_DIR}/parambs.f
	${CMAKE_CURRENT_SOURCE_DIR}/positiv.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_boozer.f
	${CMAKE_CURRENT_SOURCE_DIR}/reorganz.f
	${CMAKE_CURRENT_SOURCE_DIR}/SIMPUN.f
	${CMAKE_CURRENT_SOURCE_DIR}/smooth1.f
	${CMAKE_CURRENT_SOURCE_DIR}/sumit.f
	${CMAKE_CURRENT_SOURCE_DIR}/temp.f
	${CMAKE_CURRENT_SOURCE_DIR}/tok_fraction.f
	${CMAKE_CURRENT_SOURCE_DIR}/trig.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec0.f
	${CMAKE_CURRENT_SOURCE_DIR}/woflam.f
)
set (xbootsj_sources "${xbootsj_sources}" PARENT_SCOPE)
