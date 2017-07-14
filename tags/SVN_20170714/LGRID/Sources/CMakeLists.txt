# Append all source files to variable lgrid_sources. As new files are added this must be updated.
list (APPEND lgrid_sources
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/vector3d.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/vector3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/vertex.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/vertex.cpp
)
set (lgrid_sources "${lgrid_sources}" PARENT_SCOPE)
