cmake_minimum_required (VERSION 2.8)

# Create an empty variable to hold all the source files.
set (gtovmi_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# Define a library build from all the sources.
add_executable (xgtovmi ${gtovmi_sources})
set_target_properties (xgtovmi PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xgtovmi)

if (NOT (${NETCDF_C} STREQUAL "NETCDF_C-NOTFOUND"))
	target_link_libraries (xgtovmi ${NETCDF_C})
endif ()
if (NOT (${NETCDF_F} STREQUAL "NETCDF_F-NOTFOUND"))
	target_link_libraries (xgtovmi ${NETCDF_F})
endif ()

if (NOT PGPLOT_PATH)
	message (FATAL_ERROR "Must set -DPGPLOT_PATH")
endif ()
file (TO_NATIVE_PATH ${PGPLOT_PATH} PGPLOT_PATH)
find_library (PGPLOT_LIB NAMES pgplot PATHS ${PGPLOT_PATH})
find_file (PGPLOT_INC grpckg1.inc ${PGPLOT_PATH})
set_target_properties (xgtovmi PROPERTIES COMPILE_DEFINITIONS "PGPLOT_PATH")
get_target_property (GTOVMI_INLCUDES xgtovmi INCLUDE_DIRECTORIES)
set_target_properties (xgtovmi PROPERTIES INCLUDE_DIRECTORIES "${GTOVMI_INLCUDES};${PGPLOT_PATH}")
target_link_libraries (xgtovmi ${PGPLOT_LIB})

find_package(X11 REQUIRED)
if (X11_FOUND)
	target_link_libraries (xgtovmi ${X11_LIBRARIES})
endif ()
