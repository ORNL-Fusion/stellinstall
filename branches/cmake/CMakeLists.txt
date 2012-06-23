# cmake build system was created by Mark Cianciosa. Current this only supports Mac OS X.
#
# All options can be define by using -DOPTION_NAME="Definition of option"
#
# Options:
# EXTRA_FLAGS         -- Compiler flags added to all builds.
# EXTRA_RELEASE_FLAGS -- Compiler flags added to release builds.
# EXTRA_DEBUG_FLAGS   -- Compiler flags added to debug builds.
#
# NETCDF_LIB_PATH     -- Path to the netcdf libraries.
# NETCDF_INC_PATH     -- Path to the netcdf include files.

cmake_minimum_required (VERSION 2.8)

# Start Project and note the language used.
project (v3fit_suite Fortran)

# Add variable for the awk script.
set (awk_script Scripts/awk_cdir.awks)

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)

# Add extra compiler flags to various builds.
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_DEBUG_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS}")

# Fortran specific settings. The fits setting tells the compiler to use the C preprocessor. The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/modules)

# Set Apple specific options. This add a predefine for DARWIN and adds lapack from the Accelerate framework.
if (APPLE)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D DARWIN")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
endif ()

# Use netcdf if libraries are found. The values of NETCDF_LIB_PATH and NETCDF_INC_PATH should be set from the command line or Cmake GUI. Convert path from native format to cmake internal format.
file (TO_NATIVE_PATH ${NETCDF_LIB_PATH} NETCDF_LIB_PATH)
file (TO_NATIVE_PATH ${NETCDF_LIB_PATH} NETCDF_LIB_PATH)
find_library (NETCDF_C netcdf PATHS ${NETCDF_LIB_PATH})
find_library (NETCDF_F netcdff PATHS ${NETCDF_LIB_PATH})
find_file (NETCDF_I netcdf.inc ${NETCDF_INC_PATH})
if (NETCDF_C AND NETCDF_F AND NETCDF_I)
    message (STATUS "Using netcdf")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D NETCDF -I ${NETCDF_INC_PATH}")
endif ()

# List all the subdirectories these are in library dependancy order.
add_subdirectory (LIBSTELL)
add_subdirectory (VMEC2000)
add_subdirectory (V3FITA)
