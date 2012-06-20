cmake_minimum_required (VERSION 2.8)

# Start Project and note the language used.
project (v3fit_suite Fortran)

# Add variable for the awk script.
set (awk_script Scripts/awk_cdir.awks)

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build/lib)

# Fortran specific settings. The fits setting tells the compiler to use the C preprocessor. The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS -cpp)
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/modules)

# Set Apple specific options.
if (APPLE)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D DARWIN")
    set (CMAKE_EXE_LINKER_FLAGS "-framework Accelerate")
endif ()

# Use netcdf if libraries are found.
set (NETCDF_PATH /Users/mark/NetCDF/)
find_library (NETCDF_C netcdf PATHS ${NETCDF_PATH}/lib)
find_library (NETCDF_F netcdff PATHS ${NETCDF_PATH}/lib)
find_file (NETCDF_I netcdf.inc ${NETCDF_PATH}/include)
if (NETCDF_C AND NETCDF_F AND NETCDF_I)
    message (STATUS "Using netcdf")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D NETCDF -I ${NETCDF_PATH}/include")
endif ()

# List all the subdirectories these are in library dependancy order.
add_subdirectory (LIBSTELL)
add_subdirectory (VMEC2000)
add_subdirectory (V3FITA)
