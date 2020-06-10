# Create an empty variable to hold all the source files.
set (vmec_interactive_sources "")

# Add subdirectory for all the sources.
add_subdirectory (Sources)

# This branch assumes that gRPC and all its dependencies are already installed
# on this system, so they can be located by find_package().

# https://stackoverflow.com/a/57224542
# protobuf weirdness since the cmake file for protobuf are not avail on Arch Linux :-/
# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 
# these calls create special `PkgConfig::<MODULE>` variables
pkg_check_modules(PROTOBUF_PKG REQUIRED IMPORTED_TARGET protobuf)

# https://github.com/grpc/grpc/blob/v1.29.1/examples/cpp/helloworld/CMakeLists.txt

set(_PROTOBUF_LIBPROTOBUF protobuf::libprotobuf)
set(_REFLECTION gRPC::grpc++_reflection)

if(CMAKE_CROSSCOMPILING)
  find_program(_PROTOBUF_PROTOC protoc)
else()
  set(_PROTOBUF_PROTOC $<TARGET_FILE:protobuf::protoc>)
endif()

# Find gRPC installation
# Looks for gRPCConfig.cmake file installed by gRPC's cmake installation.
find_package(gRPC CONFIG REQUIRED)
message(STATUS "Using gRPC ${gRPC_VERSION}")

set(_GRPC_GRPCPP gRPC::grpc++)
if(CMAKE_CROSSCOMPILING)
  find_program(_GRPC_CPP_PLUGIN_EXECUTABLE grpc_cpp_plugin)
else()
  set(_GRPC_CPP_PLUGIN_EXECUTABLE $<TARGET_FILE:gRPC::grpc_cpp_plugin>)
endif()

# Proto file
get_filename_component(vi_proto "./protos/vmec_interactive.proto" ABSOLUTE)
get_filename_component(vi_proto_path "${vi_proto}" PATH)

# Generated sources
set(vi_proto_srcs "${CMAKE_CURRENT_SOURCE_DIR}/Sources/vmec_interactive.pb.cc")
set(vi_proto_hdrs "${CMAKE_CURRENT_SOURCE_DIR}/Sources/vmec_interactive.pb.h")
set(vi_grpc_srcs "${CMAKE_CURRENT_SOURCE_DIR}/Sources/vmec_interactive.grpc.pb.cc")
set(vi_grpc_hdrs "${CMAKE_CURRENT_SOURCE_DIR}/Sources/vmec_interactive.grpc.pb.h")
add_custom_command(
      OUTPUT "${vi_proto_srcs}" "${vi_proto_hdrs}" "${vi_grpc_srcs}" "${vi_grpc_hdrs}"
      COMMAND ${_PROTOBUF_PROTOC}
      ARGS --grpc_out "${CMAKE_CURRENT_SOURCE_DIR}/Sources"
        --cpp_out "${CMAKE_CURRENT_SOURCE_DIR}/Sources"
        -I "${vi_proto_path}"
        --plugin=protoc-gen-grpc="${_GRPC_CPP_PLUGIN_EXECUTABLE}"
        "${vi_proto}"
      DEPENDS "${vi_proto}")

# Define an executable and link all libraries.
add_executable (xvmec_interactive ${vmec_interactive_sources} ${vi_proto_srcs} ${vi_grpc_srcs})
add_dependencies (xvmec_interactive vmec)
target_link_libraries (xvmec_interactive vmec bmw stell)
set_target_properties (xvmec_interactive PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/xvmec_interactive)

if (MPI_Fortran_FOUND)
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xvmec_interactive PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec;${MPI_Fortran_INCLUDE_PATH}")
	else ()
		set (mpi_include "")
		FOREACH (path ${MPI_Fortran_INLUCDE_PATH})
			set (mpi_include "${mpi_include} -I${path}")
		ENDFOREACH (path)
		set_target_properties (xvmec_interactive PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec ${mpi_inlcude}")
	endif ()
else ()
	if (CMAKE_VERSION VERSION_GREATER "2.8.8")
		set_target_properties (xvmec_interactive PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/stell;${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	else ()
		set_target_properties (xvmec_interactive PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/stell -I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	endif ()
endif ()

target_link_libraries (xvmec_interactive ${NETCDF_LIBRARIES})

# actually link against protobuf and grpc
target_link_libraries(xvmec_interactive PkgConfig::PROTOBUF_PKG)
target_link_libraries(xvmec_interactive ${_REFLECTION} ${_GRPC_GRPCPP})

if (SCALAPACK_LIB)
    target_link_libraries (xvmec_interactive ${SCALAPACK_LIB})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xvmec_interactive ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xvmec_interactive ${LAPACK_LIBRARIES})
endif ()
if (MPI_Fortran_FOUND)
    target_link_libraries (xvmec_interactive ${MPI_Fortran_LIBRARIES})
endif ()
