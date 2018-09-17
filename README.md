# stellinstall Quick Start Quide

## Introduction
This guide is meant as a quick reference to get the code compiled, and generate documentation. In depth documentation will be provided in the Doxygen generated documentation. Pre generated user and developer documentation can be found at [Doxygen Pages](https://ornl-fusion.github.io/stellinstall/).

## Requirements

1. NetCDF library: Availble from your operating system package management or from [Unidata](http://www.unidata.ucar.edu).
2. MPI library: Any mpi library will work. For local machines there is [OpenMPI](https://www.open-mpi.org).
3. LAPACK/BLAS library: Should be available with your operating system or compiler.
4. SCALAPACK library: Availble from your operating system package management or from [SCALAPACK](http://www.netlib.org/scalapack/).
5. CMake: Available from your operating system’s package management or from [cmake.org](https://cmake.org).

### Optional
1. Doxygen: Available from your operating system’s package management or from [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

## Building

To build V3FIT, use the setup_cmake script. If the script does not support your setup you will either need to add support or configure cmake manually. To configure cmake, at a minimum, the `NETCDF_LIB_PATH` and `NETCDF_INC_PATH` variables must be defined with the paths to the NetCDF library and include directories. To configure these, use

`cmake –DNETCDF_INC_PATH=/path/to/include`
`      –DNETCDF_LIB_PATH=/path/to/lib`

Inside project root directory.
When the NetCDF library’s and includes have been found cmake should report, “Using netcdf”. If the configuration fails Cmake should provide warning. A likely common case of failure is cmake detecting an unsupported compiler. The Cmake build system is known to support the `ifort` and `gfortran` compilers. If such an event occurs you can force the Fortran compiler to use by setting the `CMAKE_Fortran_COMPILER` variable. Once Cmake has been successfully configured, the code can be built using the generated make file. Built executables can be found in the `build/bin` directory. If Cmake fails to find libraries or libraries are found in non-standard paths, setting the `CMAKE_PREFIX_PATH` variable can fix these problems.
If the Doxygen tool has been installed, more extensive user and programming documentation may be generated. When Cmake finds the Doxygen tool, a new target in the makefile is created. Documentation may be generated using the

`make doc`

target. HTML and RTF based documents will be created in `Documentation/internals` directory. HTML may be view by opening the `html/index.html` file in a web browser. RTF documents may be found in the `rtf/refman.rtf` file.
