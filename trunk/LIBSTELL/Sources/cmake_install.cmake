# Install script for directory: /u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Coils/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Ezcdf/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/FFTpack/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/GMRes/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Linpack/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Lsode/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Miscel/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Modules/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Optimization/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Pspline/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/QMRpack/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/SMPINT/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/SVDpack/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/Sources/Track/cmake_install.cmake")

endif()

