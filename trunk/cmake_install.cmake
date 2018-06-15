# Install script for directory: /u/gscoig/STELLOPT_corrections/trunk

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
  include("/u/gscoig/STELLOPT_corrections/trunk/LIBSTELL/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/VMEC2000/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/POINCARE/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/SIESTA/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/V3FITA/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/MAKEGRID/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/V3RFUN/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/BOOZ_XFORM/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/WOUT_CONVERTER/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/BOOTSJ/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/DESCUR/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/BMW/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/LGRID/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/u/gscoig/STELLOPT_corrections/trunk/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
