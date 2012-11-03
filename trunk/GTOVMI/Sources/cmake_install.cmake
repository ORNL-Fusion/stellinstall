# Install script for directory: /Users/mark/repo/trunk/GTOVMI/Sources

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Akima/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Eqdsk/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Ezcdf/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/General/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Graphics/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Map/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Modules/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Splining/cmake_install.cmake")
  INCLUDE("/Users/mark/repo/trunk/GTOVMI/Sources/Svd/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

