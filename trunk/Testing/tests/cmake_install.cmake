# Install script for directory: /u/gscoig/STELLOPT_corrections/trunk/Testing/tests

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
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/unit_tests/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_task_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/v3post_task_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/te_pressure_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/ne_pressure_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/reconstruction_task_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/reconstruction_task_seg_step_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/reconstruction_task_sl_step_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/param_constraints_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/v3post_lasym_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/z_offset_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_reset_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_reset_free_boundary_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/multigrid_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/aphi_test/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/lasym_boundary/cmake_install.cmake")
  include("/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/preconditioner_test/cmake_install.cmake")

endif()

