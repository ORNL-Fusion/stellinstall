enable_testing ()

#  The python interpreter may not be in the same location on all systems.
find_package (PythonInterp)

set (TEST_SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities)

#  Copy common files to the equilibrium  test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test)

#  Copy common files to the v3post test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)

#  Copy common files to the reconstruction test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)

#  Copy common files to the reconstruction seg step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_seg_step_test)

#  Copy common files to the reconstruction sl step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_sl_step_test)

#  Copy common files to the parameter constraints test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/param_constraints_test)

#  Copy common files to the te pressure test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/te_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/te_pressure_test)

#  Copy common files to the ne pressure test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/ne_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/ne_pressure_test)

#  Copy common files to the v3post lasym test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_lasym_test)

#  Copy common files to the reconstruction test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/z_offset_test)

add_subdirectory (tests)
