enable_testing ()

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

#  Run the reconstruction task. This just tests if the equilibrium executes correctly. 
#  Subsequent tests will examine the results of the equilibrium for correct 
#  results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)

add_subdirectory (tests)
