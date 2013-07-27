enable_testing ()

#  Run the built in unit test task of v3fit.
add_test (NAME              unit_tests
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_tests 
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit -test)

#  Run the equilibrium task. This just tests if the equilibrium executes
#  correctly. Subsequent tests will examine the results of the equilibrium for
#  correct results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test)
add_test (NAME              equilibrium_task_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Run the v3post task. This just tests if the equilibrium executes correctly. 
#  Subsequent tests will examine the results of the equilibrium for correct 
#  results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
add_test (NAME              v3post_task_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)
