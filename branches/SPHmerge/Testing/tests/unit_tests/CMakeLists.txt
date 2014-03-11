#  Run the built in unit test task of v3fit.
add_test (NAME    unit_tests 
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xv3fit -test)
