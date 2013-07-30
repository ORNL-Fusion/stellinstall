#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.
add_test (NAME              reconstruction_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Since the reconstruction task test was run, Check the recout file for correct values.
add_test (NAME    reconstruction_task_test_pp_ne_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_ne_af --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_test_am_aux_f
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=am_aux_f --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_test_pp_te_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_te_af --value=15.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_test_pp_sxrem_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_sxrem_af --value=1.0 --range=0.1 --index1=1)
