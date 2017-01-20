#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results. First the mdsig files must be made.
add_test (NAME    diagnostic_rotation_test_mdsig
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xv3rfun v3rfun.in)
add_test (NAME    diagnostic_rotation_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit -para=-1 -mpi=auto)

#  Since the reconstruction task test was run, Check the output files for correct values.

#  Check that the reconstruction parameter reconstructed to the correct values.
add_test (NAME    diagnostic_rotation_test_phi_offset
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=phi_offset --value=0.157 --range=0.0002)

#  Check that the signal effectiveness total.
add_test (NAME    diagnostic_rotation_test_sem_totals
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_sem_total.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit)

#  Check correlation matrix values.
add_test (NAME    diagnostic_rotation_test_check_reconstruction_correlation_matrix
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_correlation_matrix.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --param=Reconstruction)

#  Check the multiple sxr profile.
add_test (NAME    diagnostic_rotation_test_check_sxrem_model
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr1-005 --attribute=model --value=5.0E-1 --range=1.0E-1)
add_test (NAME    diagnostic_rotation_test_check_sxrem_model_2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr2-006 --attribute=model --value=1.0 --range=1.0E-1)
