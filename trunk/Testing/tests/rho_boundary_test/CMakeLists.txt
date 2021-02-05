#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.

add_test (NAME    rho_boundary_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit -para=-1)

#  Since the reconstruction task test was run, Check the output files for correct values.

#  Check that the reconstruction parameter reconstructed to the correct values.
add_test (NAME    rho_boundary_test_rho00
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.25 --range=0.0001 --index1=0 --index2=0)
add_test (NAME    rho_boundary_test_rho10
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.01 --range=0.0001 --index1=1 --index2=0)
add_test (NAME    rho_boundary_test_rho20
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.001 --range=0.0001 --index1=2 --index2=0)
add_test (NAME    rho_boundary_test_rho-21
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.005 --range=0.0001 --index1=-2 --index2=1)
add_test (NAME    rho_boundary_test_rho-11
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=-0.05 --range=0.0001 --index1=-1 --index2=1)
add_test (NAME    rho_boundary_test_rho01
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.003 --range=0.0001 --index1=0 --index2=1)
add_test (NAME    rho_boundary_test_rho11
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.1 --range=0.0001 --index1=1 --index2=1)
add_test (NAME    rho_boundary_test_rho21
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.01 --range=0.0001 --index1=2 --index2=1)
add_test (NAME    rho_boundary_test_rho_raxis_cc0
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rho_raxis_cc --value=1.0 --range=0.0001 --index1=0)
add_test (NAME    rho_boundary_test_rho_raxis_cc1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rho_raxis_cc --value=0.1 --range=0.0001 --index1=1)
add_test (NAME    rho_boundary_test_rho_raxis_cc2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rho_raxis_cc --value=0.01 --range=0.0001 --index1=2)
add_test (NAME    rho_boundary_test_rho_zaxis_cs1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rho_zaxis_cs --value=0.05 --range=0.0001 --index1=1)
add_test (NAME    rho_boundary_test_rho_zaxis_cs2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rho_zaxis_cs --value=0.005 --range=0.0001 --index1=2)

#  Check that the signal effectiveness total.
add_test (NAME    rho_boundary_test_sem_totals
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_sem_total.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit)

#  Check correlation matrix values.
add_test (NAME rho_boundary_test_check_reconstruction_correlation_matrix
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_correlation_matrix.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --param=Reconstruction)
