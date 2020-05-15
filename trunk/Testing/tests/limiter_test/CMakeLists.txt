#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.

add_test (NAME    limiter_test_lgrid
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/lgrid limiter.dat 0.7 0.001 600 -0.30 0.001 600)
add_test (NAME    limiter_test_mgrid
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/mgrid coils.simple_kp36)

add_test (NAME    limiter_test_gridtest
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.grid.v3fit -para=-1)
add_test (NAME    limiter_test_gridtest_model
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.grid.v3fit --s_name=edge_lim_1 --attribute=model --value=0.00000E+00 --range=5.0E-9)

add_test (NAME    limiter_test_isotest
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.iso.v3fit -para=-1)
add_test (NAME    limiter_test_isotest_model
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.iso.v3fit --s_name=edge_lim_1 --attribute=model --value=0.00000E+00 --range=5.0E-9)

