#  Run the reconstruction task. This just tests if the reconstruction executes correctly.
#  Subsequent tests will examine the results of the reconstruction for correct
#  results.

add_test (NAME niter_clean_wout_test
          COMMAND rm -f wout_test.vmec.nc)

if (MPI_Fortran_FOUND)
	add_test (NAME    multigrid_test
	          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} ${EXECUTABLE_OUTPUT_PATH}/xvmec input.test.vmec)
else ()
	add_test (NAME    multigrid_test
	          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xvmec input.test.vmec)
endif ()

add_test (NAME niter_check_wout_test
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_file.py --file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc)
