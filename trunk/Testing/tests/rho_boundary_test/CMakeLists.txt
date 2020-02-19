#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.

add_test (NAME    rho_boundary_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit -para=-1)

#  Since the reconstruction task test was run, Check the output files for correct values.

#  Check that the reconstruction parameter reconstructed to the correct values.
add_test (NAME    reconstruction_task_test_rho00
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.25 --range=0.0001 --index1=0 --index2=0)
add_test (NAME    reconstruction_task_test_rho10
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.01 --range=0.0001 --index1=1 --index2=0)
add_test (NAME    reconstruction_task_test_rho20
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.001 --range=0.0001 --index1=2 --index2=0)
add_test (NAME    reconstruction_task_test_rho-21
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.005 --range=0.0001 --index1=-2 --index2=1)
add_test (NAME    reconstruction_task_test_rho-11
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=-0.05 --range=0.0001 --index1=-1 --index2=1)
add_test (NAME    reconstruction_task_test_rho01
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.003 --range=0.0001 --index1=0 --index2=1)
add_test (NAME    reconstruction_task_test_rho11
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.1 --range=0.0001 --index1=1 --index2=1)
add_test (NAME    reconstruction_task_test_rho21
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.01 --range=0.0001 --index1=2 --index2=1)
add_test (NAME    reconstruction_task_test_rho_raxis_cc0
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=1.0 --range=0.0001 --index1=0)
add_test (NAME    reconstruction_task_test_rho_raxis_cc1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.1 --range=0.0001 --index1=1)
add_test (NAME    reconstruction_task_test_rho_raxis_cc2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.01 --range=0.0001 --index1=2)
add_test (NAME    reconstruction_task_test_rho_zaxis_cs1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.05 --range=0.0001 --index1=1)
add_test (NAME    reconstruction_task_test_rho_zaxis_cs2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rhoc --value=0.005 --range=0.0001 --index1=2)

#  Check that the signal effectiveness total.
add_test (NAME    reconstruction_task_test_sem_totals
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_sem_total.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit)

#  Check that the observed signal was written correctly.
add_test (NAME    reconstruction_task_test_check_observed_sxr001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr001 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_sxr002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr002 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_sxr003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr003 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_sxr004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr004 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int001 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_int002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int002 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_int003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int003 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_int004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int004 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_int005
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int005 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int006
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int006 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int007
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int007 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int008
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int008 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int009
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int009 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int010
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int010 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int011
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int011 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_int012
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int012 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm005
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm006
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm007
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm008
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm009
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=observed --value=400 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_thm010
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm011
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_thm012
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_test_check_observed_prior_1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_1 --attribute=observed --value=1 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_2 --attribute=observed --value=1 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_3
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_3 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_4
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_4 --attribute=observed --value=1 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_5
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_5 --attribute=observed --value=0.75 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_6
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_6 --attribute=observed --value=0.25 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_observed_prior_7
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_7 --attribute=observed --value=0.25 --range=1.0E-10)

#  Check that the observed sigma was written correctly.
add_test (NAME    reconstruction_task_test_check_sigma_sxr001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr001 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_sxr002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr002 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_sxr003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr003 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_sxr004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=sxr004 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int001 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_int002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int002 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_int003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int003 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_int004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int004 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_int005
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int005 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int006
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int006 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int007
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int007 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int008
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int008 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int009
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int009 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int010
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int010 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int011
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int011 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_int012
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=int012 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm001
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm002
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm003
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm004
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm005
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm006
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm007
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm008
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm009
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=sigma --value=4 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_thm010
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm011
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_thm012
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_test_check_sigma_prior_1
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_1 --attribute=sigma --value=1.00000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_2
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_2 --attribute=sigma --value=1.00000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_3
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_3 --attribute=sigma --value=1.50000E-01 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_4
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_4 --attribute=sigma --value=1.00000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_5
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_5 --attribute=sigma --value=0.0001 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_6
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_6 --attribute=sigma --value=0.0001 --range=1.0E-10)
add_test (NAME    reconstruction_task_test_check_sigma_prior_7
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --s_name=prior_7 --attribute=sigma --value=0.0001 --range=1.0E-10)

#  Check parameter constraint types for the correct values.
add_test (NAME    reconstruction_task_test_check_pp_ne_af_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=pp_ne_af --index1=2 --range_l=value)
add_test (NAME    reconstruction_task_test_check_am_aux_f_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=am_aux_f --index1=2 --range_l=value)
add_test (NAME    reconstruction_task_test_check_pp_te_af_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=pp_te_af --index1=2 --range_l=value)
add_test (NAME    reconstruction_task_test_check_pp_sxrem_af_a_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=pp_sxrem_af_a --index1=1 --index2=1 --range_l=value)
add_test (NAME    reconstruction_task_test_check_rbc00_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rbc --index1=0 --index2=0 --range_l=value --range_h=value)
add_test (NAME    reconstruction_task_test_check_rbc01_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=rbc --index1=0 --index2=1 --range_l=value)
add_test (NAME    reconstruction_task_test_check_zbs01_constraint_type
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_param_constraint_types.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --p_type=zbs --index1=0 --index2=1 --range_l=value)

#  Check correlation matrix values.
add_test (NAME reconstruction_task_test_check_reconstruction_correlation_matrix
          COMMAND ${PYTHON_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_correlation_matrix.py --file=${CMAKE_CURRENT_BINARY_DIR}/recout.input.test.v3fit --param=Reconstruction)
