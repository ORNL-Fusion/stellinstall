# CMake generated Testfile for 
# Source directory: /u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test
# Build directory: /u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(diagnostic_rotation_test_mdsig "/u/gscoig/STELLOPT_corrections/trunk/build/bin/xv3rfun" "v3rfun.in")
add_test(diagnostic_rotation_test "/u/gscoig/STELLOPT_corrections/trunk/build/bin/xv3fit" "input.test.v3fit" "-para=-1" "-mpi=auto")
add_test(diagnostic_rotation_test_phi_offset "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_param_value.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/recout.input.test.v3fit" "--p_type=phi_offset" "--value=0.157" "--range=0.0002")
add_test(diagnostic_rotation_test_sem_totals "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_sem_total.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/recout.input.test.v3fit")
add_test(diagnostic_rotation_test_check_reconstruction_correlation_matrix "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_correlation_matrix.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/recout.input.test.v3fit" "--param=Reconstruction")
add_test(diagnostic_rotation_test_check_sxrem_model "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_signal_value.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/recout.input.test.v3fit" "--s_name=sxr1-005" "--attribute=model" "--value=5.0E-1" "--range=1.0E-1")
add_test(diagnostic_rotation_test_check_sxrem_model_2 "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_signal_value.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/diagnostic_rotation_test/recout.input.test.v3fit" "--s_name=sxr2-006" "--attribute=model" "--value=1.0" "--range=1.0E-1")
