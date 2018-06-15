# CMake generated Testfile for 
# Source directory: /u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_reset_test
# Build directory: /u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_reset_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(equilibrium_reset_test "/u/gscoig/STELLOPT_corrections/trunk/build/bin/xv3fit" "input.test.v3fit" "-para=-1")
add_test(equilibrium_reset_test_jacobian_iter "/mpcdf/soft/SLES122/common/anaconda/2/5.1.0/bin/python" "/u/gscoig/STELLOPT_corrections/trunk/Testing/test_utilities/check_jacobian_iter.py" "--file=/u/gscoig/STELLOPT_corrections/trunk/Testing/tests/equilibrium_reset_test/runlog.input.test.v3fit")
