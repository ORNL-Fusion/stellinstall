
add_test(mgrid_regtest1 ${CMAKE_CURRENT_SOURCE_DIR}/regressiontest.sh -t mgrid_regtest1 -e mgrid -d ${EXECUTABLE_OUTPUT_PATH})
add_test(vmec_regtest1 ${CMAKE_CURRENT_SOURCE_DIR}/regressiontest.sh -t vmec_regtest1 -e xvmec -d ${EXECUTABLE_OUTPUT_PATH})
add_test(v3fit_regtest1 ${CMAKE_CURRENT_SOURCE_DIR}/regressiontest.sh -t v3fit_regtest1 -e xv3fit -d ${EXECUTABLE_OUTPUT_PATH})

# Set timeouts for each test
set_tests_properties(mgrid_regtest1 PROPERTIES TIMEOUT 300)
set_tests_properties(vmec_regtest1 PROPERTIES TIMEOUT 300)
set_tests_properties(v3fit_regtest1 PROPERTIES TIMEOUT 300)
