#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.
add_test (NAME              te_pressure_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit -para=-1)

#  Since the reconstruction task test was run, Check the output files for correct values.

#  Check that the observed signal was written correctly.
add_test (NAME    te_pressure_test_check_observed_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=observed --value=1872.5 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=observed --value=1872.5 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=observed --value=1872.5 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=observed --value=1872.5 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_observed_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_observed_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_observed_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_observed_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=observed --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=observed --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=observed --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_observed_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=observed --value=400 --range=1.0E-10)

#  Check that the observed sigma was written correctly.
add_test (NAME    te_pressure_test_check_sigma_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=sigma --value=1.8725 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=sigma --value=1.8725 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=sigma --value=1.8725 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=sigma --value=1.8725 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    te_pressure_test_check_sigma_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    te_pressure_test_check_sigma_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    te_pressure_test_check_sigma_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    te_pressure_test_check_sigma_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=sigma --value=4.0 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=sigma --value=4.0 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=sigma --value=4.0 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_sigma_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=sigma --value=4.0 --range=1.0E-10)

#  Check that the model signal was computed correctly.
add_test (NAME    te_pressure_test_check_model_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=model --value=1872.5 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=model --value=1872.5 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=model --value=1872.5 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=model --value=1872.5 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    te_pressure_test_check_model_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=model --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_model_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=model --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_model_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=model --value=400 --range=1.0E-10)
add_test (NAME    te_pressure_test_check_model_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=model --value=400 --range=1.0E-10)
