#  Run the v3post task. This just tests if the equilibrium executes correctly. 
#  Subsequent tests will examine the results of the v3post for correct 
#  results.
add_test (NAME              v3post_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit -para=-1)

#  Since the v3post task test was run, Check the recout file for correct values.
add_test (NAME    v3post_task_check_model_sxr001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr001 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr002 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr003 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr004 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int001 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int002 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int003 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int004 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int005 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int006 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int007 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int008 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int009 --attribute=model --value=-5.5E-5 --range=5.0E-7)
add_test (NAME    v3post_task_check_model_int010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int010 --attribute=model --value=-5.5E-5 --range=5.0E-7)
add_test (NAME    v3post_task_check_model_int011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int011 --attribute=model --value=-5.5E-5 --range=5.0E-7)
add_test (NAME    v3post_task_check_model_int012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int012 --attribute=model --value=-5.5E-5 --range=5.0E-7)
add_test (NAME    v3post_task_check_model_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=model --value=400 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=model --value=400 --range=1.0)
add_test (NAME    v3post_task_check_model_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=model --value=400 --range=1.0)
add_test (NAME    v3post_task_check_model_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=model --value=400 --range=1.0)
add_test (NAME    v3post_task_check_model_prior_1
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_1 --attribute=model --value=1 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_prior_2
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_2 --attribute=model --value=1 --range=1.0)
add_test (NAME    v3post_task_check_model_prior_3
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_3 --attribute=model --value=15 --range=1.0)
add_test (NAME    v3post_task_check_model_prior_4
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_4 --attribute=model --value=1 --range=1.0)
