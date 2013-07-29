enable_testing ()

#  Run the built in unit test task of v3fit.
add_test (NAME              unit_tests
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_tests 
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit -test)

#  Run the equilibrium task. This just tests if the equilibrium executes
#  correctly. Subsequent tests will examine the results of the equilibrium for
#  correct results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test)
add_test (NAME              equilibrium_task_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Run the v3post task. This just tests if the equilibrium executes correctly. 
#  Subsequent tests will examine the results of the equilibrium for correct 
#  results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test)
add_test (NAME              v3post_task_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Since the v3post task test was run, Check the recout file for correct values.
add_test (NAME    v3post_task_check_model_sxr001
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=sxr001 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr002
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=sxr002 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr003
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=sxr003 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_sxr004
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=sxr004 --attribute=model --value=0.5 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int001
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int001 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int002
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int002 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int003
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int003 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int004
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int004 --attribute=model --value=5.0E17 --range=1.0)
add_test (NAME    v3post_task_check_model_int005
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int005 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int006
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int006 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int007
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int007 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_int008
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=int008 --attribute=model --value=0 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm001
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm001 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm002
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm002 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm003
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm003 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm004
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm004 --attribute=model --value=15 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm005
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm005 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm006
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm006 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm007
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm007 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm008
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm008 --attribute=model --value=1.0E18 --range=1.0)
add_test (NAME    v3post_task_check_model_thm009
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm009 --attribute=model --value=400 --range=1.0E-10)
add_test (NAME    v3post_task_check_model_thm010
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm010 --attribute=model --value=400 --range=1.0)
add_test (NAME    v3post_task_check_model_thm011
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm011 --attribute=model --value=400 --range=1.0)
add_test (NAME    v3post_task_check_model_thm012
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/recout.input.test.v3fit --s_name=thm012 --attribute=model --value=400 --range=1.0)

#  Run the reconstruction task. This just tests if the equilibrium executes correctly. 
#  Subsequent tests will examine the results of the equilibrium for correct 
#  results.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test)
add_test (NAME              reconstruction_task_test
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Since the v3post task test was run, Check the recout file for correct values.
add_test (NAME    reconstruction_task_test_pp_ne_af
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/recout.input.test.v3fit --p_type=pp_ne_af --value=1.0 --range=0.1 --index1=2)
#add_test (NAME    reconstruction_task_test_am_aux_f
#          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/recout.input.test.v3fit --p_type=am_aux_f --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_test_pp_te_af
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/recout.input.test.v3fit --p_type=pp_te_af --value=15.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_test_pp_sxrem_af
          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/recout.input.test.v3fit --p_type=pp_sxrem_af --value=1.0 --range=0.1 --index1=1)
