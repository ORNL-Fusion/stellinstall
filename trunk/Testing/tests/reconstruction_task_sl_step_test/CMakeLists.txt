#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.
add_test (NAME              reconstruction_task_sl_step_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Since the reconstruction task test was run, Check the recout file for correct values.
add_test (NAME    reconstruction_task_sl_step_test_pp_ne_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_ne_af --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_sl_step_test_am_aux_f
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=am_aux_f --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_sl_step_test_pp_te_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_te_af --value=15.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_sl_step_test_pp_sxrem_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_sxrem_af --value=1.0 --range=0.1 --index1=1)
add_test (NAME    reconstruction_task_sl_step_test_sem_totals
          COMMAND ${TEST_SCRIPT_DIR}/check_sem_total.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit)
add_test (NAME    reconstruction_task_sl_step_test_jac_1-4:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=0 --ihigh=3 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_1-4:3-3
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=0 --ihigh=3 --jlow=2 jhigh=2 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_5-16:3-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=4 --ihigh=15 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_17-20:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=16 --ihigh=19 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_17-20:4-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=16 --ihigh=19 --jlow=3 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_21-28:3-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=20 --ihigh=27 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_25-28:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=24 --ihigh=27 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=28 --ihigh=28 --jlow=1 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_2-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=29 --ihigh=29 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_2-2
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=29 --ihigh=29 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_3-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=30 --ihigh=30 --jlow=0 jhigh=1 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_3-2
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=30 --ihigh=30 --jlow=3 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_sl_step_test_jac_prior_4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=31 --ihigh=31 --jlow=0 jhigh=2 --value=0.0)
