#  Run the reconstruction task. This just tests if the reconstruction executes correctly. 
#  Subsequent tests will examine the results of the reconstruction for correct 
#  results.
add_test (NAME              reconstruction_task_seg_step_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)

#  Since the reconstruction task test was run, Check the output files for correct values.

#  Check that the reconstruction parameter reconstructed to the correct values.
add_test (NAME    reconstruction_task_seg_step_test_pp_ne_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_ne_af --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_seg_step_test_am_aux_f
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=am_aux_f --value=1.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_seg_step_test_pp_te_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_te_af --value=15.0 --range=0.1 --index1=2)
add_test (NAME    reconstruction_task_seg_step_test_pp_sxrem_af
          COMMAND ${TEST_SCRIPT_DIR}/check_param_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --p_type=pp_sxrem_af --value=1.0 --range=0.1 --index1=1)
add_test (NAME    reconstruction_task_seg_step_test_sem_totals
          COMMAND ${TEST_SCRIPT_DIR}/check_sem_total.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit)

#  Check that jacobians for the the correct values.
add_test (NAME    reconstruction_task_seg_step_test_jac_1-4:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=0 --ihigh=3 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_1-4:3-3
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=0 --ihigh=3 --jlow=2 jhigh=2 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_5-16:3-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=4 --ihigh=15 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_17-20:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=16 --ihigh=19 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_17-20:4-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=16 --ihigh=19 --jlow=3 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_21-28:3-4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=20 --ihigh=27 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_25-28:1-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=24 --ihigh=27 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=28 --ihigh=28 --jlow=1 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_2-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=29 --ihigh=29 --jlow=0 jhigh=0 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_2-2
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=29 --ihigh=29 --jlow=2 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_3-1
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=30 --ihigh=30 --jlow=0 jhigh=1 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_3-2
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=30 --ihigh=30 --jlow=3 jhigh=3 --value=0.0)
add_test (NAME    reconstruction_task_seg_step_test_jac_prior_4
          COMMAND ${TEST_SCRIPT_DIR}/check_jac_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/runlog.input.test.v3fit --ilow=31 --ihigh=31 --jlow=0 jhigh=2 --value=0.0)

#  Check that the observed signal was written correctly.
add_test (NAME    reconstruction_task_seg_step_test_check_observed_sxr001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr001 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_sxr002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr002 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_sxr003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr003 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_sxr004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr004 --attribute=observed --value=0.5 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int001 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int002 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int003 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int004 --attribute=observed --value=5.0E17 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int005 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int006 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int007 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int008 --attribute=observed --value=0 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int009 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int010 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int011 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_int012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int012 --attribute=observed --value=-5.52983E-05 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=observed --value=15 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=observed --value=1.0E18 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=observed --value=400 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=observed --value=400 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_prior_1
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_1 --attribute=observed --value=1 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_prior_2
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_2 --attribute=observed --value=1 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_prior_3
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_3 --attribute=observed --value=15 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_observed_prior_4
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_4 --attribute=observed --value=1 --range=1.0)

#  Check that the observed sigma was written correctly.
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_sxr001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr001 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_sxr002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr002 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_sxr003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr003 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_sxr004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=sxr004 --attribute=sigma --value=5.00000E-03 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int001 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int002 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int003 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int004 --attribute=sigma --value=5.00000E+15 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int005 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int006 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int007 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int008 --attribute=sigma --value=1.00000E-10 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int009 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int010 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int011 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_int012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=int012 --attribute=sigma --value=5.52983E-07 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm001
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm001 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm002
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm002 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm003
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm003 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm004
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm004 --attribute=sigma --value=1.50000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm005
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm005 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm006
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm006 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm007
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm007 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm008
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm008 --attribute=sigma --value=1.0E16 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm009
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm009 --attribute=sigma --value=4 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm010
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm010 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm011
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm011 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_thm012
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=thm012 --attribute=sigma --value=4 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_prior_1
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_1 --attribute=sigma --value=1.00000E-02 --range=1.0E-10)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_prior_2
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_2 --attribute=sigma --value=1.00000E-02 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_prior_3
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_3 --attribute=sigma --value=1.50000E-01 --range=1.0)
add_test (NAME    reconstruction_task_seg_step_test_check_sigma_prior_4
          COMMAND ${TEST_SCRIPT_DIR}/check_signal_value.py --file=${CMAKE_CURRENT_SOURCE_DIR}/recout.input.test.v3fit --s_name=prior_4 --attribute=sigma --value=1.00000E-02 --range=1.0)