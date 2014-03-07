#  Run the equilibrium task. This just tests if the equilibrium executes
#  correctly. Subsequent tests will examine the results of the equilibrium for
#  correct results.
add_test (NAME              equilibrium_task_test
          COMMAND           ${EXECUTABLE_OUTPUT_PATH}/xv3fit input.test.v3fit)
