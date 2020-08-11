#  Run VMEC and SIESTA the compare the results.
add_test (NAME    siesta_init_vmec_test
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} ${EXECUTABLE_OUTPUT_PATH}/xvmec input.test.vmec)
add_test (NAME    siesta_init_siesta_test
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} ${EXECUTABLE_OUTPUT_PATH}/xsiesta)

add_test (NAME    siesta_init_curtor_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=curtor -tol=0.02)

add_test (NAME    siesta_init_chipf_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=chipf -min=0.0 -max=1,0 -tol=2.0E-6)
add_test (NAME    siesta_init_phipf_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=phipf -min=0.0 -max=1.0 -tol=1.0E-17)

add_test (NAME    siesta_init_r_1_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_init_r_2_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_init_r_3_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_init_r_4_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=r -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_init_drdu_1_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_init_drdu_2_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_init_drdu_3_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_init_drdu_4_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdu -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_init_drdv_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=drdv -min=0.0 -max=1.0 -tol=1.0E-20)

add_test (NAME    siesta_init_z_1_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_init_z_2_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_init_z_3_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_init_z_4_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=z -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_init_dzdu_1_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.0 -max=0.01 -tol=1.0E-2)
add_test (NAME    siesta_init_dzdu_2_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.01 -max=0.04 -tol=1.0E-3)
add_test (NAME    siesta_init_dzdu_3_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.04 -max=0.2 -tol=1.0E-4)
add_test (NAME    siesta_init_dzdu_4_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdu -min=0.2 -max=1.0 -tol=1.0E-5)

add_test (NAME    siesta_init_dzdv_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=dzdv -min=0.0 -max=1.0 -tol=1.0E-20)

add_test (NAME    siesta_init_jacobian_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jacobian -min=0.0 -max=1.0 -tol=1.0E-4)

add_test (NAME    siesta_init_jbsups_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.0E-20)
add_test (NAME    siesta_init_jbsupu_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=1.0E-4)
add_test (NAME    siesta_init_jbsupv_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupv -min=0.0 -max=1.0 -tol=2.0E-4)

add_test (NAME    siesta_init_vmec_divb_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=vmec_divb -min=0.0 -max=1.0 -tol=2.0E-5)
add_test (NAME    siesta_init_siesta_divb_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=siesta_divb -min=0.0 -max=1.0 -tol=1.0E-20)

add_test (NAME    siesta_init_siesta_pressure_test
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=pressure -min=0.0 -max=1.0 -tol=1.0)

add_test (NAME    siesta_init_bsups_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsups -min=0.0 -max=1.0 -tol=1.0E-20)
add_test (NAME    siesta_init_bsupu_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupu -min=0.005 -max=0.99 -tol=1.0E-3 -relative)
add_test (NAME    siesta_init_bsupv_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsupv -min=0.015 -max=0.99 -tol=2.0E-3 -relative)

add_test (NAME    siesta_init_bsubu_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubu -min=0.005 -max=1.0 -tol=2.0E-1 -relative)
add_test (NAME    siesta_init_bsubv_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=bsubv -min=0.005 -max=1.0 -tol=1.0E-5 -relative)

add_test (NAME    siesta_init_jksups_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsups -min=0.0 -max=1.0 -tol=1.0E-10)
add_test (NAME    siesta_init_jksupu_test
          COMMAND ${EXECUTABLE_OUTPUT_PATH}/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_restart_w7x.nc -test_name=jbsupu -min=0.0 -max=1.0 -tol=2.0 -relative)
