enable_testing ()

#  The python interpreter may not be in the same location on all systems.
find_package (PythonInterp)

set (TEST_SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities)

#  Copy common files to the equilibrium  test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_task_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_task_test)

#  Copy common files to the equilibrium reset test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_reset_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_reset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_reset_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_reset_test)

#  Copy common files to the equilibrium reset free boundary test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_reset_free_boundary_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_reset_free_boundary_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_reset_free_boundary_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_reset_free_boundary_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/equilibrium_reset_free_boundary_test/coils.test
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/equilibrium_reset_free_boundary_test)

#  Copy common files to the v3post test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_task_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_task_test)

#  Copy common files to the reconstruction test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/magnetic.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_test/v3rfun.in
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_test)

#  Copy common files to the reconstruction seg step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_seg_step_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_seg_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_seg_step_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_seg_step_test)

#  Copy common files to the reconstruction sl step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_sl_step_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_sl_step_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/reconstruction_task_sl_step_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/reconstruction_task_sl_step_test)

#  Copy common files to the parameter constraints test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/param_constraints_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/param_constraints_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/param_constraints_test)

#  Copy common files to the te pressure test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/te_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/te_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/te_pressure_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/te_pressure_test)

#  Copy common files to the ne pressure test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/ne_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/ne_pressure_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/ne_pressure_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/ne_pressure_test)

#  Copy common files to the v3post lasym test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_lasym_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_lasym_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/v3post_lasym_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/v3post_lasym_test)

#  Copy common files to the reconstruction test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/z_offset_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/z_offset_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/z_offset_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/z_offset_test)

#  Copy common files to the diagnostic rotation test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/thomson_rot.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/v3rfun.in
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/coils.null
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/input.unrot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/intpol_rot.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/diagnostic_rot.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/diagnostic_rotation_test/sxrem_rot.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/diagnostic_rotation_test)

#  Copy common files to the reconstruction seg step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/multigrid_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/multigrid_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/multigrid_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/multigrid_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/multigrid_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/multigrid_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/multigrid_test)

#  Copy common files to the aphi test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/aphi_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/aphi_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/aphi_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/aphi_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/aphi_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/aphi_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/aphi_test)

#  Copy common files to the reconstruction seg step test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/sxrem.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/intpol.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/test_files/thomson.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/preconditioner_test/magnetic.dot
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/preconditioner_test/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/preconditioner_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/preconditioner_test/v3rfun.in
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/preconditioner_test)

#  Copy common files to the lasym boundary test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/lasym_boundary/input.test.v3fit
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/lasym_boundary)
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/lasym_boundary/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/lasym_boundary)

#  Copy common files to the niter  test directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tests/niter_test/input.test.vmec
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/niter_test)

add_subdirectory (tests)
