!******************************************************************************
!  File v3fit.f
!*******************************************************************************
!**********************************************************************
!**                                                                  **
!**   MAIN PROGRAM:  Equilibrium Reconstruction                      **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          v3fit uses diagnostic data to reconstruct a VMEC        **
!**	     equilibrium that matches the data                           **
!**                                                                  **
!**     CALLING ARGUMENTS: This is a stand alone code                **
!**          xv3fit arg1                                             **
!**          arg1: optional input namelist file name                 **
!**                                                                  **
!**     Author: James D. Hanson, with help from:                     **
!**             Steve Hirshman                                       **
!**             Ed Lazarus                                           **
!**             Lang Lao                                             **
!**             Steve Knowlton                                       **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**********************************************************************
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!  This file uses the following modules:
!    stel_kinds
!       located in LIBSTELL/Sources/Modules
!
!    stel_constants
!       located in LIBSTELL/Sources/Modules
!
!    ezcdf
!       located in LIBSTELL/Sources/Ezcdf
!
!    read_wout_mod
!       located in LIBSTELL/Sources/Modules
!
!    safe_open_mod 
!       located in LIBSTELL/Sources/Modules
!
!    v3f_global
!       located in V3FITA/Sources
!       storage for global variables
!
!    recon_param_T
!       located in V3FITA/Sources
!       Type declarations for reconstruction parameters
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  Executing 'xv3fit -h' will printout a help message
!
!    INPUT FILES
!
!  There are many input files read in by xv3fit:
! 1) A namelist input file
!      The first argument on the execute line is the name of this file
!      If there is no argument on the execute line, the default name of the
!         namelist input file is v3fit.in
!      For example, 'xv3fita my_nl_file' reads namelist input from a file named
!         my_nl_file
!
! 2) A LIST file
!
! 3) Many mdsig.nc files o one for each diagnostic
!
! 4) A VMEC namelist input file
!
!   OUTPUT FILES
!      
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!  The module v3_global contains the declarations for the globally accessible
!   variables.
!
!*******************************************************************************
!  CODE V3FITA
!    
! SECTION I.	Main Program
! SECTION II.	Initialization Subroutines
!    initialize_read   - called from Main
!       init_command_line       - parse the command line
!       init_main_nli           - read in the main namelist input file
!       init_signal             - read in the signal information
!                                 (9/04 - reads in mdsig files)
!       init_equilibrium        - read in the initial mhd equilibrium state
!    initialize_compute - called from Main
! SECTION III.	TASK SUBROUTINES
!
!
! SECTION IV. Cleanup Subroutines 
!    cleanup           - called from Main (NYI 07-01-2005)
!    cleanup_vmec      - deallocate vmec arrays and write out files
!        (NYI 07-01-2005)
!    help               - (NYI 07-01-2005)
! SECTION V.	SUBROUTINES FOR TESTNG
! SECTION VI.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
      PROGRAM v3fit
      USE v3fit_context

      IMPLICIT NONE

!  local variables
      TYPE (v3fit_context_class), POINTER :: context => null()
      REAL (rprec)                        :: time_start, time_end

!  Start of executable code
      CALL second0(time_start)

!  Create a context and read in the v3fit namelist input file.
      context => v3fit_context_construct(commandline_parser_construct())

!  Configure parallel computation.
      CALL config_parallelism(context)

!  Run unit tests if the -test flag is set.
      IF (commandline_parser_is_flag_set(context%cl_parser,                    &
     &                                   '-test')) THEN
         my_task = 'unit_test'
      ELSE
         CALL v3fit_input_read_namelist(                                       &
     &           commandline_parser_get_string(context%cl_parser,              &
     &                                         '-file'))
      END IF

      WRITE (context%recout_iou,*) '  my_task is ', my_task
      WRITE (context%runlog_iou,*) '  my_task is ', my_task
      SELECT CASE (my_task)

         CASE ('equilibrium')
            CALL task_equilibrium(context)

         CASE ('v3post', 'vmec_v3post')
            CALL task_v3post(context)

         CASE ('reconstruct', 'reconstruct_a1')
            CALL task_reconstruct(context)

         CASE ('unit_test')
            CALL task_unit_tests

         CASE DEFAULT
            CALL err_fatal('Unkown task type ' // TRIM(my_task))

      END SELECT

!  Clean up memory
      CALL second0(time_end)

      CALL v3fit_context_write(context)

      WRITE (context%runlog_iou,'(/,a,f10.2)')                                 &
     &   '  *** V3FIT run time (s): ', time_end - time_start
      WRITE (context%recout_iou,'(/,a,f10.2)')                                 &
     &   '  *** V3FIT run time (s): ', time_end - time_start
      WRITE (*,'(/,a,f10.2)')                                                  &
     &   '  *** V3FIT run time (s): ', time_end - time_start


      CALL v3fit_context_destruct(context)

      END PROGRAM

!*******************************************************************************
!  TASK SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Solves the equilibrium.
!>
!>  Task to just solve the equilibrium. This task is choosen by setting
!>  @ref v3fit_input::my_task to 'equilibrium'.
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE task_equilibrium(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: num_steps = 1
      LOGICAL                                   :: converged

!  Start of executable code
      WRITE (*,*) ' *** In task equilibrium'
      WRITE (context%runlog_iou,*) ' *** In task equilibrium'

!  Initialize equilibrium
      CALL init_equilibrium(context)

!  Solve equilibrium
      WRITE (*,*) ' *** Starting initial equilibrium convergence'
      WRITE (context%runlog_iou,*) ' *** Starting initial ' //                 &
     &                             'equilibrium convergence'
      converged = equilibrium_converge(context%model%equilibrium,              &
     &                                 num_steps, context%runlog_iou)

      CALL v3fit_context_init_data(context)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Solves the equilibrium and calculates the modeled signals.
!>
!>  Task to just compute modeled signals from an equilibrium. This task is
!>  choosen by setting @ref v3fit_input::my_task to 'v3post'.
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE task_v3post(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: num_steps = 1
      LOGICAL                                   :: converged
      INTEGER                                   :: i, last_para_signal
      REAL (rprec)                              :: g2

!  Start of executable code
      WRITE (*,*) ' *** In task v3post'
      WRITE (context%runlog_iou,*) ' *** In task v3post'

!  Initialize v3post
      CALL init_equilibrium(context)
      CALL init_signals(context)

!  Solve equilibrium
      WRITE (*,*) ' *** Starting initial equilibrium convergence'
      WRITE (context%runlog_iou,*) ' *** Starting initial ' //                 &
     &                             'equilibrium convergence'
      converged = equilibrium_converge(context%model%equilibrium,              &
     &                                 num_steps, context%runlog_iou)

!  Calculate and display the signals.
      WRITE (*,*) ' *** Calculating modeled signals'
      WRITE (context%runlog_iou,*) ' *** Calculating modeled signals'

      g2 = 0.0

!  The combination signals cannot be computed in parallel since the depend on
!  other signals previously computed. If combinations are used, the signals need
!  to be split into two loops.
      IF (context%combination_index .eq. -1) THEN
         last_para_signal = SIZE(context%signals)
      ELSE
         last_para_signal = context%combination_index - 1
      END IF

!$OMP PARALLEL DO
!$OMP& SCHEDULE(DYNAMIC)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i)
!$OMP& REDUCTION(+:g2)
      DO i = 1, last_para_signal
         g2 = g2 + signal_get_g2(context%signals(i)%p,                         &
     &                           context%model, .false.)
      END DO
!$OMP END PARALLEL DO

      DO i = last_para_signal + 1, SIZE(context%signals)
         g2 = g2 + signal_get_g2(context%signals(i)%p,                         &
     &                           context%model, .false.)
      END DO

      WRITE (*,1000) g2
      WRITE (context%runlog_iou,1000) g2

      CALL v3fit_context_init_data(context)

1000  FORMAT('g^2 = ',es12.5)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reconstructs the equilibrium.
!>
!>  Task to perform an equilibrium reconstruction. This task is choosen by
!>  setting @ref v3fit_input::my_task to 'reconstruct'.
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE task_reconstruct(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: num_steps = 1
      INTEGER                                   :: i
      LOGICAL                                   :: eq_converged

!  Start of executable code
      WRITE (*,*) ' *** In task reconstruct'

!  Initialize the reconstruction.
      CALL init_equilibrium(context)
      CALL init_signals(context)
      CALL init_parameters(context)
      CALL init_reconstruction(context)

!  Check if the number of observed signals matches the number of signals
!  created. Otherwise the reconstruction can be put in an infinite loop when
!  performing the SVD.
      CALL assert_eq(SIZE(context%signals), n_sdata_o,                         &
     &               'task_reconstruct: number of observed signals' //         &
     &               ' does not match the number of signals created')

!  Solve equilibrium
      WRITE (*,*) ' *** Starting initial equilibrium convergence'
      WRITE (context%runlog_iou,*) ' *** Starting initial ' //                 &
     &                             'equilibrium convergence'

      eq_converged = reconstruction_eval_e(context%recon,                      &
     &                                     context%signals,                    &
     &                                     context%model, num_steps,           &
     &                                     context%runlog_iou)
      CALL assert(eq_converged,                                                &
     &            'task_reconstruct no inital convergence')
      CALL reconstruction_eval_f(context%recon, context%derived_params,        &
     &                           context%model)

      WRITE (*,*)
      WRITE (*,1002) 0
      WRITE (*,1003) reconstruction_get_g2(context%recon)

      WRITE (context%runlog_iou,*)
      WRITE (context%runlog_iou,1002) 0
      WRITE (context%runlog_iou,1003)                                          &
     &   reconstruction_get_g2(context%recon)

      CALL v3fit_context_init_data(context)

      DO i = 1, SIZE(context%recon%e, 1) - 1
         WRITE (*,*)
         WRITE (*,1002) i

         WRITE (context%runlog_iou,*)
         WRITE (context%runlog_iou,1002) i

         IF (reconstruction_step(context%recon, context%signals,               &
     &                           context%derived_params,
     &                           context%model, context%params,                &
     &                           num_steps, context%runlog_iou)) THEN

!  Step was successful. Post process the result.
            CALL reconstruction_eval_sem(context%recon,                        &
     &                                   context%params,                       &
     &                                   context%signals,                      &
     &                                   context%derived_params)
            CALL v3fit_context_write_step_data(context, .false.)

!  Test for convergence.
            IF (reconstruction_get_dg2(context%recon) .lt.                     &
     &          context%recon_stop) THEN
               WRITE (*,*)
               WRITE (*,*) ' *** Equilibrium reconstructed'
               WRITE (*,1000) context%recon_stop,                              &
     &                        reconstruction_get_dg2(context%recon)

               WRITE (context%runlog_iou,*)
               WRITE (context%runlog_iou,*) ' *** Equilibrium ' //             &
     &                                      'reconstructed'
               WRITE (context%runlog_iou,1000) context%recon_stop,             &
     &                        reconstruction_get_dg2(context%recon)
               EXIT
            END IF

            IF (i .eq. SIZE(context%recon%e, 1) - 1) THEN
              WRITE (*,*)
              WRITE (*,1001) i

              WRITE (context%runlog_iou,*)
              WRITE (context%runlog_iou,1001) i
            END IF
         ELSE
            EXIT
         END IF
      END DO

1000  FORMAT('dg2_stop = ',es12.5,' change in g^2 = ',es12.5)
1001  FORMAT('nrstep 'i4,' completed. May be reconstructed.')
1002  FORMAT('  *** Reconstruction step ',i4)
1003  FORMAT('Initial g^2 = ',es12.5)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Run unit tests
!>
!>  This task is choosen by setting @ref v3fit_input::my_task to 'unit_test'. Or
!>  by using the -test command line flag.
!-------------------------------------------------------------------------------
      SUBROUTINE task_unit_tests
      USE line_segment
      USE coordinate_utilities
      USE functions
      USE integration_path

      IMPLICIT NONE

!  local variables
      LOGICAL :: tests_failed

!  Start of executable code
      WRITE (*,*) ' *** In task unit_test'
      tests_failed = .false.

! Line segment tests
      IF (line_seg_test()) THEN
         WRITE(*,*) "Line Segment Test Passed"
      ELSE
         tests_failed = .true.
      END IF

! Coordinate tests
      IF (cood_utils_test()) THEN
         WRITE(*,*) "Coord Utilites Test Passed"
      ELSE
         tests_failed = .true.
      END IF

! Profile function tests
      IF (function_test()) THEN
         WRITE(*,*) "Functions Test Passed"
      ELSE
         tests_failed = .true.
      END IF

! Integration path tests
      IF (path_test()) THEN
         WRITE(*,*) "Integration Path Test Passed"
      ELSE
         tests_failed = .true.
      END IF

!  The CTest unit test frame work chooses a pass fail based on the exit code of
!  the program.
      IF (tests_failed) THEN
         CALL EXIT(1)
      END IF

      END SUBROUTINE

!*******************************************************************************
!  INITIALIZE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Initialize an equilibrium.
!>
!>  This constructs and initializes an equilibrium. Currently only
!>  @ref vmec_equilibium and @ref vacuum_equilibium are supported.
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE init_equilibrium(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      LOGICAL                                   :: force_solve

!  Start of executable code
      WRITE (*,*) ' *** Constructing equilibrium model'
      WRITE (context%runlog_iou,*) ' *** Constructing equilibrium model'

      force_solve =                                                            &
     &   commandline_parser_is_flag_set(context%cl_parser, '-force')

      IF (vmec_nli_filename .ne. '') THEN
         CALL init_vmec_equilibrium(context, force_solve)
         RETURN
      END IF

      IF (vacuum_nli_filename .ne. '') THEN
         CALL init_vacuum_equilibrium(context, force_solve)
         RETURN
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize a vmec equilibrium.
!>
!>  This constructs and initializes a @ref vmec_equilibrium.
!>
!>  @param[inout] context     An instance of a @ref v3fit_context object.
!>  @param[in]    force_solve If true, forces the equilbirum to resolve every
!>                            time.
!-------------------------------------------------------------------------------
      SUBROUTINE init_vmec_equilibrium(context, force_solve)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout)      :: context
      LOGICAL, INTENT(in)                            :: force_solve

!  local variables
      TYPE (pprofile), POINTER                       :: ne => null()
      TYPE (pprofile), POINTER                       :: te => null()
      TYPE (pprofile_pointer), DIMENSION(:), POINTER :: sxr => null()
      INTEGER                                        :: i
      LOGICAL                                        :: not_converged

!  Start of executable code
      ALLOCATE(ne, te)
      CALL pprofile_construct(ne, TRIM(pp_ne_ptype), pp_ne_b,                  &
     &                        pp_ne_as, pp_ne_af)
      CALL pprofile_construct(te, TRIM(pp_te_ptype), pp_te_b,                  &
     &                        pp_te_as, pp_te_af)

      ALLOCATE(sxr(num_sxrem_p))
      DO i = 1, num_sxrem_p
         ALLOCATE(sxr(i)%p)
         CALL pprofile_construct(sxr(i)%p, TRIM(pp_sxrem_ptype_a(i)),          &
     &                           pp_sxrem_b_a(i,:), pp_sxrem_as_a(i,:),        &
     &                           pp_sxrem_af_a(i,:))
      END DO

!  If the woutfile in the v3fit namelist input file is specifed, assume that the
!  equilibrium is converged.
      not_converged = TRIM(vmec_wout_input) .eq. ''

      context%model => model_construct(model_ne_type,                          &
     &   model_sxrem_type_a(1:num_sxrem_p), model_te_type,                     &
     &   ne_pp_unit, ne_min, e_pressure_fraction,                              &
     &   equilibrium_construct(vmec_construct(vmec_nli_filename,               &
     &                                        vmec_wout_input,                 &
     &                                        ne, te, sxr, phi_offset,         &
     &                                        z_offset, pol_rad_ratio,         &
     &                                        context%runlog_iou),             &
     &                         force_solve, not_converged))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize a vacuum equilibrium.
!>
!>  This constructs and initializes a @ref vacuum_equilibrium.
!>
!>  @param[inout] context     An instance of a @ref v3fit_context object.
!>  @param[in]    force_solve If true, forces the equilbirum to resolve every
!>                            time.
!-------------------------------------------------------------------------------
      SUBROUTINE init_vacuum_equilibrium(context, force_solve)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      LOGICAL, INTENT(in)                       :: force_solve

!  Start of executable code

      context%model => model_construct(model_ne_type,                          &
     &   model_sxrem_type_a(1:num_sxrem_p), model_te_type,                     &
     &   ne_pp_unit, ne_min, e_pressure_fraction,                              &
     &   equilibrium_construct(vacuum_construct(vacuum_nli_filename,           &
     &                                          context%runlog_iou),           &
     &                         force_solve))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs the @ref signal array.
!>
!>  This constructs all @ref v3fit_context::signals.
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE init_signals(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: signals_created = 0

!  Start of executable code
      IF (na_s_desc .eq. 0) THEN
         CALL err_warn('init_signals: no signals will be created. ' //         &
     &                 'na_s_desc = 0')
         RETURN
      END IF

      WRITE (*,*) ' *** Constructing signals'
      WRITE (context%runlog_iou,*) ' *** Constructing signals'

      ALLOCATE(context%signals(na_s_desc))

!  If a mdsig filename is given, construct the diagnostics contained inside it.
      IF (mdsig_list_filename .ne. '') THEN
         CALL init_magnetic_signals(context, signals_created)
      END IF

!  If a sxrem filename is given, construct the diagnostics contained inside it.
      IF (sxrch_dot_filename .ne. '') THEN
         CALL init_sxrem_signals(context, signals_created)
      END IF

!  If a intpol filename is given, construct the diagnostics contained inside it.
      IF (ipch_dot_filename .ne. '') THEN
         CALL init_intpol_signals(context, signals_created)
      END IF

!  If a thscte filename is given, construct the diagnostics contained inside it.
      IF (thscte_dot_filename .ne. '') THEN
         CALL init_thomson_signals(context, signals_created)
      END IF

!  If the radial extcur z position isn't negative, construct an extcurz signal.
      IF (extcurz_s0 .ge. 0.0) THEN
         CALL init_extcurz_signals(context, signals_created)
      END IF

      CALL init_geometric_signals(context, signals_created)

      CALL init_prior_signals(context, signals_created)

      CALL init_combination_signals(context, signals_created)

      CALL v3fit_context_resize(context)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref magnetic signals.
!>
!>  This constructs all @ref magnetic signals from the specified
!>  @ref v3fit_input::mdsig_list_filename.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_magnetic_signals(context, signals_created)
      USE v3fit_context
      USE safe_open_mod
      USE v3_utilities
      USE ezcdf

      IMPLICIT NONE
!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  local variables
      INTEGER                                   :: mdsig_list_iou
      INTEGER                                   :: mdsig_iou
      INTEGER                                   :: status
      INTEGER                                   :: mdsig_index
      CHARACTER (len=data_filename_length)      :: mdsig_path
      CHARACTER (len=data_filename_length)      :: mdsig_filename

      TYPE (magnetic_class), POINTER :: magnetic_object => null()
      LOGICAL                        :: first_run = .true.

!  Start of executable code
      status = 0
      mdsig_list_iou = 0

      WRITE (*,*) ' *** Constructing magnetic signals from ' //                &
     &            TRIM(mdsig_list_filename)
      WRITE (context%runlog_iou,*) ' *** Constructing magnetic ' //            &
     &                             'signals from ' //                          &
     &                             TRIM(mdsig_list_filename)

      CALL safe_open(mdsig_list_iou, status, TRIM(mdsig_list_filename),        &
     &               'old', 'formatted')
      CALL assert_eq(0, status, 'init_magnetic_signals: ' //                   &
     &               'failed to open ' // TRIM(mdsig_list_filename))

!  Get the directory of the mdsig path.
      mdsig_path =                                                             &
     &   mdsig_list_filename(1:INDEX(TRIM(mdsig_list_filename),                &
     &                       '/', back = .true.))

!  Loop over the file in the list file
      DO
         READ(mdsig_list_iou, *, iostat=status)                                &
     &      mdsig_index, mdsig_filename

         IF (status .ne. 0) RETURN

!  Check if there are two many signals.
         IF (signals_created + 1 .gt. v3fit_max_diagnostics) THEN
            CALL err_fatal('init_magnetic_signals: created signals' //         &
     &                     ' exceeds v3fit_max_diagnostics')
         END IF

!  Open the netcdf file for the mdsig.
         CALL cdf_open(mdsig_iou, TRIM(mdsig_path) //                          &
     &                 TRIM(mdsig_filename), 'r', status)
         CALL assert_eq(0, status, 'init_magnetic_signals: ' //                &
     &                  'failed to open ', TRIM(mdsig_path) //                 &
     &                  TRIM(mdsig_filename))

!  Construct a magnetic diagnostic signal using from the mdsig netcdf file.
         signals_created = signals_created + 1
         magnetic_object => magnetic_construct(mdsig_iou)
         context%signals(signals_created)%p => signal_construct(               &
     &      diagnostic_construct(magnetic_object),                             &
     &      mdsig_iou, sdo_data_a(signals_created),                            &
     &      sdo_sigma_a(signals_created), sdo_weight_a(signals_created))

!  If a magnetic signal is created, the magnetic cache must be used.
         IF (first_run) THEN
            CALL equilibrium_set_magnetic_cache(                               &
     &              context%model%equilibrium,                                 &
     &              magnetic_object%response_function%kp,                      &
     &              magnetic_object%response_function%kp_store,                &
     &              magnetic_object%response_function%n_field_periods)
            first_run = .false.
         END IF

         CALL cdf_close(mdsig_iou)

!  At lease one magnetic signal was made. Mark the index of it in the context.
!  This should only be run once.
         IF (context%magnetic_index .eq. -1) THEN
            context%magnetic_index = signals_created
         END IF

      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref sxrem signals.
!>
!>  This constructs all @ref sxrem signals from the specified
!>  @ref v3fit_input::sxrch_dot_filename.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_sxrem_signals(context, signals_created)
      USE v3fit_context
      USE sxrch_dot

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  Start of executable code
      WRITE (*,*) ' *** Constructing sxrem signals from ' //                   &
     &            TRIM(sxrch_dot_filename)
      WRITE (context%runlog_iou,*) ' *** Constructing sxrem ' //               &
     &                             'signals from ' //                          &
     &                             TRIM(sxrch_dot_filename)

!  The setting of the sxrem_index for sxrem happens inside sxrch_dot_parse_chord
      CALL sxrch_dot_read(sxrch_dot_filename, context%signals,                 &
     &                    signals_created, sdo_data_a, sdo_sigma_a,            &
     &                    sdo_weight_a, context%sxrem_index)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref intpol signals.
!>
!>  This constructs all @ref intpol signals from the specified
!>  @ref v3fit_input::ipch_dot_filenam.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_intpol_signals(context, signals_created)
      USE v3fit_context
      USE ipch_dot

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  Start of executable code
      WRITE (*,*) ' *** Constructing intpol signals from ' //                  &
     &            TRIM(ipch_dot_filename)
      WRITE (context%runlog_iou,*) ' *** Constructing intpol ' //              &
     &                             'signals from ' //                          &
     &                             TRIM(ipch_dot_filename)

!  The setting of the intpol_index for sxrem happens inside ipch_dot_parse_chord
      CALL ipch_dot_read(ipch_dot_filename, context%signals,                   &
     &                   signals_created, sdo_data_a, sdo_sigma_a,             &
     &                   sdo_weight_a, context%intpol_index)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref thomson signals.
!>
!>  This constructs all @ref thomson signals from the specified
!>  @ref v3fit_input::thscte_dot_filename.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_thomson_signals(context, signals_created)
      USE v3fit_context
      USE thscte_dot

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  Start of executable code
      WRITE (*,*) ' *** Constructing thomson signals from ' //                 &
     &            TRIM(thscte_dot_filename)
      WRITE (context%runlog_iou,*) ' *** Constructing thomson ' //             &
     &                             'signals from ' //                          &
     &                             TRIM(thscte_dot_filename)

!  The setting of the thomson_index for sxrem happens inside
!  thscte_dot_parse_chord
      CALL thscte_dot_read(thscte_dot_filename, context%signals,               &
     &                     signals_created, sdo_data_a, sdo_sigma_a,           &
     &                     sdo_weight_a, context%thomson_index)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs the @ref extcurz signal.
!>
!>  This constructs the extcurz signal from @ref v3fit_input::extcurz_s0 and
!>  @ref v3fit_input::extcurz_u0.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_extcurz_signals(context, signals_created)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  Start of executable code
      WRITE (*,*) ' *** Constructing extcurz signals'
      WRITE (context%runlog_iou,*) ' *** Constructing extcurz signals'

!  Check if there are two many signals.
      IF (signals_created + 1 .gt. v3fit_max_diagnostics) THEN
         CALL err_fatal('init_extcurz_signals: created signals' //             &
     &                  ' exceeds v3fit_max_diagnostics')
      END IF

      context%signals(signals_created + 1)%p => signal_construct(              &
     &   diagnostic_construct(extcurz_construct(extcurz_s0,                    &
     &                                          extcurz_u0)),                  &
     &   'extcurz', 'EXTernal CURrent along Z', 'A',                           &
     &   sdo_data_a(signals_created + 1),                                      &
     &   sdo_sigma_a(signals_created + 1),                                     &
     &   sdo_weight_a(signals_created + 1))

      signals_created = signals_created + 1

!  At lease one extcurz signal was made. Mark the index of it in the context.
!  This should only be run once.
      IF (context%extcurz_index .eq. -1) THEN
         context%extcurz_index = signals_created
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref geometric signals.
!>
!>  For now only @ref limiter signals are constructed.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_geometric_signals(context, signals_created)
      USE v3fit_context
      USE stel_constants, only: degree

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  local variables
      INTEGER                                   :: i
      CHARACTER (len=data_short_name_length)    :: s_name = ''
      CHARACTER (len=data_name_length)          :: l_name = ''
      TYPE (limiter_iso), POINTER               :: iso_funct => null()
      INTEGER                                   :: num_phi

!  Start of executable code
      WRITE (*,*) ' *** Constructing geometric signals'
      WRITE (context%runlog_iou,*) ' *** Constructing geometric signals'

!  For now, only limiter functions can be created.
      DO i = 1, n_lif
!  Check if there are two many signals
         IF (signals_created + 1 .gt. v3fit_max_signals) THEN
            CALL err_fatal('init_geometric_signals: created signals' //        &
     &                     ' exceeds v3fit_max_signals')
         END IF

!  Check if there are two many geometric signals.
         IF (i .gt. v3fit_max_geometrics) THEN
            CALL err_fatal('init_geometric_signals: created signals' //        &
     &                     ' exceeds v3fit_max_geometrics')
         END IF

         s_name = TRIM(signal_make_short_name('edge_lim_', i))
         WRITE (l_name,'(a,i4)') 'edge limiter iso function ', i

         num_phi = MAX(MIN(v3fit_max_lif_size, n_phi_lif(i)), 1)
!  Do not deallocate the iso_funct. It memory is taken over by the limiter_class
!  and freed during it's destruction.
         ALLOCATE(iso_funct)
!  The conversion from degrees to radians is performed in v3fit_input
         CALL limiter_iso_construct(iso_funct, lif_arz(i,:,:),                 &
     &           lif_rc(i), lif_zc(i), 20,                                     &
     &           lif_phi_degree(i,1:num_phi)*degree)

!  Check if the sdo_sigma is zero and replace it lif_sigma instead.
         IF (sdo_sigma_a(signals_created + 1) .eq. 0.0) THEN
            sdo_sigma_a(signals_created + 1) = lif_sigma(i)
         END IF

         context%signals(signals_created + 1)%p => signal_construct(           &
     &      geometric_construct(                                               &
     &         limiter_construct(iso_funct, lif_on_edge(i))),                  &
     &      s_name, l_name, 'm', sdo_data_a(signals_created + 1),              &
     &      sdo_sigma_a(signals_created + 1),                                  &
     &      sdo_weight_a(signals_created + 1))

!  Set the iso_funct to null. This sets the local pointer of iso_funct to null.
!  The pointer allocated above is now owned by the limiter_class instance.
         iso_funct => null()

         signals_created = signals_created + 1

!  At lease one limiter signal was made. Mark the index of it in the context.
!  This should only be run once.
         IF (context%limiter_index .eq. -1) THEN
            context%limiter_index = signals_created
         END IF

      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref prior signals.
!>
!>  For now only @ref prior_guassian signals are constructed.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_prior_signals(context, signals_created)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  local variables
      INTEGER                                   :: i

!  Start of executable code
      WRITE (*,*) ' *** Constructing prior signals'
      WRITE (context%runlog_iou,*) ' *** Constructing prior signals'

!  For now, only limiter functions can be created.
      DO i = 1, n_prior
!  Check if there are two many signals
         IF (signals_created + 1 .gt. v3fit_max_signals) THEN
            CALL err_fatal('init_prior_signals: created signals' //            &
     &                     ' exceeds v3fit_max_signals')
         END IF

!  Check if there are two many geometric signals.
         IF (i .gt. v3fit_max_priors) THEN
            CALL err_fatal('init_prior_signals: created signals' //            &
     &                     ' exceeds v3fit_max_priors')
         END IF

         context%signals(signals_created + 1)%p => signal_construct(           &
     &      prior_construct(                                                   &
     &         prior_gaussian_construct(context%model,                         &
     &                                  prior_param_name(i),                   &
     &                                  prior_indices(i,:))),                  &
     &      TRIM(signal_make_short_name('prior_', i)), prior_name(i),          &
     &      prior_units(i), sdo_data_a(signals_created + 1),                   &
     &      sdo_sigma_a(signals_created + 1),                                  &
     &      sdo_weight_a(signals_created + 1))

         signals_created = signals_created + 1

!  At lease one prior signal was made. Mark the index of it in the context.
!  This should only be run once.
         IF (context%prior_guassian_index .eq. -1) THEN
            context%prior_guassian_index = signals_created
         END IF

      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs @ref signal::combination_class signals.
!>
!>  @param[inout] context         An instance of a @ref v3fit_context object.
!>  @param[inout] signals_created Count of the current number of signals
!>                                constructed.
!-------------------------------------------------------------------------------
      SUBROUTINE init_combination_signals(context, signals_created)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context
      INTEGER, INTENT(inout)                    :: signals_created

!  local variables
      INTEGER                                   :: i, j
      TYPE (combination_class), POINTER         :: combination_signal

!  Start of executable code
      WRITE (*,*) ' *** Constructing combination signals'
      WRITE (context%runlog_iou,*) ' *** Constructing ' //                     &
     &                             'combination signals'

      DO i = 1, n_coosig
!  Check if there are two many signals
         IF (signals_created + 1 .gt. v3fit_max_signals) THEN
            CALL err_fatal('init_combination_signals: created ' //             &
     &                     'signals exceeds v3fit_max_signals')
         END IF

!  Check if there are two many geometric signals.
         IF (i .gt. v3fit_max_combinations) THEN
            CALL err_fatal('init_combination_signals: created ' //             &
     &                     'signals exceeds v3fit_max_combinations')
         END IF

!  Search for the last index. It would be better if there was a 
         combination_signal => combination_construct(n_sig_coosig(i),          &
     &                                               coosig_type(i))
         DO j = 1, n_sig_coosig(i)
            IF (coosig_indices(i,j) .gt. signals_created) THEN
               CALL err_fatal('init_combination_signals: cannot ' //           &
     &                        'create combination signal with ' //             &
     &                        'uninitalized signal')
            END IF

            CALL combination_set_signal(combination_signal,                    &
     &              context%signals(coosig_indices(i,j))%p,                    &
     &              coosig_coeff(i,j), j)
         END DO

         context%signals(signals_created + 1)%p => signal_construct(           &
     &      combination_signal, coosig_name(i), coosig_name(i),                &
     &      coosig_units(i), sdo_data_a(signals_created + 1),                  &
     &      sdo_sigma_a(signals_created + 1),                                  &
     &      sdo_weight_a(signals_created + 1))

         signals_created = signals_created + 1

!  At lease one combination signal was made. Mark the index of it in the
!  context. This should only be run once.
      IF (context%combination_index .eq. -1) THEN
         context%combination_index = signals_created
      END IF

      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Constructs the @ref parameter array.
!>
!>  This constructs all reconstruction @ref v3fit_context::params and derived
!>  @ref v3fit_context::params .
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE init_parameters(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: i

!  Start of executable code
      WRITE (*,*) ' *** Constructing reconstruction parameters'
      WRITE (context%runlog_iou,*) ' *** Constructing ' //                     &
     &                             'reconstruction parameters'

!  Initialize the reconstruction parameters.
      ALLOCATE(context%params(n_rp))
      DO i = 1, n_rp
         context%params(i)%p => param_construct(context%model,                 &
     &                                          rp_type(i),                    &
     &                                          (/ rp_index(i),                &
     &                                             rp_index2(i) /),            &
     &                                          rp_vrnc(i),                    &
     &                                          rp_range_type(i,:),            &
     &                                          rp_range_index(i,:,:),         &
     &                                          rp_range_value(i,:),           &
     &                                          SIZE(context%signals),         &
     &                                          n_rp)
      END DO

      ALLOCATE(context%derived_params(n_dp))
      DO i = 1, n_dp
         context%derived_params(i)%p => param_construct(context%model,         &
     &                                                  dp_type(i),            &
     &                                                  dp_index(i,:),         &
     &                                                  n_dp)

      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize Reconstruction
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
      SUBROUTINE init_reconstruction(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: last_para_signal

!  Start of executable code
!  The combination signals cannot be computed in parallel since the depend on
!  other signals previously computed. If combinations are used, the signals need
!  to be split into two loops.
      IF (context%combination_index .eq. -1) THEN
         last_para_signal = SIZE(context%signals)
      ELSE
         last_para_signal = context%combination_index - 1
      END IF

      context%recon_stop = dg2_stop
      context%recon => reconstruction_construct(nrstep,                        &
     &                    SIZE(context%signals),                               &
     &                    SIZE(context%derived_params),                        &
     &                    SIZE(context%params),                                &
     &                    step_type, astep_max, cut_svd, cut_eff,              &
     &                    cut_marg_eff, cut_delta_a, cut_dg2,                  &
     &                    last_para_signal)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Configure Parallelism
!>
!>  @param[inout] context An instance of a @ref v3fit_context object.
!-------------------------------------------------------------------------------
!  When using openmp the !$ allows conditional computation based on if openmp is
!  used or not.
      SUBROUTINE config_parallelism(context)
      USE v3fit_context

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: context

!  local variables
      INTEGER                                   :: num_threads

!  Forward declare the interface for OMP_GET_MAX_THREADS. This explicitly
!  defined the type for OMP_GET_MAX_THREADS to avoid a compiler error.
!$    INTERFACE
!$       FUNCTION OMP_GET_MAX_THREADS()
!$       INTEGER :: OMP_GET_MAX_THREADS
!$       END FUNCTION
!$    END INTERFACE

!  Start of executable code
      WRITE (*,*) ' *** Configuring parallelism.'

!  Configure the number of threeds to use.
      num_threads = commandline_parser_get_integer(context%cl_parser,          &
     &                                             '-para', 1)

      IF (num_threads .gt. 0) THEN
!$       CALL OMP_SET_NUM_THREADS(num_threads)
!$OMP PARALLEL
!$       num_threads = OMP_GET_MAX_THREADS()
!$OMP END PARALLEL
         WRITE (*,1000) num_threads
      END IF

1000  FORMAT('Using ',i4,' threads.')

      END SUBROUTINE