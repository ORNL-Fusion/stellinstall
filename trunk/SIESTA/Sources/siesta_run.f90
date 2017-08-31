!*******************************************************************************
!>  @file siesta_run.f90
!>  @brief Contains module @ref siesta_run.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref siesta_run_class. This module
!>  contains all the code necessary to interface setup and run SIESTA.
!*******************************************************************************

      MODULE siesta_run
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) siesta_run base class
!
!*******************************************************************************
      TYPE siesta_run_class
!>  Timer On
         REAL (rprec) :: skston
      END TYPE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct new @ref siesta_run_class object.
!>
!>  Allocates memory and initializes a @ref siesta_run_class object. Performs
!>  all the initialization needed to operate SIESTA.
!>
!>  @param[in] run_comm MPI Communicator to use.
!>  @param[in] init_mpi Instructs if MPI should be initialized.
!>  @returns A pointer to a constructed @ref siesta_run_class object.
!-------------------------------------------------------------------------------
      FUNCTION siesta_run_construct(run_comm, init_mpi, verbose)
      USE nscalingtools, ONLY: MyEnvVariables, PARSOLVER, SIESTA_COMM
      USE descriptor_mod, ONLY: LSCALAPACK, INHESSIAN, iam, nprocs,            &
                                nprow, npcol, icontxt, icontxt_1xp,            &
                                icontxt_px1, icontxt_global, isroot,           &
                                myrow, mycol
      USE blocktridiagonalsolver_s, ONLY: Initialize, GetRanks
      USE prof_mod, ONLY: profinit
      USE shared_data, ONLY: nprecon, lverbose
      USE perturbation, ONLY: init_data
      USE timer_mod, ONLY: init_data_time, sks_timers
      USE hessian, ONLY: iprec_flag, HESSPASS
      USE siesta_state, ONLY: lfirst
      USE siesta_error
      USE diagnostics_mod, ONLY: toroidal_flux0

      IMPLICIT NONE

!  Declare Arguments
      TYPE (siesta_run_class), POINTER :: siesta_run_construct
      INTEGER, INTENT(in)              :: run_comm
      LOGICAL, INTENT(in)              :: init_mpi
      LOGICAL, INTENT(in)              :: verbose

!  local variables.
      REAL (rprec)                     :: ton
      REAL (rprec)                     :: toff

!  Start of executable code
      ALLOCATE(siesta_run_construct)

      CALL siesta_error_clear_all

      lverbose = verbose
      lfirst = .true.
      toroidal_flux0=0

      SIESTA_COMM = run_comm
      iprec_flag = 0
      HESSPASS = 0

#if defined(MPI_OPT)
      LSCALAPACK = .true.
      INHESSIAN = .false.
#if defined(SKS)
      CALL MyEnvVariables
      IF (PARSOLVER) LSCALAPACK = .false.
      CALL sks_timers
      CALL Initialize(init_mpi, 0, 0)
      CALL GetRanks(iam, nprocs)
#endif
      IF (LSCALAPACK) THEN
         CALL profinit()

! setup blacs
         CALL blacs_pinfo(iam, nprocs)
         CALL blacs_setup(iam, nprocs)
         DO nprow = INT(SQRT(DBLE(nprocs))) + 1, 1, -1
            IF (nprow*npcol .eq. nprocs) THEN
               EXIT
            END IF
         END DO

         CALL blacs_get(0, 0, icontxt)
         CALL blacs_gridinit(icontxt, 'C', nprow, npcol)

         CALL blacs_get(0, 0, icontxt_1xp)
         CALL blacs_gridinit(icontxt_1xp, 'R', 1, nprow*npcol)

         CALL blacs_get(0, 0, icontxt_px1)
         CALL blacs_gridinit(icontxt_px1, 'C', nprow*npcol, 1)

         icontxt_global = icontxt_1xp

         CALL blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
         isroot = (myrow .eq. 0) .and. (mycol .eq. 0)
         IF (isroot .and. lverbose) THEN
            WRITE (*,1000) nprow, npcol
            WRITE (*,1001) icontxt_global, icontxt, icontxt_1xp, icontxt_px1
         END IF
      END IF
#else
      LSCALAPACK = .false.
      iam = 0
      nprocs = 1
#endif

      CALL second0(ton)
      siesta_run_construct%skston = ton

      CALL init_data(nprecon)
#if defined(SKS)
      CALL second0(toff)
      init_data_time = init_data_time + (toff - ton)
#endif

1000  FORMAT('nprow,npcol ',2(x,i5))
1001  FORMAT('icontxt_global,icontxt,icontxt_1xp,icontxt_px1 ',4(x,i5))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref siesta_run_class object.
!>
!>  Deallocates memory and uninitializes a @ref siesta_run_class object. This 
!>  also performs the shut down procedures for SIESTA.
!>
!>  @param[inout] this         A @ref siesta_run_class instance.
!>  @param[in]    finalize_mpi Instructs if to finalize the MPI call.
!>  @param[in]    close_wout   Instructs if to finalize the close the wout file.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_destruct(this, finalize_mpi, close_wout)
      USE timer_mod
      USE nscalingtools, ONLY: OUTPUT_TIMINGS, GetTimes, PARSOLVER,            &
                               PARFUNCTISL, PARGMRES, SKSDBG, TOFU,            &
                               finalizeremap
      USE descriptor_mod, ONLY: iam, icontxt, LSCALAPACK, myrow, mycol, nprocs
      USE shared_data, ONLY: nprecon, xc, gc, fsq_total1, mblk_size, niter,    &
                             col_scale
      USE hessian, ONLY: levmarq_param, mupar, asym_index, dealloc_hessian
      USE quantities, ONLY: fbdy, dealloc_quantities
      USE dumping_mod, ONLY: write_output
      USE diagnostics_mod, ONLY: write_profiles, dealloc_diagnostics
      USE perturbation, ONLY: lresistive, nsin, wout_file
      USE metrics, ONLY: cleanup_metric_elements, dealloc_full_lower_metrics
      USE blocktridiagonalsolver_s, ONLY: finalize
      USE prof_mod, ONLY: profstat
      USE shared_data, ONLY: lverbose

      IMPLICIT NONE

!  Declare Arguments
      TYPE (siesta_run_class), POINTER :: this
      LOGICAL, INTENT(in)              :: finalize_mpi
      LOGICAL, INTENT(in)              :: close_wout

!  local variables.
      REAL (rprec)                     :: ton
      REAL (rprec)                     :: toff
      REAL (rprec)                     :: skstoff
      INTEGER                          :: i

!  Start of executable code
      CALL second0(toff)
      skstoff = toff
      time_total = toff - this%skston
#if defined(SKS)
      total_time = skstoff - this%skston
#endif

#if defined(SKS)
      IF (OUTPUT_TIMINGS) THEN
         CALL GetTimes
      END IF
#endif

      IF (iam .EQ. 0) THEN
         DO i = 6, 33, 27
            IF (.not.lverbose .and. i .eq. 6) CYCLE
            WRITE (i, 1000) nprecon, levmarq_param, mupar, asym_index
#if defined(SKS)
            WRITE (i,'(a,i5)') ' Number processors: ', nprocs
#endif
            WRITE (i, 1001)
            WRITE (i, 1002) time_total,               fbdy(1),                 &
                            time_init,                fbdy(2),                 &
                            time_diag_prec,           fbdy(3),                 &
                            time_block_prec,          fbdy(4),                 &
                            time_factor_blocks,       fbdy(5),                 &
                            time_toijsp,              fbdy(6),                 &
                            time_tomnsp,                                       &
                            gmres_time,                                        &
                            conj_grad_time,                                    &
                            time_update_pres,         fbdy(7),                 &
                            time_update_bfield,       fbdy(8),                 &
                            time_current,             fbdy(9),                 &
                            get_force_harmonics_time, fbdy(10),                &
                            time_update_force,        fbdy(11),                &
                            time_update_upperv,       fbdy(12),                &
                            time_update_state,        fbdy(13),                &
                            time_funci,                                        &
                            time_apply_precon,                                 &
                            (diag_add_pert_time + block_add_pert_time),        &
                            time_init_state
            WRITE (i,*)
            IF (PARSOLVER) THEN
               WRITE (i,*)        'PARSOLVER      : NSCALED'
            ELSE
               WRITE (i,*)        'PARSOLVER      : ScaLAPACK'
            END IF
            WRITE (i,*)           'PARFUNCTISL    :', PARFUNCTISL
            WRITE (i,'(1x,a,L2)') 'PARGMRES       :', PARGMRES
#if defined(SKS)
            WRITE (i,*)           'OUTPUT_TIMINGS :', OUTPUT_TIMINGS
#endif
            WRITE (i,*)
            WRITE (i,'(a,3x,1p,e10.3)')' TIME DIVB  : ', time_divb
            WRITE (i,'(a,3x,1p,e10.3)')' TIME DIVJ  : ', time_divj
            WRITE (i,'(a,3x,1p,e10.3)')' TIME BGRADP: ', time_bgradp
            WRITE (i,'(a,3x,1p,e10.3)')' TIME BDOTJ : ', time_bdotj
            WRITE (i,*)
            WRITE (i,*) 'M (block size)         :', mblk_size
            WRITE (i,*) 'N (block rows)         :', nsin
            WRITE (i,*) 'P (processors)         :', nprocs
#if defined(SKS)
!            WRITE (i,'(1x,a,1p,2e12.4)')'Total time (max/min)   :',total_time_max, total_time_min
!            WRITE (i,*)
#endif
         END DO
      ENDIF

      CALL write_output(wout_file, niter, close_wout)
      CALL write_profiles(fsq_total1)  ! SPH: write jpmn, jbsupXmn, kbsubXmn, jvsupXmn profiles
      IF (iam .EQ. 0) THEN
         IF (lverbose) PRINT *,' Writing output to "siesta_profiles.txt" is finished!'
         CLOSE (UNIT=33)
      ENDIF

!CLEAN UP ARRAYS
      CALL cleanup_metric_elements

      CALL dealloc_quantities
      CALL dealloc_hessian
      CALL dealloc_diagnostics
!      IF (lresistive) THEN
         CALL dealloc_full_lower_metrics
!      END IF
      IF (ALLOCATED(xc)) THEN
         DEALLOCATE(xc, gc, col_scale)
      END IF

#if defined(MPI_OPT)
      IF (LSCALAPACK) THEN
         CALL blacs_barrier(icontxt, 'All')
         CALL blacs_gridexit(icontxt)

         CALL blacs_exit(0)
         IF ((myrow .eq. 0) .and. (mycol .eq. 0)) THEN
            CALL profstat()
         END IF

#if defined(SKS)
      ELSE
         IF (SKSDBG) THEN
            WRITE(TOFU,*) 'Called finalizeRemap and Finalize'
            FLUSH(TOFU)
         END IF
         CALL finalizeRemap
         CALL Finalize(finalize_mpi)
#endif
      END IF
#endif

      DEALLOCATE(this)

1000  FORMAT(/,' nprecon: ',i3,' LM parameter: ',1pe9.2,' mu||: ',1pe9.2,      &
             ' Symmetry Index: ',1pe9.2)
1001  FORMAT(/,'==============================',14x,'=======================', &
             /,                                                                &
             /,'  TIMING INFORMATION          ',14x,' RMS BOUNDARY FORCES',    &
             /,                                                                &
             /,'==============================',14x,'=======================')
1002  FORMAT(' Total runtime  : ',  f12.3,15x,'fs(1,m=1)  :',1pe10.2/,         &
             ' Initialization : ',0pf12.3,15x,'fs(2,m=1)  :',1pe10.2/,         &
             ' Diagonal prec  : ',0pf12.3,15x,'fs(2,m!=1) :',1pe10.2/,         &
             ' Compute blocks : ',0pf12.3,15x,'fu(1,m=1)  :',1pe10.2/,         &
             ' Factor blocks  : ',0pf12.3,15x,'fu(2,m=1)  :',1pe10.2/,         &
             ' Toijsp         : ',0pf12.3,15x,'fu(2,m!=1) :',1pe10.2/,         &
             ' Tomnsp         : ',0pf12.3,/,                                   &
             ' GMRES          : ',0pf12.3,/,                                   &
             ' Conj Gradient  : ',0pf12.3,//,                                  &
             ' SUBROUTINES     ',/,                                            &
             ' Update Pressure: ',0pf12.3,15x,'fv(1,m=0)  :',1pe10.2/,         &
             ' Update Bfield  : ',0pf12.3,15x,'fv(2,m=0)  :',1pe10.2/,         &
             ' CV Currents    : ',0pf12.3,15x,'fv(2,m!=0) :',1pe10.2/,         &
             ' Force Harmonics: ',0pf12.3,15x,'fu(ns)     :',1pe10.2/,         &
             ' Update Force   : ',0pf12.3,15x,'fu(ns-1)   :',1pe10.2/,         &
             ' Update UpperV  : ',0pf12.3,15x,'fv(ns)     :',1pe10.2/,         &
             ' Update State   : ',0pf12.3,15x,'fv(ns-1)   :',1pe10.2/,         &
             ' Funct Island   : ',0pf12.3/,                                    &
             ' Apply Precon   : ',0pf12.3/,                                    &
             ' Add Perturb    : ',0pf12.3/,                                    &
             ' Init State     : ',0pf12.3,                                     &
             /,'==============================',14x,'=======================')
      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets quanties based on VMEC values.
!>
!>  This method loads and sets variables based on the VMEC equilibrium. VMEC 
!>  controls the metric elements and coordinate system jacobian.
!>
!>  @param[inout] this A @ref siesta_run_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_set_vmec(this)
      USE perturbation, ONLY: nsin, mpolin, ntorin, wout_file
      USE metrics, ONLY: init_metric_elements
      USE timer_mod, ONLY: init_metric_elements_time, init_quantities_time
      USE quantities, ONLY: init_quantities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (siesta_run_class), INTENT(inout) :: this

!  local variables.
      REAL (rprec)                           :: ton
      REAL (rprec)                           :: toff

!  Start of executable code
#if defined(SKS)
      CALL second0(ton)
#endif
      CALL init_metric_elements(nsin, mpolin, ntorin, wout_file)
#if defined(SKS)
      CALL second0(toff)
      init_metric_elements_time = init_metric_elements_time + (toff - ton)
#endif

#if defined(SKS)
      CALL second0(ton)
#endif
      CALL init_quantities
#if defined(SKS)
      CALL second0(toff)
      init_quantities_time = init_quantities_time + (toff - ton)
#endif

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Solves the SIESTA equilibrium.
!>
!>  This method converges the siesta_run_class equilibrium. Initializes the 
!>  evolution module variables and converges the force residual using a diagonal 
!>  preconditioner. Applies the external island perturbation if possible. Loads 
!>  restart file before convergence.
!>
!>  @param[inout] this A @ref siesta_run_class instance.
!>  @returns True if the convergece was sucessful and false otherwise.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_converge(this)
      USE perturbation, ONLY: read_restart_file
      USE shared_data, ONLY: nprecon
      USE timer_mod, ONLY: time_init, init_evolution_time,                     &
                           init_quantities_time, converge_diagonal_time,       &
                           converge_blocks_time
      USE evolution, ONLY: init_evolution, converge_diagonal, converge_blocks
      USE descriptor_mod, ONLY: iam, nprocs, DIAGONALDONE
      USE shared_data, ONLY: ndims, mblk_size
      USE island_params, mpol=>mpol_i, ntor=>ntor_i
      USE perturbation, ONLY: nsin, mpolin, ntorin, ftol, lrestart, wout_file
      USE nscalingtools, ONLY: initRemap
      USE hessian, ONLY: InitHess, levmarq_param, levmarq_param0
      USE quantities, ONLY: init_fields
      USE metrics, ONLY: sqrtg

      IMPLICIT NONE

!  Declare Arguments
      TYPE (siesta_run_class), INTENT(inout) :: this

!  local variables.
      REAL (rprec)                           :: ton
      REAL (rprec)                           :: toff

!  Start of executable code
      IF (lrestart) CALL read_restart_file(nprecon)
      CALL second0(toff)
      time_init = toff - this%skston

      CALL second0(ton)
      CALL init_evolution
      mblk_size = ndims*(mpol + 1)*(2*ntor + 1)
#if defined(SKS)
      CALL initRemap(mpolin, ntorin, nsin, nprocs, iam)
#endif
      CALL InitHess
      CALL second0(toff)
      init_evolution_time = init_evolution_time + (toff - ton)

#if defined(SKS)
      CALL second0(ton)
#endif
      IF (.not.lrestart) THEN
         CALL Init_Fields
      END IF
      DEALLOCATE (sqrtg)
#if defined(SKS)
      CALL second0(toff)
      init_quantities_time = init_quantities_time + (toff - ton)
#endif

!  Converge initial residues with diagonal preconditioner
      DIAGONALDONE = .false.
#if defined(SKS)
      CALL second0(ton)
#endif
      CALL converge_diagonal(wout_file, ftol)
#if defined(SKS)
      CALL second0(toff)
      converge_diagonal_time = converge_diagonal_time + (toff - ton)
#endif
      DIAGONALDONE=.TRUE.

!  Converge using block preconditioner
#if defined(SKS)
      CALL second0(ton)
#endif
      CALL converge_blocks(wout_file, ftol)
#if defined(SKS)
      CALL second0(toff)
      converge_blocks_time = converge_blocks_time + (toff - ton)
#endif

      END SUBROUTINE

      END MODULE
