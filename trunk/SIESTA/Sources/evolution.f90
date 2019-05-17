!>  \brief Module for controlling iteration ("time") evolution of MHD convergence sequence.
!!  Updates the CONTRAVARIANT Fourier components of the "velocity" field (displacement X delt)
      MODULE evolution
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE shared_data
      USE shared_functions
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,          &
          nfp=>nfp_i, mnmax=>mnmax_i, hs=>hs_i
      USE timer_mod
      USE siesta_namelist, ONLY: eta_factor
      USE descriptor_mod, ONLY: iam, nprocs
      USE nscalingtools, ONLY: SKSDBG, PARSOLVER, PARFUNCTISL, MPI_ERR, &
                               startglobrow, endglobrow, rcounts, disp
      USE mpi_inc
      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: n_bsupu_full=1
      INTEGER, PARAMETER :: ndamp = 20              !<number smoothing steps for conj grad
      INTEGER, PARAMETER :: ncongrad_steps = 100    !<max number of conj grad steps
      LOGICAL            :: l_EPAR                  !<=TRUE when applying parallel tearing perturbations 
      LOGICAL            :: l_Gmres                 !<=TRUE performs GMRES iterations 
      LOGICAL            :: l_conjmin               !<=TRUE if conjugate gradient lowers fsq_min

      INTEGER            :: nprint_step=1
      INTEGER            :: nfun_calls              !<average number of funct_island calls per processor
      INTEGER            :: js, moff, noff

      REAL(dp)           :: tau_damp, delt_cg=1
      REAL(dp)           :: fsqprev, fsqprev1, levm_loc
      REAL(dp)           :: fsq_min, fsq_ratio, fsq_ratio1 
      REAL(dp)           :: otau(ndamp)
      REAL(dp), DIMENSION(:), ALLOCATABLE :: xcdot
            
#if defined(PETSC_OPT)
      INTEGER(KIND=SELECTED_INT_KIND(10)) Apc
      INTEGER(KIND=SELECTED_INT_KIND(10)) subApc(1048)
      INTEGER(KIND=SELECTED_INT_KIND(10)) ksp
      INTEGER(KIND=SELECTED_INT_KIND(10)) subksp(1048)
#endif
!-----------------------------------------------
      CONTAINS
      

!>  \brief Performs initial convergence of force residuals using diagonal preconditioner only.
!!  Terminates when the force residual stops changing by more than a factor of 2
      SUBROUTINE converge_diagonal (wout_file, ftol)
      USE siesta_namelist, ONLY: ladd_pert, l_output_alliter   !R. Sanchez, added Jul '10: L_OUTPUT_ALLITER Dump info on file at each timestep
      USE perturbation, ONLY: add_perturb
      USE dumping_mod, ONLY: write_output
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(IN) :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(dp), INTENT(IN)      :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL            :: l_iterate            !<call evolve until fsq_ratio stops decreasing
      REAL(dp)           :: fsq_block, t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!

!Test tri-diag structure
!     CALL TriDiag_Test (xc, gc, 1.E-8_dp)

      fsq_block = MAX(ftol, fsq_prec)
      l_iterate = .TRUE.
      l_Gmres   = .TRUE.
      nprecon_type = PREDIAG
      IF (nprecon .NE. 0) nprecon_type = PREBLOCK
      niter = 1
      levm_scale = 1    !MRC
      in_hess_nfunct=0; out_hess_nfunct=0
      
!Use column-scaled fsq_total here      

      DO WHILE (l_iterate)
         CALL second0(t1)
         CALL evolve
         CALL second0(t2)
#if defined(SKS)         
         diag_evolve_time=diag_evolve_time+(t2-t1)
#endif         
         IF (l_output_alliter) CALL write_output (wout_file, niter, .FALSE.)   !MRC Do not close wout_file

         l_iterate = (fsq_ratio1.LE.0.5_dp .AND. fsq_total1.GT.fsq_block)

         IF (ladd_pert .AND. fsq_total1.lt.100*fsq_block) THEN
            l_init_state=.TRUE.
            CALL second0(t1)
            CALL add_perturb(xc, getwmhd)
            CALL second0(t2)
            diag_add_pert_time=diag_add_pert_time+(t2-t1)
            ladd_pert = .FALSE.
         END IF
         niter = niter+1
      END DO

      END SUBROUTINE converge_diagonal


!>  \brief Performs convergence of force residuals using block preconditioner.
!!  Terminates when the force residual drops below ftol or the number of
!!  iterations exceeds a specified value (niter_max)
!!  Applies an external island perturbation if not previously done.
      SUBROUTINE converge_blocks (wout_file, ftol)
      USE siesta_namelist, ONLY: l_output_alliter, ladd_pert
      USE perturbation, ONLY: add_perturb, write_restart_file, niter_max   ! R. Sanchez, added Jul '10: Dump info on file at each timestep
      USE dumping_mod, ONLY: write_output
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(IN)  :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(dp), INTENT(IN)       :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nrow
      REAL(dp),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL            :: l_iterate            !<controls iteration sequence
      REAL(dp)           :: t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!
      l_Gmres   = .FALSE.
      nprecon_type = PREBLOCK
      nprecon = 1
      niter_max = niter+niter_max
      l_iterate = (niter.LT.niter_max) .AND. (fsq_total1.GT.ftol) 
      nrow = 0

      DO WHILE (l_iterate)
         CALL second0(t1)
         CALL evolve
         CALL second0(t2)
#if defined(SKS)
         block_evolve_time=block_evolve_time+(t2-t1)
#endif
         IF (l_output_alliter) CALL write_output (wout_file, niter, .FALSE.)   !MRC Do not close wout_file

         l_iterate = (niter.LT.niter_max) .AND. (fsq_total1.GT.ftol) 

         IF (.NOT.l_conjmin .AND. .NOT.l_gmres) THEN
            IF (IAM .EQ. 0) PRINT *,'CONJ GRADIENT UNABLE TO MAKE PROGRESS'
            l_iterate=.FALSE.
         ELSE IF (nprecon.GT.3 .AND. fsq_ratio.GT.1.E5_dp) THEN
            IF (IAM.EQ.0 .AND. lverbose) PRINT 100, fsq_ratio1
            l_iterate=.FALSE.
         END IF
!STOP IF LACK OF PROGRESS
         IF (nprecon.GT.2 .AND. ABS(1-fsq_ratio1).LT.1.E-2_dp) THEN
            levm_scale = levm_scale/3
            nrow = nrow+1
            IF (nrow.EQ.3 .AND. l_iterate) THEN
               l_iterate=.FALSE.
               IF (IAM.EQ.0 .AND. lverbose) PRINT 100, fsq_ratio1
            END IF
         ELSE 
            nrow = 0
         END IF
!In case we didn't add it in diag loop
         IF (ladd_pert .AND. fsq_total1.LT.100*fsq_prec) THEN
           l_init_state=.TRUE.
           CALL second0(t1)
           CALL add_perturb(xc, getwmhd)
           CALL second0(t2)
           block_add_pert_time=block_add_pert_time+(t2-t1)
           ladd_pert = .FALSE.
!TO OUTPUT RIGHT AFTER PERT IS APPLIED, SET l_iterate=.FALSE. HERE
         END IF
         nprecon = nprecon+1
         niter = niter+1
      END DO

100   FORMAT (' FSQ RATIO: ',1p,e10.2, ' SIESTA CAN NOT CONVERGE FURTHER!')

      END SUBROUTINE converge_blocks

!>   \brief  Initializes variables and pointers prior to calling evolve
      SUBROUTINE init_evolution
      USE hessian, ONLY: levmarq_param0, levmarq_param, l_Compute_Hessian
      USE quantities, ONLY: jvsupsmncf,jvsupumnsf,jvsupvmnsf,           &
                            jvsupsmnsf,jvsupumncf,jvsupvmncf,           &
                            fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            fsubsmnsf, fsubumncf, fsubvmncf,            &
                            fsupsmncf, fsupumnsf, fsupvmnsf,            &
                            fsupsmnsf, fsupumncf, fsupvmncf
      USE shared_functions, ONLY: init_ptrs
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat, n1
!-----------------------------------------------

!     "natural" boundary condition at origin
!     l_push_s/l_push_u can't both be .TRUE. (linearly dependent => 0 eigenvalue in Hessian)

      l_natural = .TRUE.  !.FALSE.

!     control evolve displacements (m=1 for s,u; m=0 for v) at s=0
      l_push_s = .FALSE.
      l_push_u = .FALSE.
      l_push_v = .FALSE.

      IF (l_push_u) l_push_s = .FALSE.                                    !avoid linear dependency
!      l_push_edge=.TRUE.
      l_push_edge=.FALSE.

      scale_s = hs
      scale_u = SQRT(scale_s)

      ste = 0
      l_EPAR = .FALSE.

      nfun_calls = 0
      fsqprev = -1;  fsqprev1 = -1
      fsq_min = 1.E20_dp
      niter    = 0
      wtotal0  = -1
      etak_tol = 1.E-01_dp                                                !initial Newton tolerance parameter
      delta_t = 1
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      l_getwmhd= .FALSE.                                                  !set in getwmhd function
      levmarq_param = levmarq_param0
      fsq_total = 1
      fsq_total1= 1

      CALL ASSERT(ndims.EQ.3 .OR. ndims.EQ.6,'WRONG ndims!')
      CALL ASSERT(mnmax.EQ.(1+mpol)*(2*ntor+1),'WRONG mnmax!')

      n1 = ns*mnmax
      neqs = ndims*n1
!      PRINT *,"neqs: ", neqs, __LINE__, __FILE__, ndims

      ALLOCATE(xc(neqs), col_scale(0:mpol,-ntor:ntor,ndims,ns), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocate xc failed!')
      xc = 0
      col_scale = 1
      CALL init_ptrs (xc, jvsupsmncf, jvsupumnsf, jvsupvmnsf,           &
                          jvsupsmnsf, jvsupumncf, jvsupvmncf)

      ALLOCATE(gc(neqs), gc_sup(neqs), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocate gc failed!')
      gc = 0
      CALL init_ptrs (gc_sup, fsupsmncf, fsupumnsf, fsupvmnsf,          &
                              fsupsmnsf, fsupumncf, fsupvmncf)
      CALL init_ptrs (gc,     fsubsmncf, fsubumnsf, fsubvmnsf,          &
                              fsubsmnsf, fsubumncf, fsubvmncf)
                              
      l_Compute_Hessian = .FALSE.
     
      col_scale(:,:,1,:) = scale_s
      col_scale(:,:,2,:) = scale_u
      IF (lasym) THEN
      col_scale(:,:,4,:) = scale_s
      col_scale(:,:,5,:) = scale_u
      END IF

      END SUBROUTINE init_evolution

!>   \brief Evolves one step toward MHD equilibrium
      SUBROUTINE evolve
      USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,         &
                            jvsupsmnsf, jvsupumncf, jvsupvmncf,         &
                            fsubsmncf,  fsubumnsf,  fsubvmnsf,          &
                            fsubsmnsf,  fsubumncf,  fsubvmncf
      USE gmres, ONLY: gmres_fun
      USE hessian
      USE siesta_namelist, ONLY: eta_factor, lresistive, ftol
      USE perturbation, ONLY: write_restart_file, niter_max
      USE diagnostics_mod, ONLY: bgradp_rms, toroidal_flux0, bgradp,    &
                                 tflux
      USE siesta_state, ONLY: update_state
      USE siesta_init, ONLY: init_state
#if defined(SKS)
      USE blocktridiagonalsolver_s, ONLY: RefactorHessian      
      USE timer_mod, ONLY: gmres_time
#endif
!
!     Evolves the equilibrium equations one time step and updates delta_t
!     
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: levscan(1) = (/ 0.25_dp /) 
      REAL(dp)    :: skston, skstoff
      INTEGER     :: iprint, m1, i, n1, n2
      REAL(dp)    :: v2tot, f1, f2
      REAL(dp)    :: t1, t2, lm0
      LOGICAL     :: lprint, ldiagonal, lblock
      INTEGER, DIMENSION (2) :: snum_funct, rnum_funct 
      CHARACTER(LEN=20)    :: str_resistive
!-----------------------------------------------
!
!     Update COVARIANT Fourier components of velocity field
!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD BY P. GARABEDIAN
!
      l_init_state = .TRUE.
      IF (niter .LE. 1) levm_loc = 1 ! MRC needed to reset levm_loc
!     l_Gmres = (nprecon_type .EQ. PREDIAG)  ! Turn on CG after initial diagonal preconditioner
!      l_Gmres = l_Gmres .OR. fsq_total1.LT.1.E-9_dp   ! Turn on Gmres after CG reduces residue
      l_Gmres = .TRUE.    !TURNS OFF CONJ-GRAD ALWAYS
!
!     Determine levenberg and mu|| parameter guesses based on ||F||
!
!       IF (nprecon.GE.1 .AND. (fsq_total1.LT.1.E-8_dp .OR. .NOT.l_Gmres)) THEN
      IF (fsq_total1.LT.1.E-6_dp) THEN
!      IF (fsq_total1.LT.1.E-8_dp) THEN
!        PSEUDO-TRANSIENT CONTINUATION: TURN OFF FOR CONSTANT levmarq_param, muPar option
         f2 = fsq_total1*1.E12_dp
!THIS WORKS FOR BOTH LASYM = T/F         
!         f1 = MIN(one, f2)         

!ASYM WORKS BETTER LIKE THIS - DON'T KNOW WHY
         f2 = one/(one+f2)
         IF (lasym) f2 = f2/4
         f1 = fsq_total1**f2

         f2 = 1
         IF (fsq_ratio1.GT.0.25_dp .AND. fsq_ratio1.LT.1.5_dp) f2 = 0.5_dp
         IF (fsq_ratio1.GT.1.5_dp) f2 = 1.5_dp
         IF (fsq_ratio1.GT.10._dp) f2 = 3
         levm_loc = MAX(0.1_dp,MIN(10._dp,f2*levm_loc))

!         f1 = MIN(one, f1*levm_scale)
         f1 = MIN(one, f1*levm_scale*levm_loc)

         levmarq_param = (levmarq_param0 - levm_ped)*f1 + levm_ped
         IF (levmarq_param0 .EQ. zero) levmarq_param = 0
         muPar = (muPar0-mu_ped)*f1 + mu_ped
         IF (muPar0 .EQ. zero) muPar=0

      ELSE
         levmarq_param = MAX(0.3_dp, levmarq_param0)
         muPar = muPar0
      END IF

!
!     EXCEPT FOR THE INITIAL STEP (WHERE ONLY INITIAL UNPRECONDITIONED FORCES
!     ARE NEEDED) COMPUTE FULL (OR DIAGONAL APPROXIMATION) HESSIAN
!
      IF (niter .GT. 1) THEN
         l_linearize = .TRUE.
         l_getfsq = .FALSE.
         l_ApplyPrecon = .FALSE.
         l_getwmhd = .FALSE.
         ldiagonal = (nprecon_type.EQ.PREDIAG .AND. niter.EQ.2)
         lblock    = (nprecon_type .EQ. PREBLOCK)

         IF (ldiagonal .OR. lblock) THEN

            CALL second0(skston)
            CALL Compute_Hessian_Blocks (funct_island, ldiagonal)
            CALL second0(skstoff)
            IF (ldiagonal) comp_diag_elements_time=comp_diag_elements_time+(skstoff-skston)
            IF (lblock)    compute_hessian_time=compute_hessian_time+(skstoff-skston)
            
         END IF

      END IF

!
!     Reset run-control logicals
!     (may have been reset to false in Compute_Hessian=>funct_island call)

      l_init_state = .TRUE.   
      l_ApplyPrecon = .FALSE.
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      fsq_lin = -1

!
!     Choose time-stepping algorithm
!     Update COVARIANT Fourier components the forces to evolve the contravariant components of the displacement
!
      IF (niter .EQ. 1) THEN
         l_init_state=.TRUE.
         l_getwmhd=.TRUE.
         CALL second0(skston)
         CALL funct_island
         CALL second0(skstoff)
         evolve_funct_island_time=evolve_funct_island_time+(skstoff-skston)
         l_getwmhd=.FALSE.

      ELSE IF (l_Gmres) THEN
        CALL second0(skston)
        IF (nprecon .NE. 1) THEN
           etak_tol = MIN(0.1_dp, 1.E8_dp*fsq_total)
           IF (fsq_total <= 0) etak_tol = 0.1_dp
        END IF
        CALL gmres_fun
        CALL second0(skstoff)
        gmres_time=gmres_time+(skstoff-skston)
      END IF

      IF (.NOT.l_Gmres .AND. niter.GT.1) THEN
         CALL second0(skston)
         lm0 = levmarq_param
         DO i = 1, 1
!TURN OFF UNTIL RefactorHessian PROVIDED FOR LSCALAPACK=T
!         DO i = 1, 1+SIZE(levscan)
#if defined(SKS)
            IF (i .GT. 1) THEN 
               levmarq_param = levscan(i-1)*lm0
               IF (levmarq_param .GT. levmarq_param0) levmarq_param = levmarq_param0*1.E-2_dp
               IF (levmarq_param .LE. levm_ped) levmarq_param = 100*levm_ped
               IF (iam .EQ. 0 .and. lverbose) PRINT 930,' Refactoring Hessian: LM = ',levmarq_param
               CALL RefactorHessian(levmarq_param)
            END IF
#endif
            CALL Conjugate_Grad(ftol)
            IF (l_conjmin) THEN
               IF (i .GT. 1) levm_scale = levm_scale*levmarq_param/lm0
               EXIT
            ENDIF 
         END DO
         CALL second0(skstoff)
         conj_grad_time=conj_grad_time+(skstoff-skston)
      END IF
 930  FORMAT (/,a,1pe12.3)
      v2tot = SQRT(hs*SUM(xc*xc))
!      IF (v2tot.lt.1.E-5_dp .and. .not.l_Gmres .and. niter.gt.1) THEN
!      END IF
        
      lprint = MOD(niter, nprint_step).EQ.0 .OR. niter.EQ.1
      l_update_state = .TRUE.
      CALL update_state(lprint, fsq_total1, zero)
      l_update_state=.FALSE.
        
!
!     Compute force component residues  |Fsubi|**2
!
      IF (fsq_total.LT.zero .AND. iam.EQ.0)                             &
         PRINT *, 'Fsq_Total = ', fsq_total,' < 0!'
      IF (fsqprev .GT. zero) THEN
         fsq_ratio=fsq_total/fsqprev;  fsq_ratio1=fsq_total1/fsqprev1
      ELSE
         fsq_ratio=0;  fsq_ratio1=0;
      END IF
      fsqprev=fsq_total; fsqprev1=fsq_total1

!     IF v-step is too large, recompute Hessian
!      IF (v2tot*delta_t .ge. one) fsq_ratio = 1

!
!     Convert output to REAL (VMEC) units, divide B-fields by b_factor
!     p by p_factor, WMHD by p_factor

#if defined(MPI_OPT)   
      snum_funct(1)=in_hess_nfunct; snum_funct(2)=out_hess_nfunct
      CALL MPI_ALLREDUCE(snum_funct, rnum_funct, 2, MPI_INTEGER,        &
                         MPI_MAX, SIESTA_COMM, MPI_ERR)
      nfun_calls=(rnum_funct(1)+rnum_funct(2))
#else 
      nfun_calls=(in_hess_nfunct+out_hess_nfunct)
#endif

     IF (niter .EQ. 1) THEN
        CALL bgradp(startglobrow, endglobrow)
        CALL tflux
     END IF


     IF (lprint) THEN
     
        CALL gather_array(gc)                                            !Needed to compute maximum forces for printout

        IF (niter.GT.1 .AND. iam.EQ.0) THEN
          DO iprint = 6, unit_out, unit_out-6
             IF (.NOT.lverbose .AND. iprint.EQ.6) CYCLE
             WRITE(iprint,*)'SYMMETRIC DISPLACEMENTS AND FORCES'
          END DO
          CALL WRITE_BOUNDS(jvsupsmncf, fsubsmncf, 'X^S-max:', 'F_S-max:')
          CALL WRITE_BOUNDS(jvsupumnsf, fsubumnsf, 'X^U-max:', 'F_U-max:')
          CALL WRITE_BOUNDS(jvsupvmnsf, fsubvmnsf, 'X^V-max:', 'F_V-max:')
          IF (lasym) THEN
             DO iprint = 6, unit_out, unit_out-6
                IF (.NOT.lverbose .AND. iprint.EQ.6) CYCLE
                WRITE(iprint,*)'ASYMMETRIC DISPLACEMENTS AND FORCES'
             END DO
             CALL WRITE_BOUNDS(jvsupsmnsf, fsubsmnsf, 'X^S-max:', 'F_S-max:')
             CALL WRITE_BOUNDS(jvsupumncf, fsubumncf, 'X^U-max:', 'F_U-max:')
             CALL WRITE_BOUNDS(jvsupvmncf, fsubvmncf, 'X^V-max:', 'F_V-max:')
          END IF

          IF (lverbose) THEN
             PRINT 110, bsbu_ratio_s, jsju_ratio_s,                     &
                       (bs0(i), i=2,6), (bu0(i), i=2,6) 
             IF (lasym)                                                 &
             PRINT 112, bsbu_ratio_a, jsju_ratio_a,                     &
                       (bs0(i), i=8,12), (bu0(i), i=8,12) 
             PRINT *," "
          END IF
        END IF

        DO iprint = 6, unit_out, unit_out-6
           IF ((.NOT.lverbose .AND. iprint.EQ.6) .OR. iam.NE.0) CYCLE
           IF (niter .EQ. 1) THEN
              IF (lresistive) THEN
                 str_resistive = "RESISTIVE RUN "
              ELSE
                 str_resistive = "NON-RESISTIVE RUN "
              END IF
              WRITE (iprint, 55) str_resistive, eta_factor, lasym,      &
                                 l_push_s, l_push_u, l_push_v, l_push_edge
              WRITE (iprint, 60) delta_t, etak_tol, l_Hess_sym
              WRITE (iprint, 65) levmarq_param0, mupar0
              WRITE (iprint, 70) neqs, ns, mpol, ntor, nu_i, nv_i,      &
                                 ngmres_steps
           END IF
           
           IF (fsq_lin .EQ. -1) THEN
              WRITE (iprint, 100)                                       &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,  &
                              fsq_total1, fsqvs, fsqvu,                 &
                              fsqvv, v2tot*delta_t, nfun_calls
           ELSE
              WRITE (iprint, 102)                                       &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,  &
                              fsq_total1, fsq_lin, fsqvs, fsqvu,        &
                              fsqvv, v2tot*delta_t, nfun_calls
           END IF
        END DO
     END IF

 50  FORMAT(' WMHD: ', 1pe13.6,' <BETA>: ',1pe11.4,                     &
            ' TFLUX: ',1pe11.4,' B.GRAD-P (rms): ', 1pe11.4,/,21('-'))
 55  FORMAT(1x, a,' ETA_FACTOR: ', 1pe10.2,' LASYM: ',l1,               &
            ' L_PUSH_S: ',l1, ' L_PUSH_U: ',l1,' L_PUSH_V: ',l1,        &
			' L_PUSH_EDGE: ',l1)
 60  FORMAT(' DELTA_T: ',1pe10.2,' ETA_K: ',                            &
             1pe10.2,' HESSIAN SYM: ', l1)
 65  FORMAT(' LEVMARQ_PARAM: ',1pe10.2,' MU_PAR: ',1pe10.2)
 70  FORMAT(' NEQS: ', i6,' NS: ',i4,' MPOL: ',i4,' NTOR: ',i4,         &
             ' NTHETA: ', i4,' NZETA: ',i4, ' NGMRES-STEPS: ', i4,//,   &
             ' NITER (W-W0)/W0*1E6    F2(MHD)    F2(LIN)',              &
             '    F2SUBS     F2SUBU     F2SUBV     |V|rms     NCALLS')
 100 FORMAT(1x,i5,1x, 1pe12.4, 2x, 1x,1pe10.3, 2x,9 ('-'),              &
           4(1x,1pe10.3),i9)
 102 FORMAT(1x,i5,1x, 1pe12.4, 2x, 6(1x,1pe10.3),i9)
 110 FORMAT(' |B^s-r12*B^u|/|B^s+r12*B^u| (m=1,r->0,  sym) : ',1pe10.3,/, &
            ' |J^s-r12*J^u|/|J^s+r12*J^u| (m=1,r->0,  sym) : ',1pe10.3,/, &
            ' JBSUPSH(JS=2-6,M=1)/r12  (sym): ',1p5e10.3,/,             &
            ' JBSUPUH(JS=2-6,M=1)           : ',1p5e10.3)
 112 FORMAT(' |B^s+r12*B^u|/|B^s-r12*B^u| (m=1,r->0, asym) : ',1pe10.3,/, &
            ' |J^s+r12*J^u|/|J^s-r12*J^u| (m=1,r->0, asym) : ',1pe10.3,/, &
            '-JBSUPSH(JS=2-6,M=1)/r12 (asym): ',1p5e10.3,/,             &
            ' JBSUPUH(JS=2-6,M=1)           : ',1p5e10.3)

     CALL second0(skston)
     IF (fsq_total1 .LT. fsq_min) THEN
        fsq_min = fsq_total1
        CALL write_restart_file
     END IF
     CALL second0(skstoff)
     evolve_restart_file_time=evolve_restart_file_time+(skstoff-skston)

!    Add resistive part of B perturbation
     CALL second0(skston)
     IF (lresistive .AND. nprecon.GE.1 .AND. fsq_total1.GT.fsq_res) THEN
        l_getfsq = .FALSE.
        CALL add_resistive_E
     END IF
     CALL second0(skstoff)
     evolve_add_resistive_E_time=evolve_add_resistive_E_time+(skstoff-skston)

     END SUBROUTINE evolve


!>  \brief Writes maximum values for velocity/forces and radial,m,n, locations
      SUBROUTINE Write_Bounds(vsupXmn, fsubXmn, vlabel, flabel)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)             &
                                    :: vsupXmn, fsubXmn
      CHARACTER(LEN=*), INTENT(IN)  :: vlabel, flabel
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: nm0=1, nn0=2, ns0=3
      INTEGER     :: imaxlocv(3), ilbound(3), imaxlocf(3), iprint
      REAL(dp)    :: vmax, fmax
!-----------------------------------------------

      ilbound = LBOUND(vsupXmn)-1
      imaxlocv = MAXLOC(ABS(vsupXmn)) + ilbound
      vmax = vsupXmn(imaxlocv(1),imaxlocv(2),imaxlocv(3))
      CALL ASSERT(ABS(vmax).EQ.MAXVAL(ABS(vsupXmn)),'vmax WRONG')
      imaxlocv = imaxlocv - ilbound
      imaxlocf = MAXLOC(ABS(fsubXmn)) + ilbound
      fmax = fsubXmn(imaxlocf(1),imaxlocf(2),imaxlocf(3))
      CALL ASSERT(ABS(fmax).EQ.MAXVAL(ABS(fsubXmn)),'fmax WRONG')
      imaxlocf = imaxlocf - ilbound

      DO iprint = 6, unit_out, unit_out-6
         IF ((.NOT.lverbose .AND. iprint.EQ.6) .OR. iam.NE.0) CYCLE
         WRITE(iprint, 45)                                              &
            TRIM(vlabel),vmax,imaxlocv(ns0),imaxlocv(nm0)-1,            &
                         imaxlocv(nn0)-ntor-1,                          &
            TRIM(flabel),fmax,imaxlocf(ns0),imaxlocf(nm0)-1,            &
                         imaxlocf(nn0)-ntor-1
      END DO
     
 45   FORMAT(2(1x,a,1x,1pe10.2,' AT JS: ',i4,' M: ',i4,' N: ',i4,2x))
     
      END SUBROUTINE Write_Bounds


!>  \brief Parallel routine to evolve state toward equilibrium, using Conjugate Gradients
      SUBROUTINE Conjugate_Grad (ftol)
      USE stel_kinds
      USE hessian, ONLY: mupar
      USE siesta_state, ONLY: update_state, update_state
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp) :: ftol
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.50_dp
      INTEGER, PARAMETER :: ndamp = 10
      INTEGER     :: icount, iter
      REAL(dp)    :: b1, fac, fsq_start, fsq_min0, fsq_save, fsq_conv,   &
                     wmhd_min, wmhd_start, muParS, dtau, delt_cg, otav
      REAL(dp), ALLOCATABLE  :: xcstart(:)
!-----------------------------------------------
      IF (fsq_total1.LE.ftol .OR. nprecon.EQ.0) THEN
         IF (nprecon .EQ. 0) nprecon = 1
         RETURN
      END IF

      l_getfsq = .TRUE.
      l_linearize = .FALSE.
      l_ApplyPrecon = .TRUE.
      l_init_state = .TRUE.
      l_PrintOriginForces = .FALSE.
      l_getwmhd=.TRUE.
 
! 
!     Linearize perturbations around B0 and p0
!     dB = curl(xc X B0)
!     dp = dp(p0,xc)
!
!     NON-LINEAR forces 
!     J = J0+dJ,  B = B0+dB, p = p0+dp

      ALLOCATE(xcdot(neqs), xc0(neqs), stat=icount) 
      CALL ASSERT(icount.eq.0,'ALLOCATION FAILED IN CONJ_GRAD')

      xcdot = 0
      xc = 0
      xc0 = 0

      delt_cg = 1
      dtau = 0 
      fsq_save = 0

      IF (iam .EQ. 0 .and. lverbose) PRINT 90
90    FORMAT(1x,'CONJUGATE GRADIENT CONVERGENCE SUMMARY',               &
             /,' -------------',/,1x,'ITER',7x,'FSQ_NL',10x,'||X||')

!      muParS = muPar; ! muPar=0

      DO icount = 1, 10*ncongrad_steps
!
!     LOOP OVER delt_cg TO FIND MINIMUM F(xc)
!
         CALL funct_island

         CALL update_taudamp (fsq_save, delt_cg)

         IF (icount .EQ. 1) THEN
            fsq_start  = fsq_total1
            fsq_min0   = fsq_total1
            fsq_conv  = fsq_min0
            fsq_save   = fsq_total
            wmhd_start = wtotal
            wmhd_min   = wtotal
            ALLOCATE (xcstart(SIZE(gc)))
            xcstart    = -gc
         END IF

         IF (fsq_total1 .LT. fsq_min0) THEN
            fsq_min0 = fsq_total1
            xc0 = xc
            IF (fsq_total1 .LT. ftol) EXIT
         END IF
!         ELSE IF (wtotal .lt. wmhd_min) THEN
!            wmhd_min = wtotal
!            xc0 = xc
!         END IF
         IF (iam.EQ.0 .AND. (icount.EQ.1 .OR. MOD(icount,10).EQ.0) .and. lverbose)     &
            PRINT 110, icount, fsq_total1, SQRT(SUM(xc*xc))

         IF (fsq_total1 .GT. 1.5_dp*fsq_min0) THEN
            delt_cg = 0.95_dp*delt_cg
            IF (delt_cg .LT. 0.001_dp) EXIT
            xc = xc0
            xcdot = 0
            IF (fsq_total1 .GT. 10*fsq_min0) THEN
               delt_cg = delt_cg*p5
               fsq_save = 0
               CYCLE
            END IF
         END IF

         IF (MOD(icount, 50) .EQ. 0) THEN
            IF (fsq_min0 .GT. fsq_conv/2) EXIT
            fsq_conv = fsq_min0
         END IF       

         CALL ASSERT(.NOT.l_linearize,'LINEARIZATION TURNED ON!')
         otav = SUM(otau)/ndamp
         dtau = delt_cg*otav/2

         b1 = one-dtau
         fac = one/(one+dtau)
!
!     THIS IS THE TIME-STEPPING ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD BY P. GARABEDIAN
!
!        CONJ GRAD 
         xcdot = fac*(b1*xcdot + delt_cg*gc)
         xc    = xc + delt_cg*xcdot
!
!        STEEPEST DESCENT (b1->0)
!         xc = b1*xc - delt_cg*gc

      END DO

!
!     UPDATE FIELD, PRESSURE PERTS AND UPDATE NEW NON-LINEAR STATE
!
      muParS = muPar; muPar = 0
      l_PrintOriginForces = .FALSE.

!      IF (fsq_min0 .EQ. fsq_start) xc0 = xcstart
!      CALL LineSearch(xc0, fsq_min0)

      l_conjmin = (fsq_min0 .LT. fsq_start)
      IF (l_conjmin) THEN
          xc = xc0
      ELSE
          xc = 0
      END IF

      IF (fsq_min0 .GT. 0.95_dp*fsq_start) l_conjmin=.FALSE.
      IF (.NOT.l_Gmres .AND. l_conjmin) l_PrintOriginForces = .TRUE.
	  l_ApplyPrecon = .FALSE.
      CALL funct_island

      l_PrintOriginForces = .FALSE.

      IF (muParS .NE. zero) fsq_min0 = fsq_total1

      IF (iam.EQ.0 .AND. l_conjmin) THEN
         IF (lverbose) WRITE (6, 100) icount, fsq_start, fsq_min0
         WRITE (unit_out,100) icount, fsq_start, fsq_min0
      END IF

      muPar = muParS

      DEALLOCATE (xcdot, xc0, xcstart)
 100  FORMAT(' ITER_CG: ', i3,' FSQ (START CG): ',1p,e12.3,' FSQ (END CG): ', 1pe12.3)
 110  FORMAT(i5, 2(3x,1pe12.3))

      END SUBROUTINE Conjugate_Grad

!>  \brief Updates damping parameter (tau) used by conjugate gradient routine

      SUBROUTINE update_taudamp (fsq_save, delt_cg)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(INOUT) :: fsq_save, delt_cg
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: mintau = 0.15_dp
      REAL(dp)            :: ratio
!-----------------------------------------------

      IF (fsq_save .GT. zero) THEN
         ratio = ABS(fsq_total/fsq_save)
         ratio = MIN(mintau, ratio)
         otau(1:ndamp-1) = otau(2:ndamp)
         otau(ndamp) = ratio/delt_cg
         fsq_save = fsq_total
      ELSE
         otau = mintau/delt_cg
      END IF

      END SUBROUTINE update_taudamp

      END MODULE evolution

