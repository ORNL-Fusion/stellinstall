
!>  \brief Module for controlling iteration ("time") evolution of MHD convergence sequence.
!!  Updates the CONTRAVARIANT Fourier components of the "velocity" field (displacement X delt)
      MODULE evolution
      USE stel_kinds
      USE stel_constants
      USE shared_data
      USE shared_functions, ONLY: funct_island, funct_island_par, getwmhd
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,           &
          nfp=>nfp_i, mnmax=>mnmax_i, hs=>hs_i
      USE timer_mod
      USE perturbation, ONLY: eta_factor
      USE descriptor_mod, ONLY: iam, DIAGONALDONE, nprocs,               &
                                in_hess_nfunct, out_hess_nfunct
      USE nscalingtools, ONLY: SKSDBG, TOFU, PARFUNCTISL,                &
                 startglobrow, endglobrow, leftproc, rightproc, MPI_ERR
      USE mpi_inc
      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: n_bsupu_full=1
      INTEGER, PARAMETER, PRIVATE :: ndamp = 10              !<number smoothing steps for conj grad
      INTEGER, PARAMETER, PRIVATE :: ncongrad_steps = 100    !<max number of conj grad steps
      LOGICAL, PARAMETER, PRIVATE :: l_Hess_sym=.FALSE.      !<forces symmetrization of Hessian if TRUE
      LOGICAL, PRIVATE   :: l_EPAR                  !<=TRUE when applying parallel tearing perturbations 
      LOGICAL, PRIVATE   :: l_Gmres                 !<=TRUE performs GMRES iterations 
      LOGICAL, PRIVATE   :: l_conjmin               !<=TRUE if conjugate gradient lowers fsq_min

      INTEGER, PRIVATE   :: nprint_step=1
      INTEGER, PRIVATE   :: nfun_calls
      INTEGER, PRIVATE   :: js, moff, noff

      REAL(dp), PRIVATE    :: fsqprev, fsqprev1
      REAL(dp), PRIVATE    :: fsq_min, fsq_ratio, fsq_ratio1 
      REAL(dp), PRIVATE    :: otau(ndamp)
      REAL(dp), DIMENSION(:), ALLOCATABLE :: xcdot
            
#if defined(PETSC_OPT)
      integer(kind=selected_int_kind(10)) Apc
      integer(kind=selected_int_kind(10)) subApc(1048)
      integer(kind=selected_int_kind(10)) ksp
      integer(kind=selected_int_kind(10)) subksp(1048)
#endif

#if defined(SKS)
      INTEGER, PRIVATE   :: diagnum=1, blocknum=1
      INTEGER            :: EVOLVEPASS=0
#endif
!-----------------------------------------------

      CONTAINS

!>  \brief Performs initial convergence of force residuals using diagonal preconditioner only.
!!  Terminates when the force residual stops changing by more than a factor of 2
      SUBROUTINE converge_diagonal (wout_file, ftol)
      USE perturbation, ONLY: add_perturb, ladd_pert, l_output_alliter   !R. Sanchez, added Jul '10: L_OUTPUT_ALLITER Dump info on file at each timestep
      USE dumping_mod, ONLY: write_output
      USE hessian, ONLY: levmarq_param !MRC
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(IN)  :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(dp), INTENT(IN)    :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL               :: l_iterate            !<call evolve until fsq_ratio stops decreasing
      REAL(dp)              :: fsq_block, t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!
      fsq_block = MAX(ftol, fsq_prec)
      l_iterate = .TRUE.
      nprecon_type = PREDIAG
      IF (nprecon .NE. 0) nprecon_type = PREBLOCK
      niter = 1
      levm_scale = 1    !MRC
!SPH091716      nprecon = 0

      DO WHILE (l_iterate)
         in_hess_nfunct=0; out_hess_nfunct=0
#if defined(SKS)
         CALL second0(t1)
#endif         
         CALL evolve
#if defined(SKS)
         CALL second0(t2)
         diag_evolve_time=diag_evolve_time+(t2-t1)
#endif         
         IF (l_output_alliter) CALL write_output (wout_file, niter, .false.)  ! MRC do not deallocate the woutfile here.

         l_iterate = (fsq_ratio1.LE.0.5_dp .AND. fsq_total1.GT.fsq_block)

         IF (ladd_pert .AND. fsq_total1.LT.100*fsq_block) THEN
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
!!  iterations exceeds a specified value (hard coded for now).
!!  Applies an external island perturbation if not previously done.
      SUBROUTINE converge_blocks (wout_file, ftol)
      USE PERTURBATION, ONLY: niter_max, l_output_alliter, write_restart_file, &     ! R. Sanchez, added Jul '10: Dump info on file at each timestep
                              ladd_pert, add_perturb
      USE dumping_mod, ONLY: write_output
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)  :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(dp), INTENT(in)    :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nrow
      REAL(dp),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL            :: l_iterate          !<controls iteration sequence
      REAL(dp) :: t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!
!      l_conjmin = .TRUE.                     !Set this true to initiate call to conjugate gradient routine
      nprecon_type = PREBLOCK
      niter_max = niter+niter_max
      l_iterate = (niter.LT.niter_max) .AND. (fsq_total1.GT.ftol) 
      nrow = 0

      DO WHILE (l_iterate)
         in_hess_nfunct=0; out_hess_nfunct=0
#if defined(SKS)
         CALL second0(t1)
#endif         
         nprecon = nprecon+1
         CALL evolve
#if defined(SKS)
         CALL second0(t2)
         block_evolve_time=block_evolve_time+(t2-t1)
#endif         
         IF (l_output_alliter) CALL write_output (wout_file, niter, .false.)  ! MRC do not deallocate the woutfile here.

         l_iterate = (niter.LT.niter_max) .AND. (fsq_total1.GT.ftol) 

!STOPPING CRITERIA FOR LACK OF CONVERGENCE PROGRESS

!         IF ((.NOT.l_conjmin .AND. .NOT.l_gmres) .OR.                   &
         IF (nprecon.GT.2 .AND. fsq_ratio.GT.1.E5_dp) THEN
            IF (IAM .EQ. 0 .and. lverbose) PRINT 100, fsq_ratio1
            l_iterate=.FALSE.
         END IF
         IF (nprecon.GT.1 .AND. ABS(1-fsq_ratio1).LT.1.E-2_dp) THEN
            levm_scale = levm_scale/3
            nrow = nrow+1
            IF (nrow.EQ.3 .AND. l_iterate) THEN
               l_iterate=.FALSE.
               IF (IAM .EQ. 0 .and. lverbose) PRINT 100, fsq_ratio1
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

         niter = niter+1
      END DO

100   FORMAT (' FSQ RATIO: ',1p,e10.2, ': SIESTA CAN NOT CONVERGE FURTHER')


      END SUBROUTINE converge_blocks

!>   \brief  Initializes variables and pointers prior to calling evolve
      SUBROUTINE init_evolution
      USE hessian, ONLY: levmarq_param0, levmarq_param, l_Compute_Hessian
      USE quantities, ONLY: jvsupsmncf,jvsupumnsf,jvsupvmnsf,            &
                            jvsupsmnsf,jvsupumncf,jvsupvmncf,            &
                            fsubsmncf, fsubumnsf, fsubvmnsf,             &
                            fsubsmnsf, fsubumncf, fsubvmncf

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat, n1
!-----------------------------------------------

!     "natural" boundary condition (preserves hessian symmetry): push ONLY m=1 u-flux at origin
!     and correctly calculates B^s/r0_B^u for m=1 at origin

      l_natural=.TRUE.
      IF (l_natural) THEN
         l_push_s = .FALSE.
         l_push_u = .FALSE.   !.TRUE.  
         l_push_v = .FALSE.
      ELSE
!     evolve all velocities (m=1 components) at s=0
!     results in J^s = r0 J^u (m=1,s=0) in equilibrium
         l_push_s = .FALSE.   !.TRUE.
         l_push_u = .FALSE.   !.TRUE.
         l_push_v = .FALSE.
      END IF

!      l_push_edge=.TRUE.
      l_push_edge=.FALSE.

      scale_s = hs
!      scale_s = 1        !Default - no scaling
      scale_u = SQRT(scale_s)

      ste = 0
      l_EPAR = .FALSE.

      nfun_calls = 0
      fsqprev = -1;  fsqprev1 = -1
      fsq_min = 1.E20_dp
      niter    = 0
      wtotal0  = -1
      etak_tol = 1.E-01_dp                                                  !initial Newton tolerance parameter
      delta_t = 1
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      l_getwmhd= .FALSE.                                                    !set in getwmhd function
      levmarq_param = levmarq_param0
      fsq_total = 1
      fsq_total1= 1

      IF (mnmax .NE. (1+mpol)*(2*ntor+1)) CALL ErrorAbort('WRONG mnmax VALUE!')

      n1 = ns*mnmax

      IF (lasym) THEN
         ndims = 6
      ELSE
         ndims = 3
      END IF
      neqs = ndims*n1

      ALLOCATE(xc(neqs), col_scale(0:mpol,-ntor:ntor,ndims,ns), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('Allocate xc failed!')
 
      xc = 0
      col_scale = 1
      CALL init_ptrs (xc, jvsupsmncf, jvsupumnsf, jvsupvmnsf,            &
                          jvsupsmnsf, jvsupumncf, jvsupvmncf)

      ALLOCATE(gc(neqs), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('Allocate gc failed!')
      CALL init_ptrs (gc, fsubsmncf, fsubumnsf, fsubvmnsf,               &
                          fsubsmnsf, fsubumncf, fsubvmncf)

      l_Compute_Hessian = .FALSE.
     
      col_scale(:,:,1,:) = scale_s
      col_scale(:,:,2,:) = scale_u
      IF (lasym) THEN
      col_scale(:,:,4,:) = scale_s
      col_scale(:,:,5,:) = scale_u
      END IF

!      levmarq_param0 = 10*levmarq_param0
!      IF (iam .eq. 0) PRINT* ,'LMPARAM0 = ',levmarq_param0

      END SUBROUTINE init_evolution

!>    \brief Initializes pointers for xc (dependent variables) and gc (force) arrays

      SUBROUTINE init_ptrs (xtarget, ptr1, ptr2, ptr3, ptr4, ptr5, ptr6)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), TARGET, INTENT(INOUT)     :: xtarget(0:mpol,-ntor:ntor,ndims,ns) !<input 3D, 2 3-component arrays
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr1, ptr4  !<1st (s) array subsection (lasym=F,T)
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr2, ptr5  !<2nd (u) array subsection
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr3, ptr6  !<3rd (v) array subsection

!-----------------------------------------------
      xtarget = 0
      
      ptr1 => xtarget(:,:,1,:)
      ptr2 => xtarget(:,:,2,:)
      ptr3 => xtarget(:,:,3,:)
      IF (ndims .EQ. 6) THEN
         ptr4 => xtarget(:,:,4,:)
         ptr5 => xtarget(:,:,5,:)
         ptr6 => xtarget(:,:,6,:)
      END IF

      END SUBROUTINE init_ptrs

!>   \brief Evolves one step toward MHD equilibrium
      SUBROUTINE evolve
      USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,          &
                            fsubsmncf, fsubumnsf, fsubvmnsf,             &
                            jbsupsmnsh, jbsupumnch, w_factor, signjac
      USE gmres, ONLY: gmres_fun
      USE hessian
      USE perturbation, ONLY: eta_factor, lresistive, ftol,              &
                              write_restart_file, niter_max
      USE diagnostics_mod, ONLY: bgradp_rms, toroidal_flux0, bgradp,     &
                                 bgradp_par, tflux
      USE siesta_state, ONLY: update_state, update_state_par
      USE siesta_init, ONLY: init_state, init_state_par
#if defined(SKS)
      USE blocktridiagonalsolver_s, ONLY: RefactorHessian
      USE timer_mod, ONLY: gmres_time
#endif
      USE descriptor_mod, ONLY: SIESTA_COMM
!
!     Evolves the equilibrium equations one time step and updates delta_t
!     
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: nm0=1, nn0=2, ns0=3
      REAL(dp), PARAMETER :: levscan(1) = (/ 0.25_dp /) 
      REAL(dp)    :: skston, skstoff
      INTEGER     :: iprint, m1, i
      INTEGER     :: imaxlocv(3), ilbound(3), imaxlocf(3)
      REAL(dp)    :: v2tot, f1, f2, r0
      REAL(dp)    :: vmaxs, vmaxu, vmaxv, fmaxs, fmaxu, fmaxv
      REAL(dp)    :: t1, t2, lm0
	  REAL(dp), SAVE :: levm_loc
      LOGICAL     :: lprint, ldiagonal, lblock
      CHARACTER(18) :: str_resistive
      INTEGER, DIMENSION (2) :: snum_funct, rnum_funct 
!-----------------------------------------------
#if defined(SKS)
      EVOLVEPASS=EVOLVEPASS+1
#endif
      l_init_state = .TRUE.
      IF (niter .LE. 1) levm_loc = 1 ! MRC needed to reset levm_loc

!      l_Gmres = (nprecon_type .EQ. PREDIAG)  ! Turn on CG after initial diagonal preconditioner
!      l_Gmres = l_Gmres .OR. fsq_total1.LT.1.E-9_dp   ! Turn on Gmres after CG reduces residue
      l_Gmres = .TRUE.    !TURNS OFF CONJ-GRAD
!
!     Determine levenberg and mu|| parameter guesses based on ||F||
!
!       IF (nprecon.GE.1 .AND. (fsq_total1.LT.1.E-8_dp .OR. .NOT.l_Gmres)) THEN
      IF (fsq_total1.LT.1.E-6_dp) THEN
!      IF (fsq_total1.LT.1.E-8_dp) THEN
!        PSEUDO-TRANSIENT CONTINUATION: TURN OFF FOR CONSTANT levmarq_param, muPar option
         f2 = fsq_total1*1.E12_dp
!         f2 = 0.5_dp/(one+f2)
         f2 = one/(one+f2)

!ASYM WORKS BETTER LIKE THIS?
         IF (lasym) f2 = f2/4

         f1 = fsq_total1**f2

         f2 = 1
         IF (fsq_ratio1.GT.0.25_dp .AND. fsq_ratio1.LT.1.5_dp) f2 = f2/2
         IF (fsq_ratio1.GT.1.5_dp) f2 = 1.5_dp
         IF (fsq_ratio1.GT.10._dp) f2 = 10
         levm_loc = f2*levm_loc

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
#if defined(SKS)
            CALL Compute_Hessian_Blocks (funct_island_par, l_Hess_sym, ldiagonal)
#else
            CALL Compute_Hessian_Blocks (funct_island, l_Hess_sym, ldiagonal)
#endif
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
#if defined(SKS)
         CALL funct_island_par
#else
         CALL funct_island
#endif
         CALL second0(skstoff)
         evolve_funct_island_time=evolve_funct_island_time+(skstoff-skston)
         l_getwmhd=.FALSE.
      ELSE IF (nprecon .EQ. 1) THEN
         IF (l_Gmres) THEN
           CALL second0(skston)
           CALL gmres_fun
           CALL second0(skstoff)
           gmres_time=gmres_time+(skstoff-skston)
         END IF
      ELSE 
         etak_tol = MIN(0.1_dp, 1.E8_dp*fsq_total)
         IF (fsq_total <= 0) etak_tol = 0.1_dp
         IF (l_Gmres) THEN
           CALL second0(skston)
           CALL gmres_fun
           CALL second0(skstoff)
           gmres_time=gmres_time+(skstoff-skston)
         END IF
      END IF

!      lblock = ((nprecon_type.EQ.PREBLOCK .AND. niter.GT.1) .OR. .NOT.l_Gmres)
!      IF (lblock) THEN 
#if defined(SKS)
      IF (.NOT.l_Gmres .AND. niter.GT.1) THEN
         CALL second0(skston)
         lm0 = levmarq_param
         DO i = 1, 1+SIZE(levscan)
            IF (i .GT. 1) THEN 
               levmarq_param = levscan(i-1)*lm0
               IF (levmarq_param .GT. levmarq_param0) levmarq_param = levmarq_param0*1.E-2_dp
               IF (levmarq_param .LE. levm_ped) levmarq_param = 100*levm_ped
               IF (iam .EQ. 0 .and. lverbose) PRINT 930,' Refactoring Hessian: LM = ',levmarq_param
               CALL RefactorHessian(levmarq_param)
            END IF
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
#endif

      v2tot = SQRT(hs*SUM(xc*xc))
!      IF (v2tot.lt.1.E-5_dp .and. .not.l_Gmres .and. niter.gt.1) THEN
!      END IF
        
      lprint = MOD(niter, nprint_step).EQ.0 .OR. niter.EQ.1
      l_update_state=.TRUE.
#if defined(SKS)
      IF (.NOT.l_par_state) CALL init_state_par(.FALSE.)
      CALL update_state_par(lprint, fsq_total1, zero)
#else
      IF (.NOT.l_par_state) CALL init_state(.FALSE.)
      CALL update_state(lprint, fsq_total1, zero)
#endif
      l_update_state=.FALSE.

!      IF (v2tot.lt.1.E-5_dp .and. .not.l_Gmres .and. niter.gt.1) THEN
!      END IF
        
!
!     Compute force component residues  |Fsubi|**2
!
      IF (fsq_total < zero .and. lverbose) PRINT *, 'Fsq_Total = ', fsq_total,' < 0!'
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
!     p by p_factor, WMHD by w_factor

#if defined(SKS)
      snum_funct(1)=in_hess_nfunct; snum_funct(2)=out_hess_nfunct
      CALL MPI_ALLREDUCE(snum_funct,rnum_funct,2,MPI_INTEGER,MPI_MAX,SIESTA_COMM,MPI_ERR)
      nfun_calls=nfun_calls+(rnum_funct(1)+rnum_funct(2))
#else 
      nfun_calls=nfun_calls+(in_hess_nfunct+out_hess_nfunct)
#endif

     IF (niter .EQ. 1) THEN
#if defined(SKS)
        CALL init_state_par(.TRUE.)
        CALL bgradp_par
#else
        CALL bgradp
#endif
        CALL tflux
     END IF

     IF (lprint) CALL gather_array(gc)                                      !Need to compute maximum forces for printout

     IF (iam.EQ.0 .AND. lprint) THEN

        IF (niter .GT. 1) THEN
          ilbound = LBOUND(jvsupsmncf)-1
          imaxlocv = MAXLOC(ABS(jvsupsmncf)) + ilbound
          vmaxs = jvsupsmncf(imaxlocv(1),imaxlocv(2),imaxlocv(3))
          imaxlocv = imaxlocv - ilbound
          imaxlocf = MAXLOC(ABS(fsubsmncf)) + ilbound
          fmaxs = fsubsmncf(imaxlocf(1),imaxlocf(2),imaxlocf(3))
          IF (abs(fmaxs) .ne. maxval(abs(fsubsmncf))) CALL ErrorAbort('fmaxs WRONG')
          imaxlocf = imaxlocf - ilbound
          IF ((vmaxs.NE.zero .OR. fmaxs.NE.zero) .and. lverbose)                          &
          PRINT 45, ' VS-max: ',vmaxs,imaxlocv(ns0),imaxlocv(nm0)-1,imaxlocv(nn0)-ntor-1, &
                    ' FS-max: ',fmaxs,imaxlocf(ns0),imaxlocf(nm0)-1,imaxlocf(nn0)-ntor-1
          imaxlocv = MAXLOC(ABS(jvsupumnsf)) + ilbound
          vmaxu = jvsupumnsf(imaxlocv(1),imaxlocv(2),imaxlocv(3))
          imaxlocv = imaxlocv - ilbound
          imaxlocf = MAXLOC(ABS(fsubumnsf)) + ilbound
          fmaxu = fsubumnsf(imaxlocf(1),imaxlocf(2),imaxlocf(3))
          IF (abs(fmaxu) .ne. maxval(abs(fsubumnsf))) CALL ErrorAbort('fmaxu WRONG')
          imaxlocf = imaxlocf - ilbound
          IF((vmaxu.NE.zero .OR. fmaxu.NE.zero) .and. lverbose)                           &
          PRINT 45, ' VU-max: ',vmaxu,imaxlocv(ns0),imaxlocv(nm0)-1,imaxlocv(nn0)-ntor-1, &
                    ' FU-max: ',fmaxu,imaxlocf(ns0),imaxlocf(nm0)-1,imaxlocf(nn0)-ntor-1
          imaxlocv = MAXLOC(ABS(jvsupvmnsf)) + ilbound
          vmaxv = jvsupvmnsf(imaxlocv(1),imaxlocv(2),imaxlocv(3))
          imaxlocv = imaxlocv - ilbound
          imaxlocf = MAXLOC(ABS(fsubvmnsf)) + ilbound
          fmaxv = fsubvmnsf(imaxlocf(1),imaxlocf(2),imaxlocf(3))
          IF (abs(fmaxv) .ne. maxval(abs(fsubvmnsf))) CALL ErrorAbort('fmaxv WRONG')
          imaxlocf = imaxlocf - ilbound
          IF ((vmaxv.NE.zero .OR. fmaxv.NE.zero) .and. lverbose)                          &
          PRINT 45, ' VV-max: ',vmaxv,imaxlocv(ns0),imaxlocv(nm0)-1,imaxlocv(nn0)-ntor-1, &
                    ' FV-max: ',fmaxv,imaxlocf(ns0),imaxlocf(nm0)-1,imaxlocf(nn0)-ntor-1
          IF (lverbose)PRINT 110, bsbu_ratio_s, jsju_ratio_s,                         &
                    (bs0(i), i=2,6), (bu0(i), i=2,6) 
          IF (lasym) THEN
          IF (lverbose)PRINT 112, bsbu_ratio_a, jsju_ratio_a,                         &
                    (bs0(i), i=8,12), (bu0(i), i=8,12) 
          END IF
          IF (lverbose)PRINT *," "
        END IF

 45   FORMAT(2(a,1pe10.2,' AT JS: ',i4,' M: ',i4,' N: ',i4,2x))

        DO iprint = 6, 33, 27
           IF (.not.lverbose .and. iprint .eq. 6) CYCLE
           IF (niter .EQ. 1) THEN
              IF (lresistive) THEN
                 str_resistive = "RESISTIVE RUN "
              ELSE
                 str_resistive = "NON-RESISTIVE RUN "
              END IF
              WRITE (iprint, 55) str_resistive, eta_factor, lasym,       &
                                 l_push_s, l_push_u, l_push_v, l_push_edge
              WRITE (iprint, 60) delta_t, etak_tol, l_Hess_sym
              WRITE (iprint, 65) levmarq_param0, mupar0
              WRITE (iprint, 70) neqs, ns, mpol, ntor, nu_i, nv_i,       &
                                 ngmres_steps
           END IF
            
           IF (fsq_lin .EQ. -1) THEN
              WRITE (iprint, 100)                                        &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,   &
                              fsq_total1, fsqvs, fsqvu,                  &
                              fsqvv, v2tot*delta_t, nfun_calls
           ELSE
              WRITE (iprint, 102)                                        &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,   &
                              fsq_total1, fsq_lin, fsqvs, fsqvu,         &
                              fsqvv, v2tot*delta_t, nfun_calls
           END IF
        END DO
     END IF

 50  FORMAT(' WMHD: ', 1pe13.6,' <BETA>: ',1pe11.4,                      &
            ' TFLUX: ',1pe11.4,' B.GRAD-P (rms): ', 1pe11.4,/,21('-'))
 55  FORMAT(1x, a,' ETA_FACTOR: ', 1pe10.2,' LASYM: ',l1,                &
            ' L_PUSH_S: ',l1, ' L_PUSH_U: ',l1,' L_PUSH_V: ',l1,         &
			' L_PUSH_EDGE: ',l1)
 60  FORMAT(' DELTA_T: ',1pe10.2,' ETA_K: ',                             &
             1pe10.2,' HESSIAN SYM: ', l1)
 65  FORMAT(' LEVMARQ_PARAM: ',1pe10.2,' MU_PAR: ',1pe10.2)
 70  FORMAT(' NEQS: ', i6,' NS: ',i4,' MPOL: ',i4,' NTOR: ',i4,          &
             ' NTHETA: ', i4,' NZETA: ',i4, ' NGMRES-STEPS: ', i4,//,    &
             ' NITER (W-W0)/W0*1E6    F2(MHD)    F2(LIN)',               &
             '    F2SUBS     F2SUBU     F2SUBV     |V|rms    NCALLS')
 100 FORMAT(1x,i5,1x, 1pe12.4, 2x, 1x,1pe10.3, 2x,9 ('-'),               &
           4(1x,1pe10.3),3x,i5)
 102 FORMAT(1x,i5,1x, 1pe12.4, 2x, 6(1x,1pe10.3),3x,i5)
 110 FORMAT(' |B^s-r0*B^u|/|B^s+r0*B^u| (m=1,r->0,  sym) : ',1pe10.3,/,  &
            ' |J^s-r0*J^u|/|J^s+r0*J^u| (m=1,r->0,  sym) : ',1pe10.3,/,  &
            ' JBSUPSH(JS=2-6,M=1)/R0  (sym): ',1p5e10.3,/,               &
            ' JBSUPUH(JS=2-6,M=1)          : ',1p5e10.3)
 112 FORMAT(' |B^s+r0*B^u|/|B^s-r0*B^u| (m=1,r->0, asym) : ',1pe10.3,/,  &
            ' |J^s+r0*J^u|/|J^s-r0*J^u| (m=1,r->0, asym) : ',1pe10.3,/,  &
            '-JBSUPSH(JS=2-6,M=1)/R0 (asym): ',1p5e10.3,/,               &
            ' JBSUPUH(JS=2-6,M=1)          : ',1p5e10.3)

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


!>  \brief Parallel routine to evolve state toward equilibrium, using Conjugate Gradients

      SUBROUTINE Conjugate_Grad (ftol)
      USE stel_kinds
      USE hessian, ONLY: mupar
      USE siesta_state, ONLY: update_state, update_state_par
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
      LOGICAL     :: lpar=.FALSE.
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

#if defined(SKS)      
      lpar=.TRUE.  !PARFUNCTISL
      IF (lpar) l_getwmhd=.TRUE.
#endif
 
! 
!     Linearize perturbations around B0 and p0
!     dB = curl(xc X B0)
!     dp = dp(p0,xc)
!
!     NON-LINEAR forces 
!     J = J0+dJ,  B = B0+dB, p = p0+dp

      ALLOCATE(xcdot(neqs), xc0(neqs), stat=icount) 
      IF (icount .NE. 0) CALL ErrorAbort('ALLOCATION FAILED IN CONJ_GRAD')

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
         IF (lpar) THEN
            CALL funct_island_par
         ELSE
            CALL funct_island
         END IF

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

         IF (l_linearize) CALL ErrorAbort('LINEARIZATION TURNED ON!')
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
         xcdot = fac*(b1*xcdot - delt_cg*gc)
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
      IF (lpar) THEN
         CALL funct_island_par
         CALL update_state_par(.FALSE., zero, zero)
      ELSE
         CALL funct_island
         CALL update_state(.FALSE., zero, zero)
      END IF

      l_PrintOriginForces = .FALSE.

      IF (muParS .NE. zero) fsq_min0 = fsq_total1

      IF (iam.EQ.0 .AND. l_conjmin) THEN
         IF (lverbose) WRITE (6, 100) icount, fsq_start, fsq_min0
         WRITE (33,100) icount, fsq_start, fsq_min0
      END IF

      muPar = muParS
      IF (.NOT. l_conjmin) xc = 0
 !     l_Gmres = .NOT.l_conjmin
 !     IF (fsq_min0/fsq_start .GT. 0.95_dp) l_Gmres=.TRUE.

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

