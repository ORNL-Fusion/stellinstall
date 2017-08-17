
!!  \authors S. P. Hirshman and S. K. Seal
!!  \date Jan, 2014

      MODULE gmres
      USE stel_kinds
      USE stel_constants, ONLY: zero, one
      USE shared_data, etak=>etak_tol
      USE shared_functions, ONLY: funct_island, funct_island_par, LineSearch
      USE siesta_state, ONLY: update_state, update_state_par
      USE nscalingtools, ONLY: PARFUNCTISL, MPI_ERR, startglobrow,       &
                               endglobrow, PARSOLVER
      USE timer_mod
      IMPLICIT NONE

      INTEGER, PRIVATE  :: icounter=0

      CONTAINS

!>  \brief Uses gmres to find approximate solution of Ax=b (linearized MHD equations)
!!   that minimize the nonlinear MHD forces

      SUBROUTINE gmres_fun
      USE hessian, ONLY: levmarq_param, Dmin=>mindiag, Apply_Precond,    &
                         levmarq_param0, l_diagonal_only, gather_array,  &
                         maxeigen, muPar
      USE perturbation, ONLY: lresistive, ftol
      USE descriptor_mod, ONLY: iam
      USE mpi_inc
#if defined(SKS)
      USE blocktridiagonalsolver_s, ONLY: RefactorHessian
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, iter, istat
      INTEGER :: nsteps
      REAL(dp), PARAMETER :: levscan(1) = (/16._dp/)           !!007C
!      REAL(dp), PARAMETER :: levscan(2) = (/4._dp, 16._dp/)   !!007B
      REAL(dp) :: cntl(5), facmin, res, resid, wmhd, fmhd, lm_ped,       &
                  skston, skstoff, eps, lm0, fsq_min, xmax_lin,          &
                  muparS, scale_fac
      REAL(dp), ALLOCATABLE :: work(:,:)
      REAL(dp), ALLOCATABLE :: xc1(:), xcmin(:), brhs(:)
      LOGICAL               :: bTerm
!-----------------------------------------------

!     SPH: ONLY works with 0 restart each time
!     EXTERNALLY, l_init_state = .TRUE.   

      xc = 0
      l_linearize = .FALSE.
      l_ApplyPrecon = .FALSE.    !LET GMRES APPLY PRECONDITIONER ITSELF
!      l_ApplyPrecon = .TRUE.     !TEST icntl(4) = 0


      muParS = muPar
      
!
!     STORE INITIAL UNPRECONDITIONED RESIDUE (BRHS=-GC AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      l_getfsq =.TRUE.
      l_init_state=.TRUE.
      fmhd=fsq_total1               !!initial fsq BEFORE application of resistive diffusion

      ALLOCATE (brhs(neqs), stat=istat)
      IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN GMRES_FUN'
      l_getwmhd=.TRUE.


#if defined(SKS)
      CALL second0(skston)
      CALL funct_island_par
      brhs = -gc 
      CALL gather_array(brhs)     
      CALL second0(skstoff)
      gmres_funct_island_time=gmres_funct_island_time+(skstoff-skston)
#else
      CALL funct_island
      brhs = -gc
#endif
      l_getwmhd=.FALSE.

      wmhd = wtotal
      fsq_min = HUGE(fsq_min)   !Capture one xc at least

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FAST CONVERGENCE - SPH 083016 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(xc1(neqs), xcmin(neqs), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in GMRES'
      xcmin = 0

      iter = 0
      lm0 = levmarq_param
      scale_fac = 1

920   CONTINUE
!
!     FAST CONVERGENCE GUESS
!     ITERATIVELY SOLVE (A* - eps*I) x = b, WHERE A* = A + eps*I 
!     WHERE x = x0 + x1 + ..., eps = lev_marq_param*Dmin
!     x0 = P*b   where P = (A*)^-1
!     xn = P*(lm*x[n-1])
!
      IF (.NOT.l_Diagonal_Only .AND. lcolscale) THEN

      IF (iam .EQ. 0 .and. lverbose) PRINT 925, Dmin
925   FORMAT(' DMIN: ',1pe12.3,/,1x,'FAST CONVERGENCE SUMMARY',          &
             /,' -------------',/,1x,'ITER',7x,'FSQ_NL',                 &
             10x,'||X||',9x,'MAX|X|')

      xc = 0; l_getfsq=.TRUE.
      l_linearize = .FALSE.
      IF (PARFUNCTISL) THEN
         CALL funct_island_par
      ELSE
         CALL funct_island
      END IF
      xc = -gc
      eps = levmarq_param*Dmin
      IF (.NOT.l_ApplyPrecon) CALL Apply_Precond(xc)
      xc1 = xc
      DO j = 1, 10
         IF (j .GT. 1) THEN
            gc = xc1
            CALL Apply_Precond(gc)
            xc1 = eps*gc            !nth partial sum: store for next j iteration
            xc = xc + xc1
         END IF
         IF (PARFUNCTISL) THEN
            CALL funct_island_par
         ELSE
            CALL funct_island
         END IF

         IF (fsq_total1 .GT. fsq_min) THEN
            bTerm = (J .GT. 1) 
         ELSE
            xcmin = xc
            bTerm = (fsq_total1 .GE. 0.95_dp*fsq_min)
            fsq_min = fsq_total1
         END IF
         IF (iam .EQ. 0 .and. lverbose) PRINT 912, j, fsq_total1, SQRT(SUM(xc*xc)),     &
                         MAXVAL(ABS(xc))
         IF (fsq_min.LE.ftol .OR. bTerm .OR. eps.EQ.zero) EXIT
      END DO
      ENDIF

      IF (fsq_min.GT.ftol .AND. fsq_min.GT.fmhd/10) THEN

!!!START GMRES CONVERGENCE
!
!     THIS IS STANDARD GMRES CALL
!
#if defined (SKS)
      CALL second0(skston)
#endif
      CALL gmres_wrap (xc, brhs)
#if defined (SKS)
      CALL second0(skstoff)
      gmres_wrap_time=gmres_wrap_time+(skstoff-skston)
#endif

      IF (fsq_total1 .LT. fsq_min) THEN
         xcmin = xc
         fsq_min = fsq_total1
         scale_fac = levmarq_param/lm0
      END IF
!!!END GMRES CONVERGENCE
#if defined(SKS)
!FOR NOW, THIS IS ONLY IMPLEMENTED FOR PARSOLVER=TRUE
      IF (.NOT.l_Diagonal_Only .AND. PARSOLVER) THEN
         IF (fsq_min.GE.0.95_dp*fmhd .AND. iter.LT.SIZE(levscan)) THEN
            iter = iter+1
            levmarq_param = lm0*levscan(iter)
            IF (levmarq_param .LT. levm_ped) levmarq_param = 100*levm_ped
            IF (levmarq_param .GE. levmarq_param0) levmarq_param = levmarq_param0/10._dp**iter
            IF (iam .EQ. 0 .and. lverbose) PRINT 930,' Refactoring Hessian: LM = ',levmarq_param
            CALL RefactorHessian(levmarq_param)
            GOTO 920
         END IF
      END IF
#endif
      END IF

 912  FORMAT(i5,3(3x,1pe12.3))
 930  FORMAT (/,a,1pe12.3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 1000 CONTINUE

      DEALLOCATE(xc1)
      l_getfsq = .TRUE.

 1100 CONTINUE

      xmax_lin = MAXVAL(ABS(xcmin))
      IF (iam .EQ. 0) THEN
!         PRINT 1020,' ||X||-GMRES = ', SQRT(SUM(xcmin*xcmin)), ' MAX(|X|) = ',xmax_lin
         WRITE (33, 1020) ' ||X||-GMRES = ', SQRT(SUM(xcmin*xcmin)), ' MAX(|X|) = ',xmax_lin
      END IF
 1020 FORMAT(2(a,1p,e11.3))

      fsq_lin = fsq_total
      l_linearize = .FALSE.

      IF (fsq_min .GT. ftol) THEN

!     TURN OFF NEXT 2 LINES FOR TESTING f1=levm IN EVOLUTION
      levm_scale = levm_scale*scale_fac  
      levm_scale = MIN(100._dp, levm_scale)      
!
!     SIMPLE LINE SEARCH SCAN OF |F| FOR VARIOUS |X| 
      muPar = 0
      CALL LineSearch(xcmin, fsq_min)

      END IF

      xc = xcmin
!
!     RECOMPUTE PERTURBED FIELDS, PRESSURE FOR MINIMUM STATE
!     THEN UPDATE NONLINEAR STATE BY ADDING XC (PERTURBATION)
!

      l_PrintOriginForces=.TRUE.
      l_init_state=.TRUE.
      l_getwmhd= .TRUE.
#if defined(SKS)
      CALL second0(skston)
      CALL funct_island_par
      CALL update_state_par(.FALSE., zero, zero)
      CALL second0(skstoff)
      gmres_funct_island_time=gmres_funct_island_time+(skstoff-skston)
#else
      CALL funct_island
      CALL update_state(.FALSE., zero, zero)
#endif
      
      l_PrintOriginForces=.FALSE.

      muPar = muParS
!      gc = 0  !Need this to print out max forces

      DEALLOCATE (xcmin, brhs)

      END SUBROUTINE gmres_fun

      SUBROUTINE init_gmres0 (icntl, cntl, etak, m, n)
      INTEGER, INTENT(OUT) :: icntl(9), m, n
      REAL(dp),INTENT(OUT) :: cntl(5)
      REAL(dp), INTENT(IN) :: etak
      REAL(dp)             :: skston, skstoff
      INTEGER, PARAMETER   :: noPrec=0, leftPrec=1, rightPrec=2, dblePrec=3

#if defined (SKS)
      CALL second0(skston)
#endif
      CALL init_dgmres(icntl ,cntl)
#if defined (SKS)
      CALL second0(skstoff)
      gmres_init_dgmres_time=gmres_init_dgmres_time+(skstoff-skston)
#endif

!*************************
!*  Tune some parameters
!*************************
! Tolerance
       cntl(1) = etak
!       cntl(1) = 1.E-5_dp
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      icntl(2) = 21
! Save the convergence history in file fort.20
      IF (nprecon .gt. 1) icntl(3) = 20
! Preconditioning INPUT flags (note: different than revcom flags in driver)
!      icntl(4) = leftPrec
      icntl(4) = rightPrec 
!      icntl(4) = dblePrec
!! ICGS orthogonalization
      icntl(5) = 3
! Initial guess
      icntl(6) = 0
!      icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
      icntl(7) = ngmres_steps

      icntl(8) = 1             !Default
      icntl(9) = 1             !Steps for peek at progress during rev com loop
!*********************************
! Choose the restart parameter
!*********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      n   = neqs
      m = 200

      END SUBROUTINE init_gmres0

      SUBROUTINE gmres_wrap (x0, b)
      USE perturbation, ONLY: ftol
      USE descriptor_mod, ONLY: iam, nprocs, SIESTA_COMM
      USE nscalingtools, ONLY: PARSOLVER, PARGMRES, rcounts, disp
      USE hessian, ONLY: apply_precond
      USE gmres_lib, ONLY: gmres_info, gmres_par, gmres_ser
#if defined(SKS)
      USE mpi_inc
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)    :: b(:)
      REAL(dp), INTENT(INOUT) :: x0(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE(gmres_info)   :: gi
      INTEGER :: n
!-----------------------------------------------
      x0  = 0 
      l_linearize = .TRUE.
      l_getfsq = .FALSE.         !So call to funct in matvec does not include gc0

!******************************************************
!*  Initialize GMRES control parameters and Load GMRES_INFO structure
!*******************************************************
      CALL init_gmres0(gi%icntl, gi%cntl, etak, gi%m, n)
      gi%ftol=ftol
      gi%my_comm = SIESTA_COMM
      gi%my_comm_world = SIESTA_COMM
      gi%lverbose = lverbose

!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine
#if defined(SKS)
      IF (PARSOLVER .AND. PARGMRES) THEN
         gi%endglobrow=endglobrow; gi%startglobrow=startglobrow
         gi%iam=iam; gi%nprocs=nprocs
         gi%mblk_size=mblk_size; gi%rcounts=>rcounts; gi%disp=>disp
        
         CALL gmres_par (n, gi, matvec_par, apply_precond, GetNLForce,   &
                         x0, b)
      ELSE
#endif
         CALL gmres_ser (n, gi, matvec, apply_precond, GetNLForce,       &
                         x0, b)
#if defined(SKS)
      END IF
#endif
      END SUBROUTINE gmres_wrap

#if defined (SKS)
      SUBROUTINE matvec_par (ploc, Ap, nloc)
      USE stel_kinds
      USE stel_constants, ONLY: zero
      USE island_params, ONLY: ns=>ns_i
      USE hessian, ONLY: eps_factor, mupar, gather_array
      USE blocktridiagonalsolver_s, ONLY: ParMatVec
      USE mpi_inc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                    :: nloc
      REAL(dp), INTENT(IN), DIMENSION(nloc)  :: ploc
      REAL(dp), INTENT(OUT), DIMENSION(nloc) :: Ap
!-----------------------------------------------
      LOGICAL, PARAMETER :: LPARMATVEC=.FALSE.
      INTEGER            :: myrowstart, myrowend, istat
      REAL(dp), ALLOCATABLE :: p(:)
      REAL(dp)           :: delta, skston, skstoff
!-----------------------------------------------
!SPH NOTE 032713: EVENTUALLY REMOVE CALL TO ParMatVec IF THE FUNCT_ISLAND_PAR IS FASTER
!
!     NOTE: DO NOT CALL ParMatVec WHEN nprecon_type==PREDIAG
!           SINCE IN HESSIAN, WE SET ALL THE OFF-DIAGONAL COMPONENTS OF 
!           THE A-MATRIX  (USED IN ParMatVec) TO ZERO
!
!           AT PRESENT, ParMatVec WILL NOT WORK IF mupar != 0
!
      IF (nloc .NE. (endglobrow-startglobrow+1)*mblk_size)               &
         STOP 'nloc wrong in matvec_par'

      CALL second0(skston)
      ALLOCATE (p(ns*mblk_size), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in matvec_par'

      CALL gather_array(ploc, p)

      CALL second0(skstoff)
      gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)

      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1

      CALL second0(skston)

      IF (.NOT.LPARMATVEC .OR. nprecon_type.EQ.PREDIAG                   &
          .OR. mupar.NE.zero) THEN
         delta = SQRT(EPSILON(delta))*eps_factor
         xc = delta*p

         l_linearize=.TRUE.
         l_init_state=.NOT.l_par_state
         CALL funct_island_par
         Ap = gc(myrowstart:myrowend)/delta

      ELSE 
         CALL ParMatVec(p,Ap,nloc)
      END IF

      DEALLOCATE (p)

      CALL second0(skstoff)

      ParyAx_time=ParyAx_time+(skstoff-skston)

      END SUBROUTINE matvec_par
#endif

      SUBROUTINE matvec (p, Ap, ndim)
      USE stel_kinds
      USE hessian, ONLY: eps_factor, gather_array
#if defined(MPI_OPT)
      USE mpi_inc
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: ndim
      REAL(dp), INTENT(IN)  :: p(ndim)
      REAL(dp), INTENT(OUT) :: Ap(ndim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER   :: zero=0
      INTEGER               :: nloc, myrowstart, myrowend
      REAL(dp), ALLOCATABLE :: gc1(:)
      REAL(dp)              :: delta
!-----------------------------------------------

      delta = SQRT(EPSILON(delta))*eps_factor

!     Note: GMRES without Preconditioner has pnorm = 1 (checked!)
!      pnorm = SUM(p*p)
!      delta = delta/SQRT(pnorm)
!      delta = delta*(pnorm + SUM(p*xc0))/pnorm

!IF CALLED IN PAR MODE, FIRST GATHER THE xc's
      IF (SIZE(gc) .NE. SIZE(Ap)) STOP 'gc and Ap wrong sizes'
      xc = delta*p

      IF (PARFUNCTISL) THEN
         l_init_state=.NOT.l_par_state
         l_getwmhd=.TRUE.
         CALL funct_island_par   !Note: not necessary to gather gc to all processors!
         l_getwmhd=.FALSE.
      ELSE
         l_init_state=l_par_state
         CALL funct_island
      END IF

      IF (l_linearize) THEN
         Ap = gc/delta
      ELSE
         Ap = (gc-gc0)/delta
      END IF

      END SUBROUTINE matvec


      SUBROUTINE get_etak(fk, fkm, ftol, etak)
      USE stel_kinds, ONLY: rprec, dp
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(in)    :: fk, fkm, ftol
      REAL(rprec), INTENT(inout) :: etak
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: gamma1=0.9_dp, alpha=1.5_dp, eta0 = 1.E-02_dp
      REAL(rprec)            :: etakm
!-----------------------------------------------
!     ROUTINE PROVIDED L. CHACON (11/09/06)
      etakm = etak

!     Superlinear convergence
      etak = gamma1*(fk/fkm)**alpha

!     Safeguard avoid sharp decrease in etak
      etak = MIN(eta0, MAX(etak, gamma1*etakm**alpha))

!     Safeguard avoid "oversolving"
      etak = MIN(eta0, MAX(etak, gamma1*ftol/fk))

      END SUBROUTINE get_etak
 

      SUBROUTINE GetNLForce(xcstate, fsq_nl, bnorm, ndim)
      USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ndim
      REAL(dp),INTENT(IN)  :: xcstate(ndim), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL              :: l_linloc, l_Apploc, l_getloc, l_initloc
!-----------------------------------------------
!      CALL second0(ton)
!Store state 
      l_linloc = l_linearize
      l_Apploc = l_ApplyPrecon
      l_getloc = l_getfsq
      l_initloc = l_init_state

!Set state variables
      xc(1:ndim) = bnorm*xcstate(1:ndim)                                  !undo internal gmres normalization
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      l_init_state = .NOT.l_par_state
      l_ApplyPrecon=.FALSE.

#if defined(MPI_OPT)
      CALL funct_island_par
#else
      CALL funct_island
#endif
      fsq_nl = fsq_total1
 
!Restore state variables
      l_linearize = l_linloc
      l_ApplyPrecon = l_Apploc
      l_getfsq = l_getloc
      l_init_state = l_initloc

      END SUBROUTINE GetNLForce


      SUBROUTINE qmr_fun
      USE stel_kinds
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ndim, nlen, nlim, ierr, info(4), j
      INTEGER :: revcom, colx, colb
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: vecs
      REAL(rprec) :: tol = 1.E-4_dp, gnorm
!-----------------------------------------------
!
!     STORE INITIAL POINT AND INITIAL RESIDUE (AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      ALLOCATE(xc0(neqs), gc0(neqs), vecs(neqs,9), stat=j)
      IF (j .ne. 0) STOP 'Allocation error in qmr_fun'

      ndim = SIZE(vecs,1)
      nlen = ndim
      nlim = 100

!     SPH: ONLY works with 0 restart each time??!!!     
      l_linearize = .FALSE.
      l_ApplyPrecon = .TRUE.

      xc = 0
      CALL funct_island
      xc0 = xc
      gc0 = gc

      l_linearize = .TRUE.
!
!     INITIALIZE vecs
!
      gnorm = SUM(gc*gc)
      gnorm = SQRT(gnorm)
      vecs(:ndim,2) = -gc(:ndim)/gnorm
      vecs(:ndim,3) =  gc(:ndim)/gnorm
    
      ierr = 100000
      info = 0
      info(1) = ierr

 10   CALL dutfx (ndim,nlen,nlim,vecs,tol,info)
      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
      IF (revcom .eq. 1) THEN
         CALL matvec (vecs(1,colx), vecs(1,colb), ndim)
         GO TO 10
      END IF

      xc(1:ndim) = xc0(1:ndim) + gnorm*vecs(:,1)

!
!     RECOMPUTE delta current,pressure prior to update_state call
!
      CALL funct_island
      fsq_lin = fsq_total

      l_linearize = .FALSE.
      CALL funct_island

      DEALLOCATE (xc0, gc0, vecs)

      END SUBROUTINE qmr_fun

      END MODULE gmres
