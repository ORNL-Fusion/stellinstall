!!  \author S. P. Hirshman and S. K. Seal
!!  \date Jan, 2014

      MODULE gmres
      USE stel_kinds
      USE stel_constants, ONLY: zero, one
      USE shared_data, etak=>etak_tol
      USE shared_functions, ONLY: funct_island, funct_island_par, LineSearch
      USE siesta_state, ONLY: update_state, update_state_par
      USE nscalingtools, ONLY: PARSOLVER, PARFUNCTISL, MPI_ERR, startglobrow,     &
                               endglobrow
      USE timer_mod
      IMPLICIT NONE

      CONTAINS

!>  \brief Uses gmres to find approximate solution of Ax=b (linearized MHD equations)
!!   that minimize the nonlinear MHD forces

      SUBROUTINE gmres_fun

      USE hessian, ONLY: ns, levmarq_param, Dmin=>mindiag, Apply_Precond, &
                         levmarq_param0, l_diagonal_only, gather_array, maxeigen, muPar
      USE blocktridiagonalsolver, ONLY: RefactorHessian
      USE perturbation, ONLY: lresistive, ftol
      USE descriptor_mod, ONLY: iam
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, iter, istat
      INTEGER :: nsteps
      REAL(dp), PARAMETER :: levscan(1) = (/16._dp/)           !!007C
!      REAL(dp), PARAMETER :: levscan(2) = (/4._dp, 16._dp/)   !!007B
      REAL(dp) :: cntl(5), facmin, res, resid, wmhd, fmhd, gnorm, lm_ped,  &
                  skston, skstoff, eps, lm0, fsq_min, xmax_lin, scale_fac, &
                  muparS
      REAL(dp), ALLOCATABLE :: work(:,:)
      REAL(dp), ALLOCATABLE :: xc1(:), xcmin(:), brhs(:)
      LOGICAL               :: bTerm
!-----------------------------------------------

!     SPH: ONLY works with 0 restart each time
!     EXTERNALLY, l_init_state = .TRUE.   

      xc = 0
      l_linearize = .FALSE.
      l_ApplyPrecon = .FALSE.    !GMRES APPLIES PRECONDITIONER ITSELF
      muParS = muPar
      
!
!     STORE INITIAL UNPRECONDITIONED RESIDUE (BRHS AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      ALLOCATE(brhs(neqs), stat=j)
      IF (j .NE. 0) STOP 'Allocation error in gmres_fun'

      l_getfsq =.TRUE.
      l_init_state=.TRUE.
      fmhd=fsq_total1               !!initial fsq BEFORE application of resistive diffusion

      CALL second0(skston)
      l_getwmhd=.TRUE.
      IF (PARFUNCTISL) THEN
         CALL funct_island_par
         brhs = -gc 
         CALL gather_array(brhs)     
      ELSE
         CALL funct_island
         brhs = -gc
      END IF
      l_getwmhd=.FALSE.

      gnorm = SQRT(SUM(brhs*brhs))
      IF (gnorm .EQ. 0) THEN
         IF (iam .EQ. 0) PRINT *,' GNORM = 0 in GMRES!'
         RETURN
      END IF
      brhs = brhs/gnorm

#if defined(SKS)
      CALL second0(skstoff)
      gmres_funct_island_time=gmres_funct_island_time+(skstoff-skston)
#endif

      wmhd = wtotal
      fsq_min = HUGE(fsq_min)   !Capture one xc at least

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FAST CONVERGENCE - SPH 083016 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(xc1(neqs), xcmin(neqs), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in GMRES'
      xcmin = 0

      iter = 0
      lm0 = levmarq_param

      scale_fac = 1
      IF (iam.EQ.0 .AND. .NOT.l_Diagonal_Only) PRINT '(a,1pe12.3)',' DMIN: ', dmin

920   CONTINUE
      IF (.NOT.l_Diagonal_Only .AND. lcolscale) THEN

      IF (iam .EQ. 0) PRINT 925
925   FORMAT(1x,'FAST CONVERGENCE SUMMARY',/,' -------------',/,1x,'ITER',7x,'FSQ_NL',   &
             10x,'||X||',9x,'MAX|X|')

      xc = 0; l_getfsq=.TRUE.
      l_linearize = .FALSE.
      IF (PARFUNCTISL) THEN
         CALL funct_island_par
      ELSE
         CALL funct_island
      END IF
      xc = -gc
      eps = (0 - levmarq_param)*Dmin
      CALL Apply_Precond(xc)
      xc1 = xc
      DO j = 1, 10
         IF (j .GT. 1) THEN
            gc = -xc1
            CALL Apply_Precond(gc)
            xc1 = eps*gc            !nth partial sum
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
            scale_fac = levmarq_param/lm0
         END IF
         IF (iam .EQ. 0) PRINT 912, j, fsq_total1, SQRT(SUM(xc*xc)), MAXVAL(ABS(xc))
         IF (fsq_min.LE.ftol .OR. bTerm) EXIT
      END DO
      ENDIF

      IF (fsq_min.GT.ftol .AND. fsq_min.GT.fmhd/10) THEN

!!!START GMRES CONVERGENCE
!     RHS, brhs = -gc/gnorm
      gc = brhs
!
!     THIS IS STANDARD GMRES CALL
!
#if defined (SKS)
      CALL second0(skston)
#endif
      CALL gmres_wrap (matvec, xc, gc, gnorm)
#if defined (SKS)
      CALL second0(skstoff)
      gmres_wrap_time=gmres_wrap_time+(skstoff-skston)
#endif

      IF (fsq_total1 .LT. fsq_min) THEN
         xcmin = xc*gnorm
         fsq_min = fsq_total1
         scale_fac = levmarq_param/lm0
      END IF
!!!END GMRES CONVERGENCE

      IF (.NOT.l_Diagonal_Only) THEN
         IF (fsq_min.GE.0.95_dp*fmhd .AND. iter.LT.SIZE(levscan)) THEN
            iter = iter+1
            levmarq_param = lm0*levscan(iter)
            IF (levmarq_param .LT. levm_ped) levmarq_param = 100*levm_ped
            IF (levmarq_param .GE. levmarq_param0) levmarq_param = levmarq_param0/10._dp**iter
            IF (iam .EQ. 0) PRINT 930,' Refactoring Hessian: LM = ',levmarq_param
            CALL RefactorHessian(levmarq_param)
            GOTO 920
         END IF

      END IF

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
      CALL second0(skston)
      IF (PARFUNCTISL) THEN
         CALL funct_island_par
         CALL update_state_par(.FALSE., zero, zero)
      ELSE
         CALL funct_island
         CALL update_state(.FALSE., zero, zero)
      END IF
      
#if defined (SKS)
      CALL second0(skstoff)
      gmres_funct_island_time=gmres_funct_island_time+(skstoff-skston)
#endif
      l_PrintOriginForces=.FALSE.

      muPar = muParS
!      gc = 0  !Need this to print out max forces

      DEALLOCATE (xcmin, brhs)

      END SUBROUTINE gmres_fun


      SUBROUTINE init_gmres0 (icntl, cntl, etak, m, n)
      IMPLICIT NONE
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


      SUBROUTINE gmres_wrap_par (n, m, icntl, cntl, x0, b, gnorm, info)
      USE stel_kinds, ONLY: dp
#if defined(SKS)
      USE stel_constants, ONLY: one, zero
      USE hessian, ONLY: apply_precond, gather_array, ns
      USE perturbation, ONLY: ftol
      USE nscalingtools, ONLY:  rcounts, disp
      USE descriptor_mod, ONLY: iam, nprocs
      IMPLICIT NONE
      INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: n, m
      INTEGER                 :: info(3), icntl(9)
      REAL(dp), INTENT(IN)    :: b(n), gnorm
      REAL(dp), INTENT(INOUT) :: x0(n)
      REAL(dp), INTENT(IN)    :: cntl(5)
#if defined(SKS)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: matveci=1, precondLeft=2,                    &
                            precondRight=3, dotProd=4, peek=5
      INTEGER  :: revcom, colx, coly, colz, nbscal
      INTEGER  :: irc(5), jcount, icount, icntl_s(9)
      INTEGER  :: nout, lwork, kout, kprint
      INTEGER, ALLOCATABLE :: itest(:)
      REAL(dp) :: rinfo(2), fsq_nl, fsq_min, fsq_last, delfsq, fsq_lin, fsq_min_lin
      REAL(dp), ALLOCATABLE :: work(:), xmin(:)
      INTEGER  :: nloc, myrowstart, myrowend
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: aux, tmpbuf
      REAL(dp) :: skston, skstoff, xmax, xmod
      LOGICAL, PARAMETER :: lfast=.FALSE.
!-----------------------------------------------
!
!     EASY-TO-USE WRAPPER FOR GMRES PARALLEL DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine

      fsq_min = -1; fsq_nl = fsq_total1; fsq_min_lin = -1
      delfsq = 1
      jcount = 0; icount = 0; kout = 0; kprint = 0

      nloc=(endglobrow-startglobrow+1)*mblk_size
      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1

      ALLOCATE(tmpbuf(n), aux(n), itest(nprocs), stat=nout)
      IF (nout .NE. 0) STOP 'Allocation error in gmres_wrap_par!'

      lwork = m**2 + m*(nloc+6) + 6*nloc + 1    !Additional space for peek revcom (5 -> 6)
      ALLOCATE (work(lwork), stat=nout)
      IF (nout .NE. 0) STOP 'Allocation error in gmres_wrap_par!'
      work = 0

      IF (icntl(6) .EQ. 1) THEN
         work(1:nloc)=x0(myrowstart:myrowend)
      END IF

      work(nloc+1:2*nloc) = b(myrowstart:myrowend)

      IF (iam .EQ. 0) PRINT 900
900   FORMAT(1x,'GMRES CONVERGENCE SUMMARY',/,' -------------',/,1x,'ITER',7x,'FSQ_NL',   &
             10x,'||X||',9x,'MAX|X|',9x,'FSQ_ARN')


      !****************************************
      !* Reverse communication implementation
      !****************************************

 10   CONTINUE

      CALL second0(skston)
      CALL drive_dgmres(n,nloc,m,lwork,work,irc,icntl,cntl,info,rinfo)
      CALL second0(skstoff)
      drive_dgmres_time=drive_dgmres_time+(skstoff-skston)

      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)


      IF (revcom .EQ. matveci) THEN
!!!        IF (IAM.EQ.0) PRINT *,' CALL MATVECI'
        ! perform the matrix vector product work(colz) <-- A * work(colx)
        CALL ParyAx (work(colx), work(colz), nloc)
        GOTO 10

      ELSE IF (revcom .EQ. precondLeft) THEN
        ! perform the left preconditioning work(colz) <-- M^{-1} * work(colx)
        CALL dcopy(nloc,work(colx),1,work(colz),1)
        GOTO 10

      ELSE IF (revcom .EQ. precondRight) THEN
!!!        IF (IAM.EQ.0) PRINT *,' CALL PRECONDR'        
        ! perform the right preconditioning
        CALL second0(skston)
        CALL dcopy(nloc,work(colx),1,work(colz),1)
        CALL second0(skstoff)
        dcopy_time=dcopy_time+(skstoff-skston)

        CALL second0(skston)
        aux(myrowstart:myrowend) = work(colz:colz+nloc-1)
        CALL second0(skstoff)
        gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)
        CALL apply_precond(aux)
        work(colz:colz+nloc-1)=aux(myrowstart:myrowend)
        GOTO 10

      ELSE IF (revcom .EQ. dotProd) THEN
        ! perform the scalar product (uses nbscal columns of A starting at work(colx))
        ! work(colz) <-- work(colx) work(coly)

        !MAKE SURE nbscal is the same on all processors -
        CALL MPI_ALLGATHER(nbscal,1,MPI_INTEGER,itest,1,MPI_INTEGER,MPI_COMM_WORLD,MPI_ERR)
        IF (IAM.EQ.0 .AND. ANY(itest(1:nprocs) .NE. nbscal)) THEN
           PRINT *,'itest: ',itest(1:nprocs)
           STOP 'nbscal not same!'
        END IF

        CALL second0(skston)
!       THIS IS FASTER THAN ALLGATHERV LOOP, BUT MAY LEAD TO SLIGHTLY DIFFERENT CONVERGENCE
!       SEQUENCES FOR DIFFERENT # PROCESSORS. THE DGEMV CALL IS SLIGHTLY SLOWER THAN THE LOOP
        IF (lFast) THEN

!!        CALL dgemv('C',nloc,nbscal,one,work(colx),nloc,work(coly),1,zero,aux,1)

        DO nout=1,nbscal
           aux(nout) = SUM(work(colx:colx+nloc-1)*work(coly:coly+nloc-1))
           colx = colx+nloc
        END DO
        CALL second0(skstoff)
        dgemv_time=dgemv_time+(skstoff-skston)
        CALL second0(skston)

        CALL MPI_ALLREDUCE(aux,work(colz),nbscal,MPI_REAL8,MPI_SUM,     &
                           MPI_COMM_WORLD,MPI_ERR)
        CALL second0(skstoff)
        gmres_wrap_allreduce_time=gmres_wrap_allreduce_time+(skstoff-skston)

        ELSE


!       THIS DOES NOT SCALE AS WELL WITH # PROCESSORS, BUT LEADS TO A CONVERGENCE SEQUENCE
!       THAT IS INDEPENDENT OF # PROCESSORS
!        CALL MPI_ALLGATHERV(work(coly),nloc,MPI_REAL8,aux,rcounts,      &
!                            disp,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
        CALL MPI_GATHERV(work(coly),nloc,MPI_REAL8,aux,rcounts,         &
                         disp,MPI_REAL8,0,MPI_COMM_WORLD,MPI_ERR)

        DO nout=0,nbscal-1
!           CALL MPI_ALLGATHERV(work(colx),nloc,MPI_REAL8,tmpbuf,rcounts, &
!                               disp,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
           CALL MPI_GATHERV(work(colx),nloc,MPI_REAL8,tmpbuf,rcounts,    &
                            disp,MPI_REAL8,0,MPI_COMM_WORLD,MPI_ERR)
           IF (iam .EQ. 0) work(colz+nout) = SUM(tmpbuf*aux)              !DO FOR ALL PROCS IF ALLGATHER FORM USED
           colx = colx+nloc
        END DO

        CALL MPI_BCAST(work(colz), nbscal, MPI_REAL8, 0, MPI_COMM_WORLD, MPI_ERR)

        END IF

        DO nout = 0,nbscal-1
           CALL Truncate(work(colz+nout), 15)
        END DO

        GOTO 10

      ELSE IF (revcom .EQ. peek) THEN
        fsq_last = fsq_nl                                                 !need for delfsq criteria
        CALL second0(skston)
        aux(myrowstart:myrowend) = work(colx:colx+nloc-1)
        CALL second0(skstoff)
        gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)

        CALL second0(skston)
        CALL GetNLForce(aux, fsq_nl, gnorm)                               !get nonlinear force
        CALL second0(skstoff)
        getnlforce_time=getnlforce_time+(skstoff-skston)

        fsq_lin = (rinfo(1)*gnorm)**2
        icount = icount+1
        IF (fsq_lin .LT. 1.E-30_dp) icntl(7)=info(1)
        IF (fsq_min .EQ. -1) fsq_min = fsq_nl/0.89_dp
!        IF (ngmres_type .EQ. 1) GOTO 10                                !OLDSTYLE:IGNORE LOGIC BELOW
        delfsq = (fsq_last-fsq_nl)/fsq_min                              
        IF (delfsq .LT. 0.05_dp) THEN
          jcount = jcount+1
        ELSE
          jcount = 0
        END IF

        IF (fsq_nl .LT. fsq_min) THEN
          kout = info(1); xmod = gnorm*SQRT(SUM(aux*aux))
          xmax = gnorm*MAXVAL(ABS(aux))
          IF (iam.EQ.0 .AND. (icount.EQ.1 .OR. MOD(icount,10).EQ.0)) THEN
             kprint = kout
             PRINT 910, kout, fsq_nl, xmod, xmax, fsq_lin
          END IF
          IF (fsq_nl .LE. ftol) icntl(7)=info(1)
          IF (.NOT.ALLOCATED(xmin)) ALLOCATE (xmin(nloc))
          xmin = work(colx:colx+nloc-1)
          fsq_min = fsq_nl; fsq_min_lin = fsq_lin
        ELSE IF (fsq_nl .GT. (3*fsq_min) .OR. jcount.GT.3) THEN
          ! STOPPING CRITERIA (reset max iterations to current iteration)         
          icntl(7)=info(1)
        END IF
        GOTO 10
      ENDIF

910  FORMAT(i5, 4(3x,1pe12.3))

!*******************************
! end reverse loop: dump the solution to a file for debugging
!******************************
      GOTO 100

      nout = 11
      OPEN(nout,FILE='sol_dTestgmres',STATUS='unknown')
      IF (icntl(5).eq.0) then
        WRITE(nout,*) 'Orthogonalisation : MGS'
      ELSEIF (icntl(5).eq.1) then
        WRITE(nout,*) 'Orthogonalisation : IMGS'
      ELSEIF (icntl(5).eq.2) then
        WRITE(nout,*) 'Orthogonalisation : CGS'
      ELSEIF (icntl(5).eq.3) then
        WRITE(nout,*) 'Orthogonalisation : ICGS'
      ENDIF
      WRITE(nout,*) 'Restart : ', m
      WRITE(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
      WRITE(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
      WRITE(nout,*) 'Optimal workspace = ', info(3)
      WRITE(nout,*) 'Solution : '
      DO jcount=1,n
        WRITE(nout,*) work(jcount)
      ENDDO
      WRITE(nout,*) '   '

 100  CONTINUE

      IF (ALLOCATED (xmin)) THEN
        CALL second0(skston)
        CALL gather_array(xmin, tmpbuf)
        CALL second0(skstoff)
        gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)
        DEALLOCATE (xmin)
      ELSE
        CALL second0(skston)
        CALL gather_array(work, tmpbuf)
        CALL second0(skstoff)
        gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)
      END IF

      fsq_total1 = fsq_min
      x0 = tmpbuf

      IF (iam.EQ.0 .AND. kprint.NE.kout) &
          PRINT 910, kout, fsq_min, xmod, xmax, fsq_min_lin

      DEALLOCATE (work, tmpbuf, aux, itest)

!      IF (iam .EQ. 0) PRINT 110, icount
! 110  FORMAT(1x,'Total GMRES Iterations: ',i4)
#endif 
     END SUBROUTINE gmres_wrap_par

#if defined (SKS)

      SUBROUTINE ParyAx (p, locAp, nloc)
      USE stel_kinds
      USE hessian, ONLY: ns, gather_array
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                    :: nloc
      REAL(dp), INTENT(OUT), DIMENSION(nloc) :: locAp
      REAL(dp), INTENT(IN),  DIMENSION(nloc) :: p
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
      REAL(dp), ALLOCATABLE,  DIMENSION(:)   :: fullvec
      REAL(dp)                               :: skston, skstoff
!-----------------------------------------------
   
      CALL second0(skston)
      ALLOCATE (fullvec(ns*mblk_size), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in ParyAx'

      CALL gather_array(p, fullvec)

      CALL second0(skstoff)
      gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)

      CALL second0(skston)
      CALL matvec_par(fullvec,locAp,ns*mblk_size,nloc)
      CALL second0(skstoff)

      ParyAx_time=ParyAx_time+(skstoff-skston)

      DEALLOCATE (fullvec)

      END SUBROUTINE ParyAx


      SUBROUTINE matvec_par (p, Ap, ndim, nloc)
      USE stel_kinds
      USE stel_constants, ONLY: zero
      USE hessian, ONLY: eps_factor, mupar
      USE blocktridiagonalsolver, ONLY: ParMatVec
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                    :: ndim, nloc
      REAL(dp), INTENT(IN), DIMENSION(ndim)  :: p
      REAL(dp), INTENT(OUT), DIMENSION(nloc) :: Ap
!-----------------------------------------------
      LOGICAL, PARAMETER :: LPARMATVEC=.FALSE.
      INTEGER            :: myrowstart, myrowend
      REAL(dp)           :: delta
!-----------------------------------------------
!SPH NOTE 032713: EVENTUALLY REMOVE CALL TO ParMatVec IF THE FUNCT_ISLAND_PAR IS FASTER
!
!     NOTE: DO NOT CALL ParMatVec WHEN nprecon_type==PREDIAG
!           SINCE IN HESSIAN, WE SET ALL THE OFF-DIAGONAL COMPONENTS OF 
!           THE A-MATRIX  (USED IN ParMatVec) TO ZERO
!
!           AT PRESENT, ParMatVec WILL NOT WORK IF mupar != 0
!

      IF (nloc .NE. (endglobrow-startglobrow+1)*mblk_size)              &
         STOP 'nloc wrong in matvec_par'
      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1

      IF (.NOT.LPARMATVEC .OR. nprecon_type.EQ.PREDIAG                  &
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

      END SUBROUTINE matvec_par
#endif


      SUBROUTINE gmres_wrap (yAx, x0, b, gnorm)
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: one, zero
      USE perturbation, ONLY: ftol
      USE descriptor_mod, ONLY: iam
      USE nscalingtools, ONLY: PARSOLVER, PARGMRES
      USE hessian, ONLY: apply_precond
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)    :: b(:), gnorm
      REAL(dp), INTENT(INOUT) :: x0(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: matveci=1, precondLeft=2,                    &
                            precondRight=3, dotProd=4, peek=5
      INTEGER :: m, n, icntl(9), info(3)
      INTEGER :: revcom, colx, coly, colz, nbscal
      INTEGER :: irc(5), jcount, icount
      INTEGER :: nout, lwork
      REAL(dp)  :: rinfo(2), cntl(5), fsq_nl, fsq_min, fsq_last, delfsq, &
                   fsq_lin
      REAL(dp), ALLOCATABLE :: work(:), xmin(:)
!-----------------------------------------------
      EXTERNAL yAx
!-----------------------------------------------

      x0  = 0 
      l_linearize = .TRUE.
      l_getfsq = .FALSE.         !So call to funct in matvec does not include gc0
!******************************************************
!*  Initialize GMRES control parameters
!*******************************************************
      CALL init_gmres0(icntl, cntl, etak, m, n)

!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine

      IF (PARSOLVER .AND. PARGMRES) THEN
         CALL gmres_wrap_par (n, m, icntl, cntl, x0, b, gnorm, info)
         RETURN
      END IF

      fsq_min = -1; fsq_nl = fsq_total1
      delfsq = 1
      jcount = 0; icount = 0

      lwork = m**2 + m*(n+6) + 6*n + 1    !Additional space for peek revcom (5 -> 6)
      ALLOCATE (work(lwork), stat=nout)
      IF (nout .NE. 0) STOP 'Allocation error in gmres!'
      work = 0
      IF (icntl(6) .EQ. 1) work(1:n) = x0
      work(n+1:2*n) = b(1:n)

      PRINT 900
900   FORMAT(1x,'GMRES CONVERGENCE SUMMARY',/,' -------------',/,1x,'ITER',7x,'FSQ_NL',   &
             10x,'||X||',9x,'MAX|X|',9x,'FSQ_ARN')

      !****************************************
      !  Reverse communication implementation
      !****************************************

 10   CONTINUE

      CALL drive_dgmres(n,n,m,lwork,work,irc,icntl,cntl,info,rinfo)

      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)

      IF (revcom .EQ. matveci) THEN

        ! perform the matrix vector product work(colz) <-- A * work(colx)
        CALL yAx (work(colx), work(colz), n)
        GOTO 10

      ELSE IF (revcom .EQ. precondLeft) THEN
        
        ! perform the left preconditioning work(colz) <-- M^{-1} * work(colx)
        CALL dcopy(n,work(colx),1,work(colz),1)
        GOTO 10

      ELSE IF (revcom .EQ. precondRight) THEN
        
        ! perform the right preconditioning

        CALL dcopy(n,work(colx),1,work(colz),1)
        CALL apply_precond(work(colz))
        GOTO 10

      ELSE IF (revcom .EQ. dotProd) THEN
        
        ! perform the scalar product (uses nbscal columns of A starting at work(colx))
        ! work(colz) <-- work(colx) work(coly)
!        CALL dgemv('C',n,nbscal,one,work(colx),n,work(coly),1,zero,work(colz),1)
        DO nout = 0, nbscal-1
           work(colz+nout) = SUM(work(colx:colx+n-1)*work(coly:coly+n-1))
           colx = colx+n
        END DO

        DO nout = 0, nbscal-1
          CALL truncate(work(colz+nout),15)
        END DO
        GOTO 10

      ELSE IF (revcom .EQ. peek) THEN

        fsq_last = fsq_nl                                              !need for delfsq criteria
        CALL GetNLForce(work(colx), fsq_nl, gnorm)                     !get nonlinear force
        fsq_lin = (rinfo(1)*gnorm)**2
        icount = icount+1
        IF (fsq_lin .LT. 1.E-30_dp) icntl(7)=info(1)
        IF (fsq_min .EQ. -1) fsq_min = fsq_nl/0.89_dp
        IF (ngmres_type .EQ. 1) GOTO 10                                !OLDSTYLE:IGNORE LOGIC BELOW
        delfsq = (fsq_last-fsq_nl)/fsq_min                             !SPH053014: remove ABS()
        IF (delfsq .LT. 0.05_dp) THEN
          jcount = jcount+1
        ELSE
          jcount = 0
        END IF

        IF (fsq_nl .LT. fsq_min) THEN
          IF (iam .EQ. 0)    &
          PRINT 910, info(1), fsq_nl, gnorm*SQRT(SUM(work(colx:colx+n-1)**2)), &
                                       gnorm*MAXVAL(ABS(work(colx:colx+n-1))), fsq_lin
          IF (fsq_nl .LE. ftol) icntl(7)=info(1)
          IF (.NOT.ALLOCATED(xmin)) ALLOCATE (xmin(n))  !XXX:SKS
          xmin(1:n) = work(colx:colx+n-1)               !XXX:SKS
          fsq_min = fsq_nl
        ELSE IF (fsq_nl.GT.(3*fsq_min) .OR. jcount.GT.3) THEN
          !        STOPPING CRITERIA (reset max iterations to current iteration)         
          icntl(7)=info(1)
        END IF
        GOTO 10

      ENDIF

910  FORMAT(i5, 4(3x,1pe12.3))
!*******************************
! end reverse loop: dump the solution to a file for debugging
!******************************
      GOTO 100

      nout = 11
      OPEN(nout,FILE='sol_dTestgmres',STATUS='unknown')
      IF (icntl(5) .EQ. 0) then
        WRITE(nout,*) 'Orthogonalisation : MGS'
      ELSEIF (icntl(5) .EQ. 1) then
        WRITE(nout,*) 'Orthogonalisation : IMGS'
      ELSEIF (icntl(5) .EQ. 2) then
        WRITE(nout,*) 'Orthogonalisation : CGS'
      ELSEIF (icntl(5) .EQ. 3) then
        WRITE(nout,*) 'Orthogonalisation : ICGS'
      ENDIF
      WRITE(nout,*) 'Restart : ', m
      WRITE(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
      WRITE(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
      WRITE(nout,*) 'Optimal workspace = ', info(3)
      WRITE(nout,*) 'Solution : '
      DO jcount=1,n
        WRITE(nout,*) work(jcount)
      ENDDO
      WRITE(nout,*) '   '

 100  CONTINUE

      IF (ALLOCATED (xmin)) THEN
        x0(1:n) = xmin(1:n)
        DEALLOCATE (xmin)
      ELSE
        x0(1:n) = work(1:n)
      END IF

      fsq_total1 = fsq_min
      DEALLOCATE (work)

!      IF (iam .EQ. 0) PRINT 110, icount
! 110  FORMAT(1x,'Total GMRES Iterations: ',i4)

      END SUBROUTINE gmres_wrap


      SUBROUTINE matvec (p, Ap, ndim)
      USE stel_kinds
      USE hessian, ONLY: ns, eps_factor, gather_array
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE "mpif.h"
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
      IF (ANY(xc0 .NE. zero)) STOP 'xc0 != 0'
      IF (SIZE(gc) .NE. SIZE(Ap)) STOP 'gc and Ap wrong sizes'
      xc = delta*p

      IF (PARFUNCTISL) THEN
#if defined (SKS)
      l_init_state=.NOT.l_par_state
      l_getwmhd=.TRUE.
      CALL funct_island_par
      l_getwmhd=.FALSE.
      ALLOCATE (gc1(ns*mblk_size))
      nloc=(endglobrow-startglobrow+1)*mblk_size
      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1
      CALL gather_array(gc(myrowstart:myrowend), gc1)
      gc = gc1
      DEALLOCATE (gc1)
#endif
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
      IMPLICIT NONE
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
 

      SUBROUTINE GetNLForce(xcstate, fsq_nl, gnorm)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), gnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!     REAL(dp) :: ton, toff, time_nl=0
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
      xc = gnorm*xcstate     !undo internal gmres normalization
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
      IMPLICIT NONE
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
