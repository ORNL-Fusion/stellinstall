!>  \brief Module contained subroutines for initializing the unperturbed
!>  parts of the magnetic fields, pressure and currents in real space
!>  \author S. P. Hirshman and R. Sanchez as part of the ORNL SIESTA (c) PROJECT
!>  \date Jan 12, 2007
       MODULE siesta_init
       USE v3_utilities, ONLY: assert
       USE stel_kinds
       USE descriptor_mod, ONLY: SIESTA_COMM, iam
       USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_sin, f_cos, &
                          f_du, f_dv, f_sum
       USE shared_data, ONLY: l_update_state, l_par_state
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       USE utilities, ONLY: GradientFull, to_full_mesh
       USE metrics, ONLY: tolowerf, tolowerh
       USE quantities
       USE bhtobf
       USE mpi_inc
 
       IMPLICIT NONE
       INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
       INTEGER, PRIVATE            :: nsmin, nsmax
       REAL(dp), PRIVATE           :: p00
       REAL(dp), ALLOCATABLE, PRIVATE :: pfilter(:,:,:)

       CONTAINS

       SUBROUTINE init_state(lcurrent_only, lpar_in)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
       USE shared_data, ONLY: fsq_total, l_getfsq, ste, buv_res
       USE hessian, ONLY: muPar
       USE siesta_currents, ONLY: cv_currents
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)           :: lcurrent_only
       LOGICAL, INTENT(IN), OPTIONAL :: lpar_in
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1
       INTEGER            :: js, istat,  m, n
       REAL(dp)           :: ton, toff, stes(2), rtes(2)
       INTEGER            :: i, n1, n2, endr, startr
       LOGICAL            :: ladd_pert, lpar
!-----------------------------------------------
!      NOTE: jbupXmnPh, jpmnPh are initially computed on FULL ns mesh, and are GATHERED to that 
!            FULL ns mesh in update_state

       istat = 0
       lpar = .TRUE.
       IF (PRESENT(lpar_in)) lpar=lpar_in
       IF (.NOT.lpar) THEN
          startr = startglobrow;  endr = endglobrow
          startglobrow = 1;       endglobrow = ns
       END IF
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+2)
       CALL Init_Allocate_Arrays (lpar)

!Calculate unperturbed half mesh jac*bsupX and jac*pressure on [nsmin,nsmax]
!Sum symmetric and (for lasym) asymmetric parts
!Returning from update_state, jbsupXmn, jpmn have been gathered onto all processors
       CALL GetFields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, isym)
       IF (lasym)                                                       &
       CALL GetFields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh, iasym)

! Convert jac*bsubXh to bsupXh and jac*pres to presh in real space

       bsupsijh0 = bsupsijh0/jacobh(:,:,nsmin:nsmax)
       bsupvijh0 = bsupvijh0/jacobh(:,:,nsmin:nsmax)
       pijh0 = pijh0/jacobh(:,:,nsmin:nsmax)

!      On exit, bsupXF valid on [nsmin, nsmax-1]
       CALL bhalftobfull (bsupsijh0, bsupuijh0, bsupvijh0,              &
                          bsupsijf0, bsupuijf0, bsupvijf0,              &
                          pijh0, pijf0)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       CALL cv_currents (bsupsijh0, bsupuijh0, bsupvijh0,               &
                         ksupsijf0, ksupuijf0, ksupvijf0,               &
                         l_getfsq, .TRUE.)

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          nsmax=MIN(endglobrow+1,ns)
          CALL tolowerf(bsupsijf0,bsupuijf0,bsupvijf0,                  & ! Calculate K || B
                        bsubsijf, bsubuijf, bsubvijf,                   &
                        nsmin, nsmax)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
          IF (ladd_pert) THEN
          ksubsijf(:,:,nsmin:nsmax) = bsubsijf(:,:,nsmin:nsmax)
          ksubuijf(:,:,nsmin:nsmax) = bsubuijf(:,:,nsmin:nsmax)
          ksubvijf(:,:,nsmin:nsmax) = bsubvijf(:,:,nsmin:nsmax)
          END IF
          nsmax=MIN(endglobrow+2,ns)
       END IF

     
!Compute spectrally filtered dp/du and dp/dv on half radial grid [nsmin:nsmax]
       IF (l_update_state) ALLOCATE (pfilter(ntheta, nzeta, nsmin:nsmax))
       CALL FilterPressure(isym)
       IF (lasym)                                                       &
       CALL FilterPressure(iasym)

!Compute radial pressure derivative pijf0_ds at full mesh points [nsmin:nsmax-1]
       CALL GradientFull(pijf0_ds, pijh0)

!SPH: one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
       IF (nsmin .EQ. 1) THEN
          pijf0_ds(:,:,1) = 2*pijh0(:,:,1)*ohs
          pijh0(:,:,1) = 0 
       END IF
      
!SPH11-3-16 s=1 is iso-pressure contour
       IF (nsmax .EQ. ns) pijf0(:,:,ns) = p00

!Fourier filter (Nyquist) check
       IF (.NOT.l_update_state) RETURN
!pijh is unfiltered
       n1=MAX(1,startglobrow); n2=MIN(ns,endglobrow)
       IF (n1 .EQ. 1) pfilter(:,:,1) = 0
       stes(1) = SUM((pijh0(:,:,n1:n2)-pfilter(:,:,n1:n2))**2)
       stes(2) = SUM((pijh0(:,:,n1:n2)+pfilter(:,:,n1:n2))**2)
#if defined(MPI_OPT)
       CALL MPI_ALLREDUCE(stes, rtes, 2, MPI_REAL8, MPI_SUM,            &
                          SIESTA_COMM, MPI_ERR)
       IF (rtes(2) .NE. zero) ste(1) = SQRT(rtes(1)/rtes(2))
#else
       IF (stes(2) .NE. zero) ste(1) = SQRT(stes(1)/stes(2))
#endif
       DEALLOCATE (pfilter)

       IF (.NOT.lpar) THEN
          startglobrow=startr;  endglobrow=endr 
       END IF

       END SUBROUTINE init_state


!>  \brief subroutine for allocating unperturbed array
       SUBROUTINE Init_Allocate_Arrays(lpar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: lpar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL                 :: lrealloc, lalloc
       INTEGER                 :: istat, n1, n2, n3, nsmin, nsmax,      &
                                  nsmin1, nsmax1
!-----------------------------------------------
       lalloc = ALLOCATED(jvsupsijf) 
       lrealloc = .TRUE.
       istat = 0
       nsmin =MAX(1,startglobrow);   nsmax =MIN(endglobrow+1,ns)
       nsmin1=MAX(1,startglobrow-1); nsmax1=MIN(endglobrow+2,ns)

       IF (lpar) THEN
          n1 = ntheta; n2 = nzeta; n3 = nsmax1-nsmin1+1
          l_par_state=.TRUE.
       ELSE
          n1 = ntheta; n2 = nzeta; n3 = ns
          l_par_state=.FALSE.
       END IF

       IF (lalloc) THEN
          IF (SIZE(bsupsijh0,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (jvsupsijf, jvsupuijf, jvsupvijf,                  &
                      bsupsijf0, bsupuijf0, bsupvijf0,                  &
                      bsupsijh0, bsupuijh0, bsupvijh0,                  &
                      bsubsijf,  bsubuijf,  bsubvijf,                   &
                      ksupsijf0, ksupuijf0, ksupvijf0,                  &
                      ksubsijf,  ksubuijf,  ksubvijf,                   &
                      pijf0, pijh0, pijf0_ds,                           &
                      pijh0_du, pijh0_dv,  stat=istat)
          CALL ASSERT(istat.EQ.0,'Deallocation error #1 in INIT_ALLOCATE_ARRAYS')
          END IF
       END IF

       IF (lrealloc) THEN
          n3 = ns
          ALLOCATE (jvsupsijf(n1,n2,n3), jvsupuijf(n1,n2,n3),           & ! Full mesh quantities (real space)
                    jvsupvijf(n1,n2,n3))                                  ! jacobian*(V_s , V_u , V_v)
          IF (lpar) THEN
          ALLOCATE (ksupsijf0(n1,n2,nsmin1:nsmax1),                     & ! K^s, K^u, K^v
                    ksupuijf0(n1,n2,nsmin1:nsmax1),                     &
                    ksupvijf0(n1,n2,nsmin1:nsmax1),                     &
                    bsupsijf0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (full)
                    bsupuijf0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijf0(n1,n2,nsmin1:nsmax1),                     &
                    bsupsijh0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (half)
                    bsupuijh0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijh0(n1,n2,nsmin1:nsmax1),                     &
                    ksubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    ksubuijf(n1,n2,nsmin1:nsmax),                       & ! K^s, K^u, K^v
                    ksubvijf(n1,n2,nsmin1:nsmax),                       &
                    bsubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    bsubuijf(n1,n2,nsmin1:nsmax),                       & ! B_s, B_u, B_v
                    bsubvijf(n1,n2,nsmin1:nsmax),                       &
                    pijh0(n1,n2,nsmin1:nsmax1),                         &
                    pijh0_du(n1,n2,nsmin1:nsmax1),                      &
                    pijh0_dv(n1,n2,nsmin1:nsmax1),                      &
                    pijf0(n1,n2,nsmin1:nsmax1),                         &
                    pijf0_ds(n1,n2,nsmin1:nsmax1), stat=istat)
          ELSE
          ALLOCATE (ksupsijf0(n1,n2,n3), ksupuijf0(n1,n2,n3),           &
                    ksupvijf0(n1,n2,n3),                                &
                    ksubsijf(n1,n2,n3),                                 &
                    ksubuijf(n1,n2,n3),                                 &
                    ksubvijf(n1,n2,n3),                                 &
                    bsubsijf(n1,n2,n3),                                 &
                    bsubuijf(n1,n2,n3),                                 &
                    bsubvijf(n1,n2,n3),                                 &
                    pijf0(n1,n2,n3), pijh0(n1,n2,n3),                   &
                    pijh0_du(n1,n2,n3), pijh0_dv(n1,n2,n3),             &
                    pijf0_ds(n1,n2,n3), bsupsijf0(n1,n2,n3),            &
                    bsupuijf0(n1,n2,n3),bsupvijf0(n1,n2,n3),            &
                    bsupsijh0(n1,n2,n3),                                &
                    bsupuijh0(n1,n2,n3),bsupvijh0(n1,n2,n3),            &
                    stat=istat)
          END IF
          CALL ASSERT(istat.EQ.0, 'Allocation error #1 in Init_Allocate_Arrays')
          jvsupsijf = 0; jvsupuijf = 0; jvsupvijf = 0                     !Need in add_resistivity loop
       END IF


       lalloc = ALLOCATED(ksupsmnsf) 
       lrealloc = .TRUE.

       IF (lpar) THEN
          n1 = mpol+1; n2 = 2*ntor+1; n3 = nsmax-nsmin1+1
       ELSE
          n1 = mpol+1; n2 = 2*ntor+1; n3 = ns
       END IF

       IF (lalloc) THEN
          IF (SIZE(ksupsmnsf,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (ksupsmnsf, ksupumncf, ksupvmncf,                  &
                      djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,   &
                      stat=istat)
          CALL ASSERT(istat.EQ.0,'Deallocation error #2 in Init_Allocate_Arrays')
          IF (lasym) THEN
          DEALLOCATE (ksupsmncf, ksupumnsf, ksupvmnsf,                  &
                      djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh,   &
                      stat=istat)
          CALL ASSERT(istat.EQ.0,'Deallocation error #3 in Init_Allocate_Arrays')
          END IF
          END IF
       END IF

       IF (lrealloc) THEN
          IF (lpar) THEN
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       & 
                       ksupumncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       ksupvmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       djpmnch(0:mpol,-ntor:ntor,nsmin:nsmax),          &
                       djbsupsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupumnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupvmnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       stat=istat)
             IF (lasym) THEN
             ALLOCATE (ksupsmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       & 
                       ksupumnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       ksupvmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       djpmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),          &
                       djbsupsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       stat=istat)
             END IF
                       
          ELSE
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,n3),                 & 
                       ksupumncf(0:mpol,-ntor:ntor,n3),                 &
                       ksupvmncf(0:mpol,-ntor:ntor,n3),                 & 
                       djpmnch(0:mpol,-ntor:ntor,n3),                   &
                       djbsupsmnsh(0:mpol,-ntor:ntor,n3),               &
                       djbsupumnch(0:mpol,-ntor:ntor,n3),               &
                       djbsupvmnch(0:mpol,-ntor:ntor,n3),               &
                       stat=istat)
             IF (lasym) THEN
             ALLOCATE (ksupsmncf(0:mpol,-ntor:ntor,n3),                 & 
                       ksupumnsf(0:mpol,-ntor:ntor,n3),                 &
                       ksupvmnsf(0:mpol,-ntor:ntor,n3),                 & 
                       djpmnsh(0:mpol,-ntor:ntor,n3),                   &
                       djbsupsmnch(0:mpol,-ntor:ntor,n3),               &
                       djbsupumnsh(0:mpol,-ntor:ntor,n3),               &
                       djbsupvmnsh(0:mpol,-ntor:ntor,n3),               &
                       stat=istat)
             END IF
          END IF
          CALL ASSERT(istat.EQ.0, 'Allocation error #2 in Init_Allocate_Arrays')
          djpmnch = 0; djbsupsmnsh = 0; djbsupumnch = 0; djbsupvmnch = 0 
          IF (lasym) THEN
          djpmnsh = 0; djbsupsmnch = 0; djbsupumnsh = 0; djbsupvmnsh = 0 
          ENDIF
       END IF


       END SUBROUTINE Init_Allocate_Arrays
       

!>  \brief Subroutine for getting magnetic field components, pressure in
!>   real space (half radial grid)
       SUBROUTINE GetFields(jbsupsmnh, jbsupumnh, jbsupvmnh,            &
                            jpmnh, ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                       &
             jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
       INTEGER, INTENT(IN) :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: fours, fouru, fourv, icomb
!-----------------------------------------------
       IF (ipar .EQ. isym) THEN
          icomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
       ELSE
          icomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
       END IF

       CALL toijsp_par(jbsupsmnh(:,:,nsmin:nsmax), bsupsijh0, icomb,    &
                       fours)
       CALL toijsp_par(jbsupumnh(:,:,nsmin:nsmax), bsupuijh0, icomb,    &
                       fouru)
       CALL toijsp_par(jbsupvmnh(:,:,nsmin:nsmax), bsupvijh0, icomb,    &
                       fourv)
       CALL toijsp_par(jpmnh(:,:,nsmin:nsmax), pijh0, icomb,            &
                       fourv)

       END SUBROUTINE GetFields


!>  \brief Subroutine for getting angular derivatives of the pressure
       SUBROUTINE FilterPressure(ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: four, icomb, mfour, nfour, istat
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: pmnh
!-----------------------------------------------
       IF (ipar .EQ. isym) THEN
          four = f_cos; icomb = f_none
          mfour = f_du; nfour = f_dv
       ELSE
          four = f_sin; icomb=f_sum
          mfour = IOR(f_du, f_sum);  nfour = IOR(f_dv, f_sum)
       END IF

!pmnh contains Nyquist-filtered (to MPOL, NTOR) Fourier harmonics of p on half mesh
       ALLOCATE (pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation error in FilterPressure')

       CALL tomnsp_par(pijh0, pmnh, four)
       IF (nsmax.EQ.ns .AND. ipar.EQ.isym) p00 = pmnh(0, 0, ns)

!spectrally filtered angle derivatives of p in real space (half mesh)
       CALL toijsp_par(pmnh, pijh0_du, mfour, four)
       CALL toijsp_par(pmnh, pijh0_dv, nfour, four)

       IF (l_update_state)                                              &
       CALL toijsp_par(pmnh, pfilter, icomb, four)                        ! p (fourier filtered)

       DEALLOCATE (pmnh)

       END SUBROUTINE FilterPressure

!>  \brief Subroutine for checking that p>0
       SUBROUTINE CheckPressure
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp) :: pmin, pming
!-----------------------------------------------
       pmin = MINVAL(pijh0)
       pming = pmin
#if defined(MPI_OPT)
       IF (PARSOLVER) THEN
          CALL MPI_ALLREDUCE(pmin, pming, 1 , MPI_REAL8,                &
               MPI_MIN, SIESTA_COMM, MPI_ERR)
       END IF
#endif
       
       IF (iam.EQ.0 .AND. lverbose .AND. pming.LT.zero) THEN
          PRINT 100, 'pijh < 0 in init_state.  ===>',                   &
              ' Recommend lowering f2 exponent in EVOLUTION'
       END IF
 100   FORMAT(1x,2a)

!SPH101017 This may stabilize convergence (lower f2 exponent or raise lev_marq)
       IF (pmin .LT. zero) WHERE (pijh0 .LT. zero) pijh0 = zero
       
       END SUBROUTINE CheckPressure

       END MODULE siesta_init
