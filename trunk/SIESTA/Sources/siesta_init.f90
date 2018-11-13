!>  \brief Module contained subroutines for initializing the unperturbed
!>  parts of the magnetic fields, pressure and currents in real space
!>  \author S. P. Hirshman and R. Sanchez as part of the ORNL SIESTA (c) PROJECT
!>  \date Jan 12, 2007
       MODULE siesta_init
       USE v3_utilities, ONLY: assert
       USE stel_kinds
       USE descriptor_mod, ONLY: SIESTA_COMM, iam
       USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_sin, f_cos, &
                          f_du, f_dv, f_sum, orthonorm
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       USE quantities
       USE shared_data, ONLY: l_update_state, l_pedge
       USE mpi_inc
 
       IMPLICIT NONE
       PRIVATE
       INTEGER, PARAMETER :: isym=0, iasym=1
       INTEGER, PARAMETER :: m0=0, m1=1
       INTEGER            :: nsmin, nsmax
       REAL(dp), ALLOCATABLE :: pfilter(:,:,:)
       REAL(dp)           :: pedge
       
       PUBLIC  :: init_state
       

       CONTAINS

       SUBROUTINE init_state(lcurrent_only, lpar_in)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
       USE bhtobf
       USE utilities, ONLY: GradientFull
       USE metrics, ONLY: tolowerf
       USE shared_data, ONLY: fsq_total, l_getfsq, ste, buv_res,        &
                              l_update_state
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
       INTEGER            :: istat, n1, n2, nmax, endr, startr
       REAL(dp)           :: ton, toff, stes(2), rtes(2)
       LOGICAL            :: ladd_pert0, lpar
!-----------------------------------------------
!NOTE: jbupXmnPh, jpmnPh are initially computed on FULL ns mesh, and are GATHERED to that 
!      FULL ns mesh in update_state

       istat = 0
       lpar = .TRUE.
       IF (PRESENT(lpar_in)) lpar=lpar_in
       IF (.NOT.lpar) THEN
          startr = startglobrow;  endr = endglobrow
          startglobrow = 1;       endglobrow = ns
       END IF
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+2)
       
!Unnecessary to recompute magnetic fields or currents: already stored at this point!
       ladd_pert0 = ALLOCATED(buv_res)
       IF (ladd_pert0) GOTO 100
       
       CALL Init_Allocate_Arrays (lpar)


!Calculate unperturbed half mesh jac*bsupX and jac*pressure on [nsmin,nsmax]
!Sum symmetric and (for lasym) asymmetric parts
!Returning from update_state, jbsupXmn, jpmn have been gathered onto all processors
       CALL GetFields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, isym)
       IF (lasym)                                                       &
       CALL GetFields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh, iasym)
       
!Convert jac*bsubXh to bsupXh and jac*pres to presh in real space

!On exit, bsupXf, pf0 valid on [nsmin, nsmax-1]
       CALL bhalftobfull (bsupsijh0, bsupuijh0, bsupvijh0,              &
                          bsupsijf0, bsupuijf0, bsupvijf0,              &
                          pijh0, pijf0)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       CALL cv_currents (bsupsijh0, bsupuijh0, bsupvijh0,               &
                         ksupsijf0, ksupuijf0, ksupvijf0,               &
                         l_getfsq, .TRUE.)

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
 100   CONTINUE
 
       IF (ladd_pert0 .OR. muPar.NE.zero) THEN
          nmax=MIN(endglobrow+1,ns)
          CALL tolowerf(bsupsijf0,bsupuijf0,bsupvijf0,                  & ! Calculate K || B
                        bsubsijf, bsubuijf, bsubvijf,                   &
                        nsmin, nmax)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
          IF (ladd_pert0) THEN
             ksubsijf(:,:,nsmin:nmax) = bsubsijf(:,:,nsmin:nmax)
             ksubuijf(:,:,nsmin:nmax) = bsubuijf(:,:,nsmin:nmax)
             ksubvijf(:,:,nsmin:nmax) = bsubvijf(:,:,nsmin:nmax)
             RETURN
          END IF
       END IF
     
!Compute spectrally filtered dp/du and dp/dv on half radial grid [nsmin:nsmax]
       IF (l_update_state) ALLOCATE (pfilter(ntheta, nzeta, nsmin:nsmax))
       CALL FilterPressure(isym)
       IF (lasym)                                                       &
       CALL FilterPressure(iasym)

!Compute radial pressure derivative pijf0_ds at full mesh points [nsmin:nsmax-1]
       CALL GradientFull(pijf0_ds, pijh0)

!SPH: one-sided (m=1) derivative at origin yields factor of 2 (1/(hs/2))
       IF (nsmin .EQ. 1) THEN
          pijf0_ds(:,:,1) = 2*(pijh0(:,:,2)-pijh0(:,:,1))*ohs
       END IF
      
!SPH11-3-16 preserve s=1 as iso-pressure contour
       IF (nsmax.EQ.ns .AND. l_pedge) THEN
          pijf0(:,:,ns) = pedge
          pijf0_ds(:,:,ns) = 0
       END IF

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

       CALL CheckFFT(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, isym)
       IF (lasym)                                                       &
       CALL CheckFFT(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh, iasym)

       IF (.NOT.lpar) THEN
          startglobrow=startr;  endglobrow=endr 
       END IF

       END SUBROUTINE init_state


!>  \brief subroutine for allocating unperturbed array
       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE shared_data, ONLY: l_par_state
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
       n1 = ntheta; n2 = nzeta
       
       IF (lpar) THEN
          n3 = nsmax1-nsmin1+1
          l_par_state=.TRUE.
       ELSE
          n3 = ns
          l_par_state=.FALSE.
       END IF

       IF (lalloc) THEN
          IF (SIZE(bsupsijh0,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (jvsupsijf, jvsupuijf, jvsupvijf,                  &
                      bsupsijf0, bsupuijf0, bsupvijf0,                  &
                      bsupsijf,  bsupuijf,  bsupvijf,                     &
                      bsupsijh0, bsupuijh0, bsupvijh0,                  &
                      bsubsijf,  bsubuijf,  bsubvijf,                   &
                      ksupsijf0, ksupuijf0, ksupvijf0,                  &
                      ksupsijf,  ksupuijf,  ksupvijf,                     &
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
          ALLOCATE (bsupsijh0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (half)
                    bsupuijh0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijh0(n1,n2,nsmin1:nsmax1),                     &
                    bsupsijf0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (full)
                    bsupuijf0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijf0(n1,n2,nsmin1:nsmax1),                     &
                    ksupsijf0(n1,n2,nsmin1:nsmax),                      & ! K^s, K^u, K^v 
                    ksupuijf0(n1,n2,nsmin1:nsmax),                      &
                    ksupvijf0(n1,n2,nsmin1:nsmax),                      &
                    ksubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    ksubuijf(n1,n2,nsmin1:nsmax),                       & ! K_s, K_u, K_v
                    ksubvijf(n1,n2,nsmin1:nsmax),                       &
                    bsubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    bsubuijf(n1,n2,nsmin1:nsmax),                       & ! B_s, B_u, B_v
                    bsubvijf(n1,n2,nsmin1:nsmax),                       &
                    ksupsijf(n1,n2,nsmin:nsmax),                           & 
                    ksupuijf(n1,n2,nsmin:nsmax),                           & 
                    ksupvijf(n1,n2,nsmin:nsmax),                           & 
                    bsupsijf(n1,n2,nsmin:nsmax),                           &
                    bsupuijf(n1,n2,nsmin:nsmax),                           &
                    bsupvijf(n1,n2,nsmin:nsmax),                           &
                    pijh0(n1,n2,nsmin1:nsmax1),                         &
                    pijh0_du(n1,n2,nsmin1:nsmax),                       &
                    pijh0_dv(n1,n2,nsmin1:nsmax),                       &
                    pijf0(n1,n2,nsmin1:nsmax1),                         &
                    pijf0_ds(n1,n2,nsmin1:nsmax), stat=istat)
          ELSE
          ALLOCATE (ksupsijf0(n1,n2,n3), ksupuijf0(n1,n2,n3),           &
                    ksupvijf0(n1,n2,n3), ksupsijf(n1,n2,n3),               &
                    ksupuijf(n1,n2,n3),  ksupvijf(n1,n2,n3),               &
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
                    bsupsijf(n1,n2,n3), bsupuijf(n1,n2,n3),               &
                    bsupvijf(n1,n2,n3), bsupsijh0(n1,n2,n3),              &
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
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) ::        &
             jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
       INTEGER, INTENT(IN) :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: fours, fouru, fourv, fcomb
!-----------------------------------------------
       IF (ipar .EQ. isym) THEN
          fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
       ELSE
          fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
       END IF

       CALL toijsp_par(jbsupsmnh(:,:,nsmin:nsmax), bsupsijh0, fcomb,    &
                       fours)
       CALL toijsp_par(jbsupumnh(:,:,nsmin:nsmax), bsupuijh0, fcomb,    &
                       fouru)
       CALL toijsp_par(jbsupvmnh(:,:,nsmin:nsmax), bsupvijh0, fcomb,    &
                       fourv)
       CALL toijsp_par(jpmnh(:,:,nsmin:nsmax), pijh0, fcomb,            &
                       fourv)

!
!check origin values
       IF (nsmin .EQ. 1) THEN
          CALL ASSERT(ALL(jbsupsmnh(m1,:,1).EQ.jbsupsmnh(m1,:,2)), 'bmnsh(1) != bmnsh(2)')
          CALL ASSERT(ALL(jbsupsmnh(m0,:,1).EQ.zero), 'bmnsh(m0:,1) != 0')
          CALL ASSERT(ALL(jbsupsmnh(m1+1:,:,1).EQ.zero), 'bmnsh(m2:,1) != 0')
          CALL ASSERT(ALL(jbsupumnh(m0:m1,:,1).EQ.jbsupumnh(m0:m1,:,2)), 'bmnuh(1) != bmnuh(2)')
          CALL ASSERT(ALL(jbsupumnh(m1+1:,:,1).EQ.zero), 'bmnuh(m2:,1) != 0')
          CALL ASSERT(ALL(jbsupvmnh(m0,:,1).EQ.jbsupvmnh(m0,:,2)), 'bmnvh(m0,1) != bmnvh(m0,2)')
          CALL ASSERT(ALL(jbsupvmnh(m1:,:,1).EQ.zero), 'bmnvh(m1:,1) != 0')
          CALL ASSERT(ALL(jpmnh(m0,:,1).EQ.jpmnh(m0,:,2)), 'pmnh(m0,1) != pmnh(m0,2)')
          CALL ASSERT(ALL(jpmnh(m1:,:,1).EQ.zero), 'pmnh(m1:,1) != 0')
       END IF          

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
       INTEGER  :: four, fcomb, mfour, nfour, istat
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: pmnh
!-----------------------------------------------
       IF (ipar .EQ. isym) THEN
          four = f_cos; fcomb = f_none
       ELSE
          four = f_sin; fcomb = f_sum
       END IF

       mfour = IOR(f_du, fcomb);  nfour = IOR(f_dv, fcomb)

!pmnh contains Nyquist-filtered (to MPOL, NTOR) Fourier harmonics of p on half mesh
       ALLOCATE (pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation error in FilterPressure')

       CALL tomnsp_par(pijh0, pmnh, four)
       
       IF (ipar.EQ.isym .AND. nsmax.EQ.ns) pedge = pmnh(m0, 0,ns)*orthonorm(0,0)
       
!spectrally filtered angle derivatives of p in real space (half mesh)
       CALL toijsp_par(pmnh, pijh0_du, mfour, four)
       CALL toijsp_par(pmnh, pijh0_dv, nfour, four)
       
       IF (l_update_state)                                              &
       CALL toijsp_par(pmnh, pfilter, fcomb, four)                       ! p (fourier filtered)

       DEALLOCATE (pmnh)

       END SUBROUTINE FilterPressure
       
       
       SUBROUTINE CheckFFT(jbsupsmnh, jbsupumnh, jbsupvmnh,             &
                           jpmnh, ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN) ::           &
             jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
       INTEGER, INTENT(IN) :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp), PARAMETER :: rtol = 1.E-14_dp
       INTEGER  :: fours, fouru, fourv, fcomb, sparity, fouru1, fourv1, &
                   m, mp, n, np
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: FFT_TEST, pmnh, temp
       REAL(dp) :: rtest
!-----------------------------------------------
#undef _TEST_INIT
!#define _TEST_INIT
#if defined(_TEST_INIT)
       IF (ipar .EQ. isym) THEN
          fours = f_sin; fouru = f_cos; fourv = f_cos; sparity=1
          fouru1 = f_sin; fourv1 = f_sin          
       ELSE
          fours = f_cos; fouru = f_sin; fourv = f_sin; sparity=-1
          fouru1 = f_cos; fourv1 = f_cos          
       END IF

!check independent variables
       ALLOCATE (fft_test(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                 temp(ntheta,nzeta,nsmin:nsmax),                        &
                 pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=m)
       CALL ASSERT(m.EQ.0,' Allocation failed in CheckFFT')
       temp = bsupsijh0*jacobh(:,:,nsmin:nsmax)
       CALL tomnsp_par(temp, fft_test, fours)
       rTest = MAXVAL(ABS(fft_test - jbsupsmnh(:,:,nsmin:nsmax)))
       CALL ASSERT(rtest.LT.rtol, ' jbsupsmnh FFT error in siesta_init')
       temp = bsupuijh0*jacobh(:,:,nsmin:nsmax)
       CALL tomnsp_par(temp, fft_test, fouru)
       rTest = MAXVAL(ABS(fft_test - jbsupumnh(:,:,nsmin:nsmax)))
       CALL ASSERT(rTest.LT.rtol,' jbsupumnh FFT error in siesta_init')
       temp = bsupvijh0*jacobh(:,:,nsmin:nsmax)
       CALL tomnsp_par(temp, fft_test, fourv)
       rTest = MAXVAL(ABS(fft_test - jbsupvmnh(:,:,nsmin:nsmax)))
       CALL ASSERT(rTest.LT.rtol,' jbsupvmnh FFT error in siesta_init')
       temp = pijh0*jacobh(:,:,nsmin:nsmax)
       CALL tomnsp_par(temp, fft_test, fourv)
       rTest = MAXVAL(ABS(fft_test - jpmnh(:,:,nsmin:nsmax)))
       CALL ASSERT(rTest.LT.rtol,' jpmnh FFT error in siesta_init')
            
! check angle derivatives of pressure            
       CALL tomnsp_par(pijh0, pmnh, fouru)
       CALL tomnsp_par(pijh0_du, fft_test, fouru1)
       DO m = 0, mpol
          mp = -sparity*m
          CALL ASSERT(MAXVAL(ABS(fft_test(m,:,:) - mp*pmnh(m,:,:))) .LE. rtol,   &
                     ' pijh0_du FFT error in siesta_init')
       END DO
       
       CALL tomnsp_par(pijh0_dv, fft_test, fourv1)
       DO n = -ntor, ntor
          np = -sparity*n*nfp
          CALL ASSERT(MAXVAL(ABS(fft_test(:,n,:) - np*pmnh(:,n,:))) .LE. rtol,   &
                     ' pijh0_dv FFT error in siesta_init')
       END DO
            
       DEALLOCATE(fft_test, pmnh, temp)
#endif                       
       
       END SUBROUTINE CheckFFT

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
