       MODULE siesta_init
       USE stel_kinds
       USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_sin, f_cos,  &
                          f_du, f_dv, f_sum
       USE shared_data, ONLY: l_update_state
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       IMPLICIT NONE
       INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
       INTEGER, PRIVATE            :: nsmin, nsmax, ns_span
       REAL(dp), ALLOCATABLE, PRIVATE :: pfilter(:,:,:)

       CONTAINS

       SUBROUTINE init_state_par(lcurrent_only)
!     
!     WRITTEN 01-12-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
!#if defined(SKS)
       USE shared_data, ONLY: fsq_total, l_getfsq, ste, buv_res
       USE hessian, ONLY: muPar
       USE siesta_currents, ONLY: cv_currents
       USE quantities
       USE mpi_inc
       USE descriptor_mod, ONLY: SIESTA_COMM
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1
       INTEGER            :: js, istat,  m, n, nblock
       REAL(dp)           :: ton, toff, stes(2), rtes(2), pres_edge
       INTEGER            :: i, n1, n2
       LOGICAL            :: ladd_pert, lcurr
!-----------------------------------------------
!      NOTE: jbupXmnPh, jpmnPh are initially computed on FULL ns mesh, and are GATHERED to that 
!            FULL ns mesh in update_state

       istat = 0
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       ns_span = nsmax-nsmin+1

       CALL Init_Allocate_Arrays (.TRUE.)

!Calculate unperturbed half mesh jac*bsupX and jac*pressure
!Sum symmetric and (for lasym) asymmetric parts
       CALL GetFields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, isym)
       IF (lasym)                                                        &
       CALL GetFields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh, iasym)

! Convert jac*bsupX to bsupX in real space. Jacobian for bsupu term is removed
! in bhalftobfull.
       bsupsijh0 = bsupsijh0/jacobh(:,:,nsmin:nsmax)
       bsupvijh0 = bsupvijh0/jacobh(:,:,nsmin:nsmax)

!On exit, bsupXijf0 valid on [nsmin, nsmax-1]
       CALL bhalftobfull_par (bsupsijh0, bsupuijh0, bsupvijh0,           &
                              bsupsijf0, bsupuijf0, bsupvijf0,           &
                              nsmin, nsmax)

!Convert jac*pijh to pijh in real space
       pijh = pijh/jacobh(:,:,nsmin:nsmax)

!On exit, pijf valid on [nsmin:nsmax-1]
       nblock=SIZE(pijh,1)*SIZE(pijh,2)
       CALL to_full_mesh_par(pijh, pijf, nblock, nsmin, nsmax)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents (bsupsijh0, bsupuijh0, bsupvijh0,                &
                         ksupsijf0, ksupuijf0, ksupvijf0,                &
                         l_getfsq, .TRUE.)

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          nsmax=MIN(endglobrow+1,ns)
! Calculate K || B. Bsub contains the vector direction. A scalar scale factor is 
! included in the helical perturbation to scale B to K. K|| = sigma*B
          CALL tolowerf_par(bsupsijf0,bsupuijf0,bsupvijf0,               &
                            bsubsijf, bsubuijf, bsubvijf,                &
                            nsmin, nsmax)  
!IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
          IF (ladd_pert) THEN
          ksubsijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)
          ksubuijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)
          ksubvijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)
          END IF
          nsmax=MIN(endglobrow+2,ns)
       END IF
 
!Compute real pressure p Fourier coefficients on half grid
!and compute dp/du and dp/dv
       IF (l_update_state) ALLOCATE (pfilter(ntheta, nzeta, nsmin:nsmax))
       CALL FilterPressure(isym)
       IF (lasym)                                                        &
       CALL FilterPressure(iasym)

!      compute radial derivative at full mesh points
!      On Exit, pijf_ds valid on [nsmin:nsmax-1]
       CALL GradientFull_par(pijf_ds, pijh, nblock, nsmin, nsmax)

!SPH:  one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
       IF (nsmin .EQ. 1) THEN
          pijf_ds(:,:,1) = 2*pijh(:,:,1)*ohs
          pijh(:,:,1) = 0
       END IF
      
!SPH   AVERAGE pijh_du,dv to FULL mesh and set origin value (SIN SYMMETRY)
!      On Exit, pijf_du,v valid on [nsmin:nsmax-1], with extension to nsmax=ns
       CALL to_full_mesh_par(pijh_du, pijf_du, nblock, nsmin, nsmax)
       CALL to_full_mesh_par(pijh_dv, pijf_dv, nblock, nsmin, nsmax)
!SPH11-3-16 s=1 is iso-pressure contour
       IF (nsmax .EQ. ns) THEN
          CALL SurfAverage(pres_edge, pijh(:,:,ns), ns, ns)
          pijf(:,:,ns) = pres_edge
          pijf_du(:,:,ns) = 0
          pijf_dv(:,:,ns) = 0
       END IF

!Fourier filter (Nyquist) check
       IF (.NOT.l_update_state) RETURN
!pijh is unfiltered
       n1=MAX(1,startglobrow); n2=MIN(endglobrow,ns)
       IF (n1 .EQ. 1) pfilter(:,:,1) = 0
       stes(1) = SUM((pijh(:,:,n1:n2)-pfilter(:,:,n1:n2))**2)
       stes(2) = SUM((pijh(:,:,n1:n2)+pfilter(:,:,n1:n2))**2)
#if defined(MPI_OPT)
       CALL MPI_ALLREDUCE(stes,rtes,2,MPI_REAL8,MPI_SUM,                 &
                          SIESTA_COMM,MPI_ERR)
       IF (rtes(2) .NE. zero) ste(1) = SQRT(rtes(1)/rtes(2))
#else
       IF (stes(2) .NE. zero) ste(1) = SQRT(stes(1)/stes(2))
#endif
       DEALLOCATE (pfilter)
!#endif
       END SUBROUTINE init_state_par

       SUBROUTINE init_state(lcurrent_only)
!     
!     WRITTEN 01-12-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
       USE quantities, ONLY: ns
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: srow, erow
!-----------------------------------------------
       srow = startglobrow; erow = endglobrow
       startglobrow = 1;  endglobrow = ns

       CALL init_state_par(lcurrent_only)

       startglobrow = srow; endglobrow = erow

       END SUBROUTINE init_state


       SUBROUTINE GetFields(jbsupsmnh, jbsupumnh, jbsupvmnh,             &
                            jpmnh, iparity)
       USE quantities, ONLY: bsupsijh0, bsupuijh0, bsupvijh0, pijh
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                        &
!       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns) ::                      &
             jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
       INTEGER, INTENT(IN) :: iparity
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: fours, fouru, fourv, icomb
!-----------------------------------------------
       IF (iparity .EQ. isym) THEN
          icomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
       ELSE
          icomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
       END IF

       CALL toijsp_par(jbsupsmnh(:,:,nsmin:nsmax), bsupsijh0, icomb,     &
                       fours, 1, ns_span)
       CALL toijsp_par(jbsupumnh(:,:,nsmin:nsmax), bsupuijh0, icomb,     &
                       fouru, 1, ns_span)
       CALL toijsp_par(jbsupvmnh(:,:,nsmin:nsmax), bsupvijh0, icomb,     &
                       fourv, 1, ns_span)
       CALL toijsp_par(jpmnh(:,:,nsmin:nsmax), pijh, icomb,              &
                       fourv, 1, ns_span)

       END SUBROUTINE GetFields

       SUBROUTINE FilterPressure(iparity)
       USE quantities, ONLY: mpol, ntor, pijh, pijh_du, pijh_dv
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: iparity
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: four, icomb, mfour, nfour, istat
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: pmnh
!-----------------------------------------------
       IF (iparity .EQ. isym) THEN
          four = f_cos; icomb = f_none
          mfour = f_du; nfour = f_dv
       ELSE
          four = f_sin; icomb=f_sum
          mfour = IOR(f_du,f_sum);  nfour = IOR(f_dv, f_sum)
       END IF

!pmnh contains filtered (to MPOL, NTOR) Fourier harmonics of p on half mesh
       ALLOCATE (pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

       CALL tomnsp_par(pijh, pmnh(:,:,nsmin:nsmax), four,                &
                       1, ns_span)

!      spectrally filtered angle derivatives
       CALL toijsp_par(pmnh, pijh_du, mfour, four, 1, ns_span)            ! dp/du in real space (half mesh)
       CALL toijsp_par(pmnh, pijh_dv, nfour, four, 1, ns_span)            ! dp/dv in real space (half mesh, include js=1 extrap)

       IF (l_update_state)                                               &
       CALL toijsp_par(pmnh, pfilter, icomb, four, 1, ns_span)            ! p (fourier filtered)

       DEALLOCATE (pmnh)

       END SUBROUTINE FilterPressure

       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE quantities
       USE shared_data, ONLY: l_par_state
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)     :: lpar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL                 :: lrealloc, lalloc
       INTEGER                 :: istat, n1, n2, n3, nsmin, nsmax,       &
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
          DEALLOCATE (jvsupsijf, jvsupuijf, jvsupvijf,                   &
                      bsupsijf0, bsupuijf0, bsupvijf0,                   &
                      bsupsijh0, bsupuijh0, bsupvijh0,                   &
                      bsubsijf,  bsubuijf,  bsubvijf,                    &
                      ksupsijf0, ksupuijf0, ksupvijf0,                   &
                      ksubsijf,  ksubuijf,  ksubvijf,                    &
                      pijf, pijh, pijf_ds, pijf_du, pijf_dv,             &
                      pijh_du, pijh_dv, stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          END IF
       END IF

       IF (lrealloc) THEN
          n3 = ns
          ALLOCATE (jvsupsijf(n1,n2,n3), jvsupuijf(n1,n2,n3),            & ! Full mesh quantities (real space)
                    jvsupvijf(n1,n2,n3))                                   ! V_s (cos), V_u (sin), V_v (sin)
          IF (lpar) THEN
          ALLOCATE (ksupsijf0(n1,n2,nsmin1:nsmax1),                      & ! K^s (sin), K^u (cos), K^v (cos)
                    ksupuijf0(n1,n2,nsmin1:nsmax1),                      &
                    ksupvijf0(n1,n2,nsmin1:nsmax1),                      &
                    bsupsijf0(n1,n2,nsmin1:nsmax1),                      & ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijf0(n1,n2,nsmin1:nsmax1),                      &
                    bsupvijf0(n1,n2,nsmin1:nsmax1),                      &
                    bsupsijh0(n1,n2,nsmin1:nsmax1),                      & ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijh0(n1,n2,nsmin1:nsmax1),                      &
                    bsupvijh0(n1,n2,nsmin1:nsmax1),                      &
                    ksubsijf(n1,n2,nsmin1:nsmax),                        & ! Full mesh quantities (real space)
                    ksubuijf(n1,n2,nsmin1:nsmax),                        & ! K^s (sine), K^u (cosine), K^v (cosine)
                    ksubvijf(n1,n2,nsmin1:nsmax),                        &
                    bsubsijf(n1,n2,nsmin1:nsmax),                        & ! Full mesh quantities (real space)
                    bsubuijf(n1,n2,nsmin1:nsmax),                        & ! K^s (sine), K^u (cosine), K^v (cosine)
                    bsubvijf(n1,n2,nsmin1:nsmax),                        &
                    pijh(n1,n2,nsmin1:nsmax1),                           &
                    pijf(n1,n2,nsmin1:nsmax1),                           &
                    pijf_ds(n1,n2,nsmin1:nsmax1),                        &
                    pijf_du(n1,n2,nsmin1:nsmax1),                        &
                    pijf_dv(n1,n2,nsmin1:nsmax1),                        &
                    pijh_du(n1,n2,nsmin1:nsmax1),                        &
                    pijh_dv(n1,n2,nsmin1:nsmax1), stat=istat)
          ELSE
          ALLOCATE (ksupsijf0(n1,n2,n3), ksupuijf0(n1,n2,n3),            &
                    ksupvijf0(n1,n2,n3),                                 &
                    ksubsijf(n1,n2,n3),                                  &
                    ksubuijf(n1,n2,n3),                                  &
                    ksubvijf(n1,n2,n3),                                  &
                    bsubsijf(n1,n2,n3),                                  &
                    bsubuijf(n1,n2,n3),                                  &
                    bsubvijf(n1,n2,n3),                                  &
                    pijf(n1,n2,n3), pijh(n1,n2,n3),                      &
                    pijf_ds(n1,n2,n3), pijf_du(n1,n2,n3),                &
                    pijf_dv(n1,n2,n3), pijh_du(n1,n2,n3),                &
                    pijh_dv(n1,n2,n3), bsupsijf0(n1,n2,n3),              &
                    bsupuijf0(n1,n2,n3),bsupvijf0(n1,n2,n3),             &
                    bsupsijh0(n1,n2,n3),                                 &
                    bsupuijh0(n1,n2,n3),bsupvijh0(n1,n2,n3),             &
                    stat=istat)
          END IF
          IF (istat .EQ. 0) THEN
             jvsupsijf = 0; jvsupuijf = 0; jvsupvijf = 0                   ! Need in add_resistivity loop
          ELSE
             STOP 'Allocation error #1 in Allocate_Arrays'
          END IF
       END IF


       lalloc = ALLOCATED(ksupsmnsf) 
       lrealloc = .TRUE.

       IF (lpar) THEN
          n1 = mpol+1; n2 = 2*ntor+1; n3 = nsmax-nsmin1+1
       ELSE
          n1 = mpol+1; n2 = 2*ntor+1; n3 = ns
       END IF

       IF (lalloc) THEN
          IF (SIZE(ksupsmnsf,3).EQ.n3 .AND. SIZE(djpmnch,3).EQ.n3)       &
             lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (ksupsmnsf, ksupumncf, ksupvmncf,                   &
                      djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,    &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          IF (lasym) THEN
          DEALLOCATE (ksupsmncf, ksupumnsf, ksupvmnsf,                   &
                      djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh,    &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          END IF
          END IF
       END IF

       IF (lrealloc) THEN
          IF (lpar) THEN
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),        & 
                       ksupumncf(0:mpol,-ntor:ntor,nsmin1:nsmax),        &
                       ksupvmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),        &
                       djpmnch(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                       djbsupsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       djbsupumnch(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       djbsupvmnch(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       stat=istat)
             IF (lasym) THEN
              ALLOCATE(ksupsmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),        & 
                       ksupumnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),        &
                       ksupvmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),        &
                       djpmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                       djbsupsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       djbsupumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       djbsupvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),       &
                       stat=istat)
             END IF

          ELSE
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,n3),                  & 
                       ksupumncf(0:mpol,-ntor:ntor,n3),                  &
                       ksupvmncf(0:mpol,-ntor:ntor,n3),                  & 
                       djpmnch(0:mpol,-ntor:ntor,n3),                    &
                       djbsupsmnsh(0:mpol,-ntor:ntor,n3),                &
                       djbsupumnch(0:mpol,-ntor:ntor,n3),                &
                       djbsupvmnch(0:mpol,-ntor:ntor,n3),                &
                       stat=istat)
             IF (lasym) THEN
              ALLOCATE(ksupsmncf(0:mpol,-ntor:ntor,n3),                  & 
                       ksupumnsf(0:mpol,-ntor:ntor,n3),                  &
                       ksupvmnsf(0:mpol,-ntor:ntor,n3),                  & 
                       djpmnsh(0:mpol,-ntor:ntor,n3),                    &
                       djbsupsmnch(0:mpol,-ntor:ntor,n3),                &
                       djbsupumnsh(0:mpol,-ntor:ntor,n3),                &
                       djbsupvmnsh(0:mpol,-ntor:ntor,n3),                &
                       stat=istat)
             END IF
          END IF
          IF (istat .EQ. 0) THEN
             djpmnch = 0; djbsupsmnsh = 0; djbsupumnch = 0; djbsupvmnch = 0 
             IF (lasym) THEN
             djpmnsh = 0; djbsupsmnch = 0; djbsupumnsh = 0; djbsupvmnsh = 0
             END IF
          ELSE
             STOP 'Allocation error #2 in Allocate_Arrays'
          END IF
       END IF

       END SUBROUTINE Init_Allocate_Arrays

       END MODULE siesta_init
