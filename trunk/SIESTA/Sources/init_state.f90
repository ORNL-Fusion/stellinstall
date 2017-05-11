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
       USE stel_kinds
       USE quantities
       USE diagnostics_mod
       USE evolution, ONLY: ste
       USE shared_data, ONLY: fsq_total, l_getfsq, l_update_state
       USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_sin, f_cos,  &
                          f_sum, f_full
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       USE perturbation, ONLY: buv_res
       USE hessian, ONLY: muPar
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1
       INTEGER            :: js, istat,  m, n, nblock
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: pmnch, pmnsh
       REAL(dp)    :: ton, toff, stes(2), rtes(2), pres_edge
       INTEGER     :: i, nsmin, nsmax, ns_span, n1, n2
       LOGICAL     :: ladd_pert, lcurr
!-----------------------------------------------
!      NOTE: jbupXmnPh, jpmnPh are initially computed on FULL ns mesh, and are GATHERED to that 
!            FULL ns mesh in update_state

       istat = 0
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       ns_span = nsmax-nsmin+1

       CALL Init_Allocate_Arrays (.TRUE.)

!Calculate unperturbed half mesh jac*bsupX
!SPH030410: need to get values at js=1, where UPDATE_BFIELD stored origin values
       CALL toijsp_par(jbsupsmnsh(:,:,nsmin:nsmax), bsupsijh0, f_none, f_sin, f_full, 1, ns_span)
       CALL toijsp_par(jbsupumnch(:,:,nsmin:nsmax), bsupuijh0, f_none, f_cos, f_full, 1, ns_span)
       CALL toijsp_par(jbsupvmnch(:,:,nsmin:nsmax), bsupvijh0, f_none, f_cos, f_full, 1, ns_span)

       IF (l_asym) THEN

       CALL toijsp_par(jbsupsmnch(:,:,nsmin:nsmax), bsupsijh0, f_sum, f_cos, f_full, 1, ns_span)
       CALL toijsp_par(jbsupumnsh(:,:,nsmin:nsmax), bsupuijh0, f_sum, f_sin, f_full, 1, ns_span)
       CALL toijsp_par(jbsupvmnsh(:,:,nsmin:nsmax), bsupvijh0, f_sum, f_sin, f_full, 1, ns_span)

       END IF

!Convert jac*bsubX to bsupX in real space
       bsupsijh0 = bsupsijh0/jacobh(:,:,nsmin:nsmax)
       bsupvijh0 = bsupvijh0/jacobh(:,:,nsmin:nsmax)


!On exit, bsupXijPf0 valid on [nsmin, nsmax-1] (P=parity, c, s)
       CALL bhalftobfull_par (bsupsijh0, bsupuijh0, bsupvijh0,           &
                              bsupsijf0, bsupuijf0, bsupvijf0,           &
                              nsmin, nsmax)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       END IF

       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents_par (bsupsijh0(:,:,nsmin:nsmax),                 &
                             bsupuijh0(:,:,nsmin:nsmax),                 &
                             bsupvijh0(:,:,nsmin:nsmax),                 &
                             ksupsijf0(:,:,nsmin:nsmax),                 &
                             ksupuijf0(:,:,nsmin:nsmax),                 &
                             ksupvijf0(:,:,nsmin:nsmax),                 &
                             l_getfsq, lcurr )
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       END IF

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          nsmax=MIN(endglobrow+1,ns)
          CALL tolowerf_par(bsupsijf0,bsupuijf0,bsupvijf0,               & ! Calculate K || B
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

!Calculate pressure in real space (half mesh) - need js=1 point here
       CALL toijsp_par(jpmnch(:,:,nsmin:nsmax), pijh, f_none, f_cos, f_full, 1, ns_span)
       IF (l_asym) THEN
       CALL toijsp_par(jpmnsh(:,:,nsmin:nsmax), pijh, f_sum, f_sin, f_full, 1, ns_span)
       END IF

!      Get real pressure from pijh = jacobh * real-pressure
       pijh = pijh/jacobh(:,:,nsmin:nsmax)

!      pressure on full mesh (extrap to origin stored in pijh(1)
!      On Exit, pijf valid on [nsmin:nsmax-1]
       nblock=SIZE(pijh,1)*SIZE(pijh,2)
       CALL to_full_mesh_par(pijh, pijf, nblock, nsmin, nsmax)            ! Get pressure on full mesh

!      Get real pressure p on half grid needed for dp/du and dp/dv
       ALLOCATE (pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

!      pmnc(s)h contains filtered (to MPOL, NTOR) Fourier harmonics of p on half mesh
       CALL tomnsp_par(pijh, pmnch(:,:,nsmin:nsmax), f_cos, f_full, 1, ns_span)
       IF (l_asym) THEN
       ALLOCATE (pmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'
       CALL tomnsp_par(pijh, pmnsh(:,:,nsmin:nsmax), f_sin, f_full, 1, ns_span)
       END IF

       IF (l_update_state) THEN
          pijh_du = pijh                                                  ! unfiltered; use pijh_dv for work array
          CALL toijsp_par(pmnch, pijh_dv, f_none, f_cos, f_half, 1, ns_span) ! p (fourier filtered)
          IF (l_asym) THEN
          CALL toijsp_par(pmnsh, pijh_dv, f_sum, f_sin, f_half, 1, ns_span)  ! p (fourier filtered)
          END IF
          n1=MAX(1,startglobrow); n2=MIN(endglobrow,ns)
          IF (n1 .EQ. 1) THEN
             pijh_du(:,:,1) = 0;   pijh_dv(:,:,1) = 0
          END IF
          stes(1) = SUM((pijh_du(:,:,n1:n2)-pijh_dv(:,:,n1:n2))**2)
          stes(2) = SUM((pijh_du(:,:,n1:n2)+pijh_dv(:,:,n1:n2))**2)
          CALL MPI_ALLREDUCE(stes,rtes,2,MPI_REAL8,MPI_SUM,              &
                             MPI_COMM_WORLD,MPI_ERR)
          IF (rtes(2) .NE. zero) ste(1) = SQRT(rtes(1)/rtes(2))
       END IF

!      compute radial derivative at full mesh points
!      On Exit, pijf_ds valid on [nsmin:nsmax-1]
       CALL GradientFull_par(pijf_ds, pijh, nblock, nsmin, nsmax)

!SPH:  one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
!!       CALL getorigin(pijh, m1, 0)
       IF (nsmin .EQ. 1) THEN
          pijf_ds(:,:,1) = 2*pijh(:,:,1)*ohs
          pijh(:,:,1) = 0 
       END IF

!      spectrally filtered angle derivatives
       CALL toijsp_par(pmnch, pijh_du, f_du, f_cos, f_half, 1, ns_span)   ! dp/dtheta in real space (half mesh)
       CALL toijsp_par(pmnch, pijh_dv, f_dv, f_cos, f_full, 1, ns_span)   ! dp/dzeta in real space (half mesh, include js=1 extrap)
       DEALLOCATE(pmnch)                 
       IF (l_asym) THEN
       m = IOR(f_du,f_sum);  n = IOR(f_dv, f_sum)
       CALL toijsp_par(pmnsh, pijh_du, m, f_sin, f_half, 1, ns_span)       ! dp/dtheta in real space (half mesh)
       CALL toijsp_par(pmnsh, pijh_dv, n, f_sin, f_full, 1, ns_span)       ! dp/dzeta in real space (half mesh, include js=1 extrap)
       DEALLOCATE(pmnsh)                 
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
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE quantities, ONLY: ns
       IMPLICIT NONE
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

       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE quantities
       USE shared_data, ONLY: l_par_state
       USE nscalingtools, ONLY: startglobrow, endglobrow
       IMPLICIT NONE
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
                    jvsupvijf(n1,n2,n3))                                   ! V_s (cosine), V_u (sine), V_v (sine)
          IF (lpar) THEN
          ALLOCATE (ksupsijf0(n1,n2,nsmin1:nsmax1),                      & ! K^s (sine), K^u (cosine), K^v (cosine)
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
          IF (SIZE(ksupsmnsf,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (ksupsmnsf, ksupumncf, ksupvmncf,                   &
                      djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,    &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          IF (l_asym) THEN
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
             IF (l_asym) THEN
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
             IF (l_asym) THEN
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
             IF (l_asym) THEN
             djpmnsh = 0; djbsupsmnch = 0; djbsupumnsh = 0; djbsupvmnsh = 0
             END IF
          ELSE
             STOP 'Allocation error #2 in Allocate_Arrays'
          END IF
       END IF

       END SUBROUTINE Init_Allocate_Arrays
