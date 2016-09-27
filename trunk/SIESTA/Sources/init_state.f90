      SUBROUTINE init_state_par(lcurrent_only)
!     
!     WRITTEN 01-12-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
#if defined(SKS)
       USE stel_kinds
       USE quantities
       USE diagnostics_mod
       USE evolution, ONLY: ste
       USE shared_data, ONLY: fsq_total, l_getfsq, l_update_state
       USE fourier, ONLY: toijsp_par, tomnsp_par,                              &
                          fourier_parity_even, fourier_parity_odd,             &
                          fourier_full_grid, fourier_half_grid,                &
                          fourier_du, fourier_dv, fourier_none
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
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1
       LOGICAL, PARAMETER :: lpar=.TRUE.
       INTEGER            :: istat, m, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: pmnch, pmnsh, tempij,     &
                    pijh_du, pijh_dv 
       REAL(rprec) :: ton, toff, stes(2), rtes(2)
       INTEGER     :: i, nsmin, nsmax, ns_span, n1, n2, js
       LOGICAL     :: ladd_pert, lcurr
!-----------------------------------------------
       istat = 0
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       ns_span = nsmax-nsmin+1

       CALL Init_Allocate_Arrays (lpar)

! Calculate unperturbed half mesh jac*bsupX
! SPH030410: need to get values at js=1, where UPDATE_BFIELD stored origin values
       CALL toijsp_par(jbsupsmnsh(:,:,nsmin:nsmax), bsupsijh0, 0, 0, 1, ifull, 1, ns_span)
       CALL toijsp_par(jbsupumnch(:,:,nsmin:nsmax), bsupuijh0, 0, 0, 0, ifull, 1, ns_span)
       CALL toijsp_par(jbsupvmnch(:,:,nsmin:nsmax), bsupvijh0, 0, 0, 0, ifull, 1, ns_span)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(tempij(LBOUND(bsupsijh0, 1):UBOUND(bsupsijh0, 1),           &
                          LBOUND(bsupsijh0, 2):UBOUND(bsupsijh0, 2),           &
                          LBOUND(bsupsijh0, 3):UBOUND(bsupsijh0, 3)), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

          CALL toijsp_par(jbsupsmnch(:,:,nsmin:nsmax), tempij,                 &
                          fourier_none, fourier_none, fourier_parity_even, fourier_full_grid, 1, ns_span)
          bsupsijh0 = bsupsijh0 + tempij

          CALL toijsp_par(jbsupumnsh(:,:,nsmin:nsmax), tempij,                 &
                          fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid, 1, ns_span)
          bsupuijh0 = bsupuijh0 + tempij

          CALL toijsp_par(jbsupvmnsh(:,:,nsmin:nsmax), tempij,                 &
                          fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid, 1, ns_span)
          bsupvijh0 = bsupvijh0 + tempij

          DEALLOCATE(tempij)
       END IF
#endif

! Convert jac*bsubX to bsupX in real space


       bsupsijh0 = bsupsijh0/jacobh(:,:,nsmin:nsmax)
       bsupvijh0 = bsupvijh0/jacobh(:,:,nsmin:nsmax)

!      On exit, bsupXF valid on [nsmin, nsmax-1]
       CALL bhalftobfull_par (bsupsijh0, bsupuijh0, bsupvijh0,       &
                              bsupsijf0, bsupuijf0, bsupvijf0,       &
                              nsmin, nsmax)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       END IF

       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents_par (bsupsijh0(:,:,nsmin:nsmax),               &
                             bsupuijh0(:,:,nsmin:nsmax),               &
                             bsupvijh0(:,:,nsmin:nsmax),               &
                             ksupsijf0(:,:,nsmin:nsmax),               &
                             ksupuijf0(:,:,nsmin:nsmax),               &
                             ksupvijf0(:,:,nsmin:nsmax),               &
                             l_getfsq, lcurr )
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       END IF

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          nsmax=MIN(endglobrow+1,ns)
          CALL tolowerf_par(bsupsijf0,bsupuijf0,bsupvijf0,           & ! Calculate K || B
                            bsubsijf, bsubuijf, bsubvijf,            &
                            nsmin, nsmax)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
          IF (ladd_pert) THEN
          ksubsijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)
          ksubuijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)
          ksubvijf(:,:,nsmin:nsmax) = jacobf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)
          END IF
          nsmax=MIN(endglobrow+2,ns)
       END IF

! Calculate pressure in real space (half mesh) - need js=1 point here
       CALL toijsp_par(jpmnch(:,:,nsmin:nsmax), pijh, 0, 0, 0, ifull, 1, ns_span)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(tempij(LBOUND(pijh, 1):UBOUND(pijh, 1),           &
                          LBOUND(pijh, 2):UBOUND(pijh, 2),           &
                          LBOUND(pijh, 3):UBOUND(pijh, 3)), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

          CALL toijsp_par(jpmnsh(:,:,nsmin:nsmax), pijh,                       &
                          fourier_none, fourier_none, fourier_parity_odd, ifull, fourier_full_grid, ns_span)

          DEALLOCATE(tempij)
       END IF
#endif

!      Get real pressure from pijh = jacobh * real-pressure
       pijh = pijh/jacobh(:,:,nsmin:nsmax)
 
!      pressure on full mesh (extrap to origin stored in pijh(1)
!      On Exit, pijf valid on [nsmin:nsmax-1]
       nblock=SIZE(pijh,1)*SIZE(pijh,2)
       CALL to_full_mesh_par(pijh, pijf, nblock, nsmin, nsmax)           ! Get pressure on full mesh

!      Get true pressure p on half grid needed for dp/du and dp/dv
       ALLOCATE (pmnch(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                 pijh_du(ntheta,nzeta,nsmin:nsmax),                     &
                 pijh_dv(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

!      pmnch now contains filtered (to MPOL, NTOR) fourier harmonics of p on half mesh
       CALL tomnsp_par(pijh, pmnch(:,:,nsmin:nsmax), 0, ifull, 1, ns_span)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

          CALL tomnsp_par(pijh, pmnsh(:,:,nsmin:nsmax), fourier_parity_odd, fourier_full_grid, 1, ns_span) ! Not sure why this in on full grid.
       END IF
#endif

       IF (l_update_state) THEN
          pijh_du = pijh                                                ! unfiltered
          CALL toijsp_par(pmnch, pijh_dv, 0, 0, 0, ihalf, 1, ns_span)    ! p (fourier filtered)

#if defined(MRC)
          IF (l_force_asym) THEN
             ALLOCATE(tempij(LBOUND(pijh_dv, 1):UBOUND(pijh_dv, 1),           &
                             LBOUND(pijh_dv, 2):UBOUND(pijh_dv, 2),           &
                             LBOUND(pijh_dv, 3):UBOUND(pijh_dv, 3)), stat=istat)
             IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

             CALL toijsp_par(pmnsh, tempij, fourier_none, fourier_none,        &
                             fourier_parity_odd, fourier_half_grid, 1, ns_span)
             pijh_dv = pijh_dv + tempij
          END IF
#endif

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
       CALL toijsp_par(pmnch, pijh_du, 1, 0, 0, ihalf, 1, ns_span)      ! dp/dtheta in real space (half mesh)
       CALL toijsp_par(pmnch, pijh_dv, 0, 1, 0, ifull, 1, ns_span)      ! dp/dzeta in real space (half mesh, include js=1 extrap)

       DEALLOCATE(pmnch)

#if defined(MRC)
       IF (l_force_asym) THEN
          WRITE (*,*) UBOUND(pijh_du, 3), LBOUND(pijh_du, 3), UBOUND(tempij, 3), LBOUND(tempij, 3)

          CALL toijsp_par(pmnsh, tempij, fourier_du, fourier_none,             &
                          fourier_parity_odd, fourier_half_grid, 1, ns_span)
          pijh_du = pijh_du + tempij

          CALL toijsp_par(pmnsh, tempij, fourier_none, fourier_dv,             &
                          fourier_parity_odd, fourier_full_grid, 1, ns_span)
          pijh_dv = pijh_dv + tempij

          DEALLOCATE(tempij)
          DEALLOCATE(pmnsh)
       END IF
#endif

!SPH   AVERAGE pijh_du,dv to FULL mesh and set origin value (SIN SYMMETRY)
!      On Exit, pijf_du,v valid on [nsmin:nsmax-1]
       CALL to_full_mesh_par(pijh_du, pijf_du, nblock, nsmin, nsmax)
       CALL to_full_mesh_par(pijh_dv, pijf_dv, nblock, nsmin, nsmax)

       DEALLOCATE(pijh_du, pijh_dv)                 
#endif
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
       USE stel_kinds
       USE quantities
       USE diagnostics_mod
       USE evolution, ONLY: ste
       USE shared_data, ONLY: fsq_total, l_getfsq, l_update_state
       USE perturbation, ONLY: buv_res
       USE hessian, ONLY: muPar
       USE fourier, ONLY: toijsp_par, tomnsp_par,                              &
                          fourier_parity_even, fourier_parity_odd,             &
                          fourier_full_grid, fourier_half_grid,                &
                          fourier_du, fourier_dv, fourier_none
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1
       LOGICAL, PARAMETER :: lpar=.FALSE.
       INTEGER            :: istat, m, nblock                                   !, mpol_s, ntor_s
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: pmnch, pmnsh, tempij,     &
                    pijh_du, pijh_dv 
       LOGICAL            :: ladd_pert, lcurr
!-----------------------------------------------
       istat = 0

       CALL Init_Allocate_Arrays (lpar)

! Calculate unperturbed half mesh jac*bsupX
! SPH030410: need to get values at js=1, where UPDATE_BFIELD stored origin values

       CALL toijsp(jbsupsmnsh, bsupsijh0, 0, 0, 1, ifull)
       CALL toijsp(jbsupumnch, bsupuijh0, 0, 0, 0, ifull)
       CALL toijsp(jbsupvmnch, bsupvijh0, 0, 0, 0, ifull)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(tempij(LBOUND(bsupsijh0, 1):UBOUND(bsupsijh0, 1),           &
                          LBOUND(bsupsijh0, 1):UBOUND(bsupsijh0, 1),           &
                          LBOUND(bsupsijh0, 1):UBOUND(bsupsijh0, 1)), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

          CALL toijsp(jbsupsmnch, tempij, fourier_none, fourier_none,          &
                      fourier_parity_even, fourier_full_grid)
          bsupsijh0 = bsupsijh0 + tempij

          CALL toijsp(jbsupumnsh, tempij, fourier_none, fourier_none,          &
                      fourier_parity_odd, fourier_full_grid)
          bsupuijh0 = bsupuijh0 + tempij

          CALL toijsp(jbsupvmnsh, tempij, fourier_none, fourier_none,          &
                      fourier_parity_odd, fourier_full_grid)
          bsupvijh0 = bsupvijh0 + tempij
       END IF
#endif

! Convert jac*bsubX to bsupX in real space

       bsupsijh0 = bsupsijh0/jacobh
       bsupvijh0 = bsupvijh0/jacobh

       CALL bhalftobfull (bsupsijh0, bsupuijh0, bsupvijh0,           &
                          bsupsijf0, bsupuijf0, bsupvijf0)
       IF (ANY(bsupsijf0(:,:,ns) .NE. zero)) PRINT *,'bsupsijf0(ns) != 0'


!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents (bsupsijh0, bsupuijh0, bsupvijh0,            &
                         ksupsijf0, ksupuijf0, ksupvijf0,            &
                         l_getfsq, lcurr)

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          CALL tolowerf(bsupsijf0, bsupuijf0, bsupvijf0,             & ! Calculate K || B
                        bsubsijf,  bsubuijf,  bsubvijf)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB
          IF (ladd_pert) THEN
             ksubsijf = jacobf*bsubsijf
             ksubuijf = jacobf*bsubuijf
             ksubvijf = jacobf*bsubvijf
          END IF
       END IF

! Calculate pressure in real space (half mesh) - need js=1 point here
       CALL toijsp(jpmnch, pijh, 0, 0, 0, ifull)

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL toijsp(jpmnsh, tempij, fourier_none, fourier_none,              &
                      fourier_parity_odd, fourier_full_grid)
          pijh = pijh + tempij
       END IF
#endif
              
!      Get real pressure from pijh = jacobh * real-pressure
       pijh = pijh/jacobh

!      pressure on full mesh (extrap to origin stored in pijh(1)
       nblock=SIZE(pijh,1)*SIZE(pijh,2)
       CALL to_full_mesh(pijh, pijf, nblock)                               ! Get pressure on full mesh
       pijf(:,:,1) = pijh(:,:,1)   
       pijf(:,:,ns) = pijh(:,:,ns)                                    ! 1pt "extrap"

!      Get true pressure p on half grid needed for dp/du and dp/dv
       ALLOCATE (pmnch(0:mpol,-ntor:ntor,ns),                           &
                 pijh_du(ntheta,nzeta,ns),                              &
                 pijh_dv(ntheta,nzeta,ns), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in INIT_STATE'

!      work1 now contains filtered (to MPOL, NTOR) fourier harmonics of p on half mesh
       CALL tomnsp(pijh, pmnch, 0, ifull)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsh(0:mpol,-ntor:ntor,ns), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation error in INIT_STATE'

          CALL tomnsp(pijh, pmnsh, fourier_parity_odd, fourier_full_grid)
       END IF
#endif

!      spectral filtering of p (compute spectral truncation error diagnostic)
       IF (l_update_state) THEN
          pijh_du = pijh                                                 ! unfiltered
          CALL toijsp(pmnch, pijh_dv, 0, 0, 0, ihalf)                     ! p (fourier filtered)
#if defined(MRC)
          IF (l_force_asym) THEN
             CALL toijsp(pmnsh, tempij, fourier_none, fourier_none,            &
                         fourier_parity_odd, fourier_half_grid)
             pijh_dv = pijh_dv + tempij
          END IF
#endif
          pijh_du(:,:,1) = 0;   pijh_dv(:,:,1) = 0
          ste(1) = SUM((pijh_du+pijh_dv)**2)
          IF (ste(1) .NE. zero) ste(1) = SQRT(SUM((pijh_du-pijh_dv)**2)                     &
                                       /      ste(1))
       END IF

!      compute radial derivative at full mesh points
       CALL GradientFull(pijf_ds, pijh, nblock)
!!       CALL rad_deriv_sp(pijh, pijf_ds)
!SPH:  one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
!!       CALL getorigin(pijh, m1, 0)
       pijf_ds(:,:,1) = 2*pijh(:,:,1)*ohs
       pijh(:,:,1) = 0 

!      spectrally filtered angle derivatives
       CALL toijsp(pmnch, pijh_du, 1, 0, 0, ihalf)                     ! dp/dtheta in real space (half mesh)
       CALL toijsp(pmnch, pijh_dv, 0, 1, 0, ifull)                     ! dp/dzeta in real space (half mesh, include js=1 extrap)
       DEALLOCATE(pmnch)

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL toijsp(pmnsh, tempij, fourier_du, fourier_none,                    &
                      fourier_parity_odd, fourier_half_grid)
          pijh_du = pijh_du + tempij

          CALL toijsp(pmnsh, tempij, fourier_none, fourier_du,                    &
                      fourier_parity_odd, fourier_half_grid)
          pijh_dv = pijh_dv + tempij

          DEALLOCATE(tempij)
          DEALLOCATE(pmnsh)
       END IF
#endif

!SPH   AVERAGE pijh_du,dv to FULL mesh and set origin value (SIN SYMMETRY)
       CALL to_full_mesh(pijh_du, pijf_du, nblock)
       pijf_du(:,:,ns) = pijh_du(:,:,ns)
       pijf_du(:,:,1) = pijh_du(:,:,2)        !!This is zero, but mult by v^u is finite for m=1

       CALL to_full_mesh(pijh_dv, pijf_dv, nblock)
       pijf_dv(:,:,ns) = pijh_dv(:,:,ns)      ! evaluate edge forces at Ns - 1/2 point
       pijf_dv(:,:,1) = pijh_dv(:,:,1)

       pijh(:,:,1) = 0 

       DEALLOCATE(pijh_du, pijh_dv)                 

       END SUBROUTINE init_state

       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE quantities
       USE shared_data, ONLY: l_par_state
       USE nscalingtools, ONLY: startglobrow, endglobrow
       IMPLICIT NONE
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
          DEALLOCATE (jvsupsijf, jvsupuijf, jvsupvijf,               &
                      bsupsijf0, bsupuijf0, bsupvijf0,               &
                      bsupsijh0, bsupuijh0, bsupvijh0,               &
                      bsubsijf,  bsubuijf,  bsubvijf,                &
                      ksupsijf0, ksupuijf0, ksupvijf0,               &
                      ksubsijf,  ksubuijf,  ksubvijf,                &
                      pijf, pijh, pijf_ds, pijf_du, pijf_dv,         &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'

          END IF
       END IF

       IF (lrealloc) THEN
          n3 = ns

          ALLOCATE (jvsupsijf(n1,n2,n3), jvsupuijf(n1,n2,n3),          &  ! Full mesh quantities (real space)
                    jvsupvijf(n1,n2,n3))                                  ! V_s (cosine), V_u (sine), V_v (sine)
          IF (lpar) THEN
          ALLOCATE (ksupsijf0(n1,n2,nsmin1:nsmax1),                    &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    ksupuijf0(n1,n2,nsmin1:nsmax1),                    &
                    ksupvijf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupsijf0(n1,n2,nsmin1:nsmax1),                    &  ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupvijf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupsijh0(n1,n2,nsmin1:nsmax1),                    &  ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijh0(n1,n2,nsmin1:nsmax1),                    &
                    bsupvijh0(n1,n2,nsmin1:nsmax1),                    &
                    ksubsijf(n1,n2,nsmin1:nsmax),                      &  ! Full mesh quantities (real space)
                    ksubuijf(n1,n2,nsmin1:nsmax),                      &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    ksubvijf(n1,n2,nsmin1:nsmax),                      &
                    bsubsijf(n1,n2,nsmin1:nsmax),                      &  ! Full mesh quantities (real space)
                    bsubuijf(n1,n2,nsmin1:nsmax),                      &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    bsubvijf(n1,n2,nsmin1:nsmax),                      &
                    pijh(n1,n2,nsmin1:nsmax1),                         &
                    pijf(n1,n2,nsmin1:nsmax1),                         &
                    pijf_ds(n1,n2,nsmin1:nsmax1),                      &
                    pijf_du(n1,n2,nsmin1:nsmax1),                      &
                    pijf_dv(n1,n2,nsmin1:nsmax1), stat=istat)
          ELSE
          ALLOCATE (ksupsijf0(n1,n2,n3), ksupuijf0(n1,n2,n3),          &
                    ksupvijf0(n1,n2,n3),                               &
                    ksubsijf(n1,n2,n3),                                &
                    ksubuijf(n1,n2,n3),                                &
                    ksubvijf(n1,n2,n3),                                &
                    bsubsijf(n1,n2,n3),                                &
                    bsubuijf(n1,n2,n3),                                &
                    bsubvijf(n1,n2,n3),                                &
                    pijf(n1,n2,n3), pijh(n1,n2,n3),                    &
                    pijf_ds(n1,n2,n3), pijf_du(n1,n2,n3),              &
                    pijf_dv(n1,n2,n3), bsupsijf0(n1,n2,n3),            &
                    bsupuijf0(n1,n2,n3),bsupvijf0(n1,n2,n3),           &
                    bsupsijh0(n1,n2,n3),                               &
                    bsupuijh0(n1,n2,n3),bsupvijh0(n1,n2,n3),           &
                    stat=istat)
          END IF
          IF (istat .EQ. 0) THEN

             jvsupsijf = 0; jvsupuijf = 0; jvsupvijf = 0                ! Need in add_resistivity loop
          ELSE
             STOP 'Allocation error in Allocate_Arrays'
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
            DEALLOCATE (ksupsmnsf, ksupumncf, ksupvmncf,                  &
                        djbsupsmnsh, djbsupumnch, djbsupvmnch,            &
                        djpmnch, stat=istat)
#if defined(MRC)
            IF (l_force_asym) THEN
               DEALLOCATE (ksupsmncf, ksupumnsf, ksupvmnsf,                  &
                           djbsupsmnch, djbsupumnsh, djbsupvmnsh,            &
                           djpmnsh, stat=istat)
            END IF
#endif
            IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
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
#if defined(MRC)
             IF (l_force_asym) THEN
                ALLOCATE (ksupsmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                          ksupumnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                          ksupvmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                          djpmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),          &
                          djbsupsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                          djbsupumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                          djbsupvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                          stat=istat)
             END IF
#endif
          ELSE
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,n3),                 &
                       ksupumncf(0:mpol,-ntor:ntor,n3),                 &
                       ksupvmncf(0:mpol,-ntor:ntor,n3),                 &
                       djpmnch(0:mpol,-ntor:ntor,n3),                   &
                       djbsupsmnsh(0:mpol,-ntor:ntor,n3),               &
                       djbsupumnch(0:mpol,-ntor:ntor,n3),               &
                       djbsupvmnch(0:mpol,-ntor:ntor,n3),               &
                       stat=istat)
#if defined(MRC)
             IF (l_force_asym) THEN
                ALLOCATE (ksupsmncf(0:mpol,-ntor:ntor,n3),        &
                          ksupumnsf(0:mpol,-ntor:ntor,n3),        &
                          ksupvmnsf(0:mpol,-ntor:ntor,n3),        &
                          djpmnsh(0:mpol,-ntor:ntor,n3),          &
                          djbsupsmnch(0:mpol,-ntor:ntor,n3),      &
                          djbsupumnsh(0:mpol,-ntor:ntor,n3),      &
                          djbsupvmnsh(0:mpol,-ntor:ntor,n3),      &
                          stat=istat)
             END IF
#endif
          END IF
          IF (istat .EQ. 0) THEN
             djpmnch = 0; djbsupsmnsh = 0; djbsupumnch = 0; djbsupvmnch = 0
#if defined(MRC)
             IF (l_force_asym) THEN
                djpmnsh = 0; djbsupsmnch = 0; djbsupumnsh = 0; djbsupvmnsh = 0
             END IF
#endif
          ELSE
             STOP 'Allocation error in Allocate_Arrays'
          END IF
       END IF


       END SUBROUTINE Init_Allocate_Arrays
