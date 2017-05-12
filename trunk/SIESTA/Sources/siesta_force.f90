!*******************************************************************************
!>  \file siesta_force.f90
!>  \brief Contains the \ref siesta_force module.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  \brief Update the JxB - Grad(p) force balance.
!>  \author S. P. Hirshman and R. Sanchez
!>  \date August 21, 2006
!>
!>  The MHD is in equilibrium when the force of the magnetic pressure balances
!>  the plasma pressure.
!>
!>     F = J X B - Grad(p) = 0                                               (1)
!>
!>  where J is the current density. Updates for B and p where computed in
!>  \ref siesta_bfield and \ref siesta_pressure.
!>
!>     B_i+1 = B_i + dB                                                     (2a)
!>
!>     p_i+1 = p_i + dp                                                     (2b)
!>
!>  These values must be updated first to update the MHD force. Note that the
!>  magnetic field and pressure terms contain a Jacobian term in them. New
!>  current density term can be found from the updated magnetic field.
!>
!>     J = Curl(B)/mu0                                                       (3)
!>
!>  Actual computation of the contravariant current components is implimented in
!>  \ref cv_currents.f90. The curvilinear coordinates, the J X B term becomes
!>
!>     (J X B)_s = (K^u*B^v - K^v*B^u)                                     (4a)
!>
!>     (J X B)_u = (K^v*B^s - K^s*B^v)                                     (4b)
!>
!>     (J X B)_v = (K^s*B^u - K^u*B^s)                                     (4c)
!>
!>  where K^i = sqrt(g) J^i.
!*******************************************************************************
      MODULE siesta_force

      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i, hs_i, dnorm_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd
      USE hessian, ONLY: l_Compute_Hessian

      CONTAINS
!#if defined(SKS)
!>  \brief Parallel subroutine to update MHD force
      SUBROUTINE update_force_par
!     
!     PURPOSE: UPDATES FORCE using the advanced values of B and p obtained
!              from calls to UPDATE_PRES and UPDATE_BFIELD
!              Updates values of FSUB?MNF (on the FULL radial mesh) and stores in QUANTITIES MODULE
!
!              Linearized force is Flin ~ delta_v has ONLY the linear part (not DC part and not non-linear)
!
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, l
       INTEGER, PARAMETER :: m1=1, m0=0
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         pmncf_ds, pmncf, pmnch, jpmnch_i,                               &
         bsupsmnsh_i, bsupumnch_i, bsupvmnch_i,                          &
         pmnsf_ds, pmnsf, pmnsh, jpmnsh_i,                               &
         bsupsmnch_i, bsupumnsh_i, bsupvmnsh_i,                          &
         ksupsijf, ksupuijf, ksupvijf,                                   &
         bsupsijf1,bsupuijf1,bsupvijf1,                                  &
         KxBsmnc, KxBumns, KxBvmns,                                      &
         KxBsmns, KxBumnc, KxBvmnc

       REAL(dp), POINTER, DIMENSION(:,:,:) :: KxBsij, KxBuij,            &
                                              KxBvij, jpijh
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::                &
         bsupsijh, bsupuijh, bsupvijh
       REAL(dp)    :: ton, toff, skston, skstoff, xfact, wp0, pres00_edge
       INTEGER     :: nsmin, nsmax, ns_span, nblock
!-----------------------------------------------
       CALL second0(skston)
       nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow+1)
       ns_span = nsmax-nsmin+1
!
!     Allocate space for TOTAL values of BSUPXMN?H_i and pmnch_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments are calculated in UPDATE_BFIELD and UPDATE_PRES
!     
       ALLOCATE(bsupsmnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                bsupumnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                bsupvmnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &  
                jpmnch_i   (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbed fields here (add unperturbed force gc0 in funct_island)
!
       xfact = 1
       IF (l_linearize) xfact = 0
       bsupsmnsh_i = xfact*jbsupsmnsh(:,:,nsmin:nsmax) + djbsupsmnsh
       bsupumnch_i = xfact*jbsupumnch(:,:,nsmin:nsmax) + djbsupumnch
       bsupvmnch_i = xfact*jbsupvmnch(:,:,nsmin:nsmax) + djbsupvmnch
       jpmnch_i    = xfact*jpmnch(:,:,nsmin:nsmax)     + djpmnch

!SPH:  Convert to real contravariant components of B

       ALLOCATE(bsupsijh(ntheta,nzeta,nsmin:nsmax),                      &
                bsupuijh(ntheta,nzeta,nsmin:nsmax),                      &
                bsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE_PAR'

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

       CALL toijsp_par(bsupsmnsh_i, bsupsijh, f_none, f_sin, f_full, 1, ns_span)   ! Calculate BSUP on half mesh
       CALL toijsp_par(bsupumnch_i, bsupuijh, f_none, f_cos, f_full, 1, ns_span)
       CALL toijsp_par(bsupvmnch_i, bsupvijh, f_none, f_cos, f_full, 1, ns_span)

       DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i, stat=istat)

       IF (.NOT.l_asym) GOTO 1000

       ALLOCATE(bsupsmnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                bsupumnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                bsupvmnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                jpmnsh_i   (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

! bsup*mn*h terms contain the jacobian.
       bsupsmnch_i = xfact*jbsupsmnch(:,:,nsmin:nsmax) + djbsupsmnch
       bsupumnsh_i = xfact*jbsupumnsh(:,:,nsmin:nsmax) + djbsupumnsh
       bsupvmnsh_i = xfact*jbsupvmnsh(:,:,nsmin:nsmax) + djbsupvmnsh
       jpmnsh_i    = xfact*jpmnsh(:,:,nsmin:nsmax)     + djpmnsh

! bsup*ijh terms still contain the jacobian.
       CALL toijsp_par(bsupsmnch_i, bsupsijh, f_sum, f_cos, f_full, 1, ns_span)
       CALL toijsp_par(bsupumnsh_i, bsupuijh, f_sum, f_sin, f_full, 1, ns_span)
       CALL toijsp_par(bsupvmnsh_i, bsupvijh, f_sum, f_sin, f_full, 1, ns_span)

       DEALLOCATE(bsupsmnch_i, bsupumnsh_i, bsupvmnsh_i, stat=istat)

 1000  CONTINUE

!
!      Calculate contravariant components of J = curl B (Ksup = gsqrt*J)
!      and update magnetic energy if nonlinear force is being computed
!                             
       ALLOCATE(ksupsijf(ntheta,nzeta,nsmin:nsmax),                      & ! Full mesh quantities (real space)
                ksupuijf(ntheta,nzeta,nsmin:nsmax),                      & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijf(ntheta,nzeta,nsmin:nsmax),                      &
                bsupsijf1(ntheta,nzeta,nsmin:nsmax),                     &
                bsupuijf1(ntheta,nzeta,nsmin:nsmax),                     &
                bsupvijf1(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat.ne.0) STOP 'Allocation failed in UPDATE_FORCE'

!
!      Calculate updated FULL-MESH components of BSUPX in real-space by averaging
!
       bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
       bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
       CALL bhalftobfull_par (bsupsijh, bsupuijh, bsupvijh,              &
                              bsupsijf1,bsupuijf1,bsupvijf1,             &
                              nsmin, nsmax)

!Origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents_par(bsupsijh, bsupuijh, bsupvijh,                &
                            ksupsijf, ksupuijf, ksupvijf,                &
                           (.NOT.l_linearize .OR. l_getwmhd), .FALSE.)

!      alias work arrays for use as temp work space
       KxBsij=>bsupsijh;  KxBuij=>bsupuijh;  KxBvij=>bsupvijh

       nsmax=MIN(endglobrow,ns)
       ns_span = nsmax-nsmin+1

       IF (l_linearize) THEN
!      dJ X B0 (note here, ksupXij?f is the perturbed current)
          KxBsij(:,:,nsmin:nsmax) =  bsupvijf0(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax) - &
            bsupuijf0(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax)
          KxBuij(:,:,nsmin:nsmax) =  bsupsijf0(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax) - &
            bsupvijf0(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax)
          KxBvij(:,:,nsmin:nsmax) =  bsupuijf0(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax) - &
            bsupsijf0(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax)

!      J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!      Now, bsupXij?f1 = perturbed B-field

          KxBsij(:,:,nsmin:nsmax) = KxBsij(:,:,nsmin:nsmax) +            &
              (bsupvijf1(:,:,nsmin:nsmax)*ksupuijf0(:,:,nsmin:nsmax) -   &
               bsupuijf1(:,:,nsmin:nsmax)*ksupvijf0(:,:,nsmin:nsmax))
          KxBuij(:,:,nsmin:nsmax) = KxBuij(:,:,nsmin:nsmax) +            &
              (bsupsijf1(:,:,nsmin:nsmax)*ksupvijf0(:,:,nsmin:nsmax) -   &
               bsupvijf1(:,:,nsmin:nsmax)*ksupsijf0(:,:,nsmin:nsmax))
          KxBvij(:,:,nsmin:nsmax) = KxBvij(:,:,nsmin:nsmax) +            &
              (bsupuijf1(:,:,nsmin:nsmax)*ksupsijf0(:,:,nsmin:nsmax) -   &
               bsupsijf1(:,:,nsmin:nsmax)*ksupuijf0(:,:,nsmin:nsmax))

       ELSE

          KxBsij(:,:,nsmin:nsmax) = bsupvijf1(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax) - &
             bsupuijf1(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax)           ! WORK1 has cosine parity (IPARITY = 0)
          KxBuij(:,:,nsmin:nsmax) = bsupsijf1(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax) - &
             bsupvijf1(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax)           ! WORK2 has sine parity (IPARITY = 1)
          KxBvij(:,:,nsmin:nsmax) = bsupuijf1(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax) - &
             bsupsijf1(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax)           ! WORK3 has sine parity (IPARITY = 1)
       END IF

       IF (nsmax .EQ. ns) THEN
          IF (ANY(bsupsijf1(:,:,ns) .NE. zero)) STOP 'bsupsijf1 != 0'
       END IF

       DEALLOCATE(ksupsijf, ksupuijf, ksupvijf,                          &
                  bsupsijf1,bsupuijf1,bsupvijf1, stat=istat)

       ALLOCATE(KxBsmnc(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                KxBumns(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                KxBvmns(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'
       
       CALL tomnsp_par(KxBsij, KxBsmnc, f_cos, f_full, 1, ns_span)          ! Harmonics corresponding to..
       CALL tomnsp_par(KxBuij, KxBumns, f_sin, f_full, 1, ns_span)          ! Lorentz Force on full grid
       CALL tomnsp_par(KxBvij, KxBvmns, f_sin, f_full, 1, ns_span)

       IF (l_asym) THEN
       ALLOCATE(KxBsmns(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                KxBumnc(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                KxBvmnc(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL tomnsp_par(KxBsij, KxBsmns, f_sin, f_full, 1, ns_span)          ! Harmonics corresponding to..
       CALL tomnsp_par(KxBuij, KxBumnc, f_cos, f_full, 1, ns_span)          ! Lorentz Force on full grid
       CALL tomnsp_par(KxBvij, KxBvmnc, f_cos, f_full, 1, ns_span)
       END IF

       DEALLOCATE(bsupuijh, bsupvijh, stat=istat)                           ! aliased arrays now

       nsmax=MIN(endglobrow+1,ns)
       ns_span = nsmax-nsmin+1
!
!      Calculate radial derivative of pressure gradient 
!
       jpijh => bsupsijh
       ALLOCATE(pmncf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                pmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL toijsp_par(jpmnch_i, jpijh, f_none, f_cos, f_half, 1, ns_span)  ! Calculate jac*pressure in real space (half mesh)
       DEALLOCATE (jpmnch_i, stat=istat)
 
       IF (l_asym) THEN
       ALLOCATE(pmnsf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
                pmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL toijsp_par(jpmnsh_i, jpijh, f_sum, f_sin, f_half, 1, ns_span)   ! Calculate jac*pressure in real space (half mesh)
       DEALLOCATE (jpmnsh_i, stat=istat)
       END IF
!
!      UPDATE THERMAL (PRESSURE) ENERGY
!
       IF (l_getwmhd) THEN
          IF (INHESSIAN) STOP 'l_getwmhd must be set to FALSE in Hessian'
          nsmax=MIN(endglobrow,ns)
          wp = 0
          DO l = 1, nu_i
             wp = wp + SUM(jpijh(l,:,MAX(2,nsmin):nsmax))*cosmui(l,m0)
          END DO
          wp0 = signjac*(twopi*twopi * hs_i)*wp
          nsmax=MIN(endglobrow+1,ns)
#if defined(MPI_OPT)
          CALL MPI_ALLREDUCE(wp0,wp,1,MPI_REAL8, MPI_SUM,                &
                             MPI_COMM_WORLD,MPI_ERR)
#else
          wp = wp0
#endif
       END IF

!SPH   NOW CONVERT THIS TO REAL PRESSURE
!SPH   Evolved pressure is jacobh * real-pressure
       jpijh = jpijh/jacobh(:,:,nsmin:nsmax)
       ALLOCATE(pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL tomnsp_par(jpijh, pmnch, f_cos, f_half, 1, ns_span)

       IF (l_asym) THEN
       ALLOCATE(pmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL tomnsp_par(jpijh, pmnsh, f_sin, f_half, 1, ns_span)
       END IF

       DEALLOCATE(bsupsijh, stat=istat)
       
       nblock=SIZE(pmnch,1)*SIZE(pmnch,2)
       CALL to_full_mesh_par(pmnch, pmncf, nblock, nsmin, nsmax)
!SPH11-7-16: s=1 is an iso-pressure contour
       IF (nsmax .EQ. ns) THEN
          pres00_edge = pmncf(0,0,ns)
          pmncf(:,:,ns) = 0;   pmncf(0,0,ns) = pres00_edge
       END IF
       CALL GradientFull_par(pmncf_ds, pmnch, nblock, nsmin, nsmax)

!
!      Calculate harmonics of MHD force (F = JxB - grad p) on full mesh
!      For l_linearize=T:  PERTURBED F (ignore F(xc=0) part)
!      For l_linearize=F:  NONLINEAR F


!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
       IF (nsmin .EQ. 1) THEN
          pmncf(:,:,1) = 0
          pmncf(m0:m1,:,1) = pmnch(m0:m1,:,2)
          pmncf_ds(m1,:,1) = 2*ohs*pmnch(m1,:,2)
       END IF

       IF (l_asym) THEN
       CALL to_full_mesh_par(pmnsh, pmnsf, nblock, nsmin, nsmax)
!SPH11-7-16: s=1 is an iso-pressure contour
       IF (nsmax .EQ. ns) THEN
          pmnsf(:,:,ns) = 0
       END IF
       CALL GradientFull_par(pmnsf_ds, pmnsh, nblock, nsmin, nsmax)

       IF (nsmin .EQ. 1) THEN
           pmnsf(:,:,1) = 0
           pmnsf(m0:m1,:,1) = pmnsh(m0:m1,:,2)
           pmnsf_ds(m1,:,1) = 2*ohs*pmnsh(m1,:,2)
       END IF
       END IF

       nsmax=MIN(endglobrow,ns)
       KxBsmnc = KxBsmnc - pmncf_ds(:,:,nsmin:nsmax)
       CALL get_force_harmonics_par(pmncf, KxBsmnc, KxBumns, KxBvmns,    &
                                    fsubsmncf, fsubumnsf, fsubvmnsf, 0)
       DEALLOCATE(KxBsmnc, KxBumns, KxBvmns, pmncf, pmncf_ds, pmnch,     &
                  stat=istat)

       IF (l_asym) THEN
       KxBsmns = KxBsmns - pmnsf_ds(:,:,nsmin:nsmax)
       CALL get_force_harmonics_par(pmnsf, KxBsmns, KxBumnc, KxBvmnc,    &
                                    fsubsmnsf, fsubumncf, fsubvmncf, 1)
       DEALLOCATE(KxBsmns, KxBumnc, KxBvmnc, pmnsf, pmnsf_ds, pmnsh,     &
                  stat=istat)
       END IF

       CALL second0(skstoff)
       time_update_force = time_update_force + (skstoff-skston)

      END SUBROUTINE update_force_par
!#endif

!>  \brief Serial subroutine to update MHD force
      SUBROUTINE update_force

      CALL update_force_par

      END SUBROUTINE update_force

      END MODULE siesta_force
