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
!>     F = K X B - Grad(p) = 0                                               (1)
!>
!>  where K is the current density. Updates for B and p where computed in
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
!>     K = Curl(B)/mu0                                                       (3)
!>
!>  Actual computation of the contravariant current components is implimented in 
!>  \ref cv_currents.f90. The curvilinear coordinates, the J X B term becomes
!>
!>     (K X B)_s = J(K^u*B^v - K^v*B^u)                                     (4a)
!>
!>     (K X B)_u = J(K^v*B^s - K^s*B^v)                                     (4b)
!>
!>     (K X B)_v = J(K^s*B^u - K^u*B^s)                                     (4c)
!>
!*******************************************************************************
      MODULE siesta_force

      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i, hs_i, dnorm_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd, xPosDef
      USE hessian, ONLY: l_Compute_Hessian

      CONTAINS
#if defined(SKS)
!*******************************************************************************
!>  \brief Parallel subroutine to update MHD force
!*******************************************************************************
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
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m1=1, m0=0, pcos=0, psin=1
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    &
         pmncf_ds, pmncf, pmnch, jpmnch_i,                              &
         bsupsmnsh_i, bsupumnch_i, bsupvmnch_i,                         &
#if defined(MRC)
         pmnsf_ds, pmnsf, pmnsh, jpmnsh_i,                              &
         bsupsmnch_i, bsupumnsh_i, bsupvmnsh_i,                         &
#endif
         ksupsijf, ksupuijf, ksupvijf,                                  &
         bsupsijf1,bsupuijf1,bsupvijf1,                                 &
#if defined(MRC)
         jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,                         &
#endif
         jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc, tempij

       REAL(rprec), POINTER, DIMENSION(:,:,:) :: jkxbsubsij,            &
                                                 jkxbsubuij,            &
                                                 jkxbsubvij,            &
                                                 jpijh
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::            &
         bsupsijh, bsupuijh, bsupvijh
       REAL(rprec) :: ton, toff, skston, skstoff, xfact, wp0
       INTEGER     :: nsmin, nsmax, ns_span, nblock
!-----------------------------------------------
       CALL second0(skston)
       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       ns_span = nsmax-nsmin+1
!
!     Allocate space for TOTAL values of BSUPXMN?H_i and pmnch_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments are calculated in UPDATE_BFIELD and UPDATE_PRES
!     
       ALLOCATE(bsupsmnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                bsupumnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                bsupvmnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &  
                jpmnch_i   (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbation here (add gc0 back in funct_island)
!
       xfact = 1
       IF (l_linearize) xfact = 0
       bsupsmnsh_i = xfact*jbsupsmnsh(:,:,nsmin:nsmax) + djbsupsmnsh
       bsupumnch_i = xfact*jbsupumnch(:,:,nsmin:nsmax) + djbsupumnch
       bsupvmnch_i = xfact*jbsupvmnch(:,:,nsmin:nsmax) + djbsupvmnch
       jpmnch_i    = xfact*jpmnch(:,:,nsmin:nsmax)     + djpmnch

!SPH:  Convert to real covariant components

       ALLOCATE(bsupsijh(ntheta,nzeta,nsmin:nsmax),                   &
                bsupuijh(ntheta,nzeta,nsmin:nsmax),                    &
                bsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE_PAR'

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

       CALL toijsp_par(bsupsmnsh_i, bsupsijh, 0, 0, psin, ifull, 1, ns_span)             ! Calculate BSUP on half mesh
       CALL toijsp_par(bsupumnch_i, bsupuijh, 0, 0, pcos, ifull, 1, ns_span)
       CALL toijsp_par(bsupvmnch_i, bsupvijh, 0, 0, pcos, ifull, 1, ns_span)

       DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i, stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(bsupsmnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                   bsupumnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                   bsupvmnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                   jpmnsh_i   (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          bsupsmnch_i = xfact*jbsupsmnch(:,:,nsmin:nsmax) + djbsupsmnch
          bsupumnsh_i = xfact*jbsupumnsh(:,:,nsmin:nsmax) + djbsupumnsh
          bsupvmnsh_i = xfact*jbsupvmnsh(:,:,nsmin:nsmax) + djbsupvmnsh
          jpmnsh_i    = xfact*jpmnsh(:,:,nsmin:nsmax)     + djpmnsh

          ALLOCATE(tempij(ntheta,nzeta,nsmin:nsmax))

          CALL toijsp_par(bsupsmnch_i, tempij, fourier_none, fourier_none, fourier_parity_even, fourier_full_grid, 1, ns_span)
          bsupsijh = bsupsijh + tempij

          CALL toijsp_par(bsupumnsh_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid, 1, ns_span)
          bsupuijh = bsupuijh + tempij

          CALL toijsp_par(bsupvmnsh_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid, 1, ns_span)
          bsupvijh = bsupvijh + tempij

          DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i)

       END IF
#endif

       bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
       bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)

!
!      Calculate contravariant components of J = curl B (Ksup = gsqrt*J)
!      and update magnetic energy if nonlinear force is being computed
!                             
       ALLOCATE(ksupsijf(ntheta,nzeta,nsmin:nsmax),                     & ! Full mesh quantities (real space)
                ksupuijf(ntheta,nzeta,nsmin:nsmax),                     & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijf(ntheta,nzeta,nsmin:nsmax),                     &
                bsupsijf1(ntheta,nzeta,nsmin:nsmax),                    &
                bsupuijf1(ntheta,nzeta,nsmin:nsmax),                    &
                bsupvijf1(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat.ne.0) STOP 'Allocation failed in UPDATE_FORCE'

!
!      Calculate updated FULL-MESH components of BSUPX in real-space by averaging
!
       CALL bhalftobfull_par (bsupsijh, bsupuijh, bsupvijh,          &
                              bsupsijf1,bsupuijf1,bsupvijf1,         &
                              nsmin, nsmax)

!Origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents_par(bsupsijh, bsupuijh, bsupvijh,            &
                            ksupsijf, ksupuijf, ksupvijf,            &
                           (.NOT.l_linearize .OR. l_getwmhd), .FALSE.)

!      alias work arrays for use as temp work space
       jkxbsubsij=>bsupsijh
       jkxbsubuij=>bsupuijh
       jkxbsubvij=>bsupvijh

       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
       ns_span = nsmax-nsmin+1

       IF (l_linearize) THEN
!      dJ X B0 (note here, ksupXij?f is the perturbed current)
          jkxbsubsij(:,:,nsmin:nsmax) =  bsupvijf0(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax) - & ! J(K^u*B^v - K^v*B^u)
            bsupuijf0(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax)
          jkxbsubuij(:,:,nsmin:nsmax) =  bsupsijf0(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax) - & ! J(K^v*B^s - K^s*B^v)
            bsupvijf0(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax)
          jkxbsubvij(:,:,nsmin:nsmax) =  bsupuijf0(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax) - & ! J(K^s*B^u - K^u*B^s)
            bsupsijf0(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax)

!      J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!      Now, bsupXij?f1 = perturbed B-field
          xfact = 1
!!        IF (l_Compute_Hessian) xfact = xfact-xPosDef
          xfact = xfact-xPosDef

          jkxbsubsij(:,:,nsmin:nsmax) = jkxbsubsij(:,:,nsmin:nsmax) + xfact*                &
                     (bsupvijf1(:,:,nsmin:nsmax)*ksupuijf0(:,:,nsmin:nsmax) -   &
                      bsupuijf1(:,:,nsmin:nsmax)*ksupvijf0(:,:,nsmin:nsmax))
          jkxbsubuij(:,:,nsmin:nsmax) = jkxbsubuij(:,:,nsmin:nsmax) + xfact*                &
                     (bsupsijf1(:,:,nsmin:nsmax)*ksupvijf0(:,:,nsmin:nsmax) -   &
                      bsupvijf1(:,:,nsmin:nsmax)*ksupsijf0(:,:,nsmin:nsmax))
          jkxbsubvij(:,:,nsmin:nsmax) = jkxbsubvij(:,:,nsmin:nsmax) + xfact*                &
                     (bsupuijf1(:,:,nsmin:nsmax)*ksupsijf0(:,:,nsmin:nsmax) -   &
                      bsupsijf1(:,:,nsmin:nsmax)*ksupuijf0(:,:,nsmin:nsmax))

       ELSE

          jkxbsubsij(:,:,nsmin:nsmax) = bsupvijf1(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax) - & ! J(K^u*B^v - K^v*B^u)
             bsupuijf1(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax)
          jkxbsubuij(:,:,nsmin:nsmax) = bsupsijf1(:,:,nsmin:nsmax)*ksupvijf(:,:,nsmin:nsmax) - & ! J(K^v*B^s - K^s*B^v)
             bsupvijf1(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax)
          jkxbsubvij(:,:,nsmin:nsmax) = bsupuijf1(:,:,nsmin:nsmax)*ksupsijf(:,:,nsmin:nsmax) - & ! J(K^s*B^u - K^u*B^s)
             bsupsijf1(:,:,nsmin:nsmax)*ksupuijf(:,:,nsmin:nsmax)

       END IF

       DEALLOCATE(ksupsijf, ksupuijf, ksupvijf,                      &
                  bsupsijf1,bsupuijf1,bsupvijf1, stat=istat)

       ALLOCATE(jkxbsubsmnc(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                jkxbsubumns(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                jkxbsubvmns(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       
       CALL tomnsp_par(jkxbsubsij, jkxbsubsmnc, pcos, ifull, 1, ns_span)             ! Harmonics corresponding to..
       CALL tomnsp_par(jkxbsubuij, jkxbsubumns, psin, ifull, 1, ns_span)             ! Lorentz Force on full grid
       CALL tomnsp_par(jkxbsubvij, jkxbsubvmns, psin, ifull, 1, ns_span)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(jkxbsubsmns(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                   jkxbsubumnc(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                   jkxbsubvmnc(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL tomnsp_par(jkxbsubsij, jkxbsubsmns, fourier_parity_odd, fourier_full_grid, 1, ns_span)             ! Harmonics corresponding to..
          CALL tomnsp_par(jkxbsubuij, jkxbsubumnc, fourier_parity_even, fourier_full_grid, 1, ns_span)             ! Lorentz Force on full grid
          CALL tomnsp_par(jkxbsubvij, jkxbsubvmnc, fourier_parity_even, fourier_full_grid, 1, ns_span)
       END IF
#endif

       DEALLOCATE(bsupuijh, bsupvijh, stat=istat)                       ! aliased work2-3 deallocated now

       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       ns_span = nsmax-nsmin+1
!
!      Calculate radial derivative of pressure gradient 
!
       jpijh => bsupsijh

       ALLOCATE(pmncf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                pmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL toijsp_par(jpmnch_i, jpijh, 0, 0, pcos, ihalf, 1, ns_span)                    ! Calculate jac*pressure in real space (half mesh)

       DEALLOCATE (jpmnch_i, stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                   pmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL toijsp_par(jpmnsh_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_half_grid, 1, ns_span)                    ! Calculate jac*pressure in real space (half mesh)
          jpijh = jpijh + tempij

          DEALLOCATE (jpmnsh_i, stat=istat)
          DEALLOCATE (tempij)
       END IF
#endif

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
          wp0 = 2*signjac*(twopi*pi * hs_i)*wp                            !Factor of 2 if 0-pi only in theta
          nsmax=MIN(endglobrow+1,ns)
          CALL MPI_ALLREDUCE(wp0,wp,1,MPI_REAL8, MPI_SUM,               &
                             MPI_COMM_WORLD,MPI_ERR)
       END IF

!SPH   NOW CONVERT THIS TO REAL PRESSURE
!SPH   Evolved pressure is jacobh * real-pressure
       jpijh = jpijh/jacobh(:,:,nsmin:nsmax)

       ALLOCATE(pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL tomnsp_par(jpijh, pmnch, pcos, ihalf, 1, ns_span)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL tomnsp_par(jpijh, pmnsh, fourier_parity_odd, fourier_half_grid, 1, ns_span)
       END IF
#endif

       DEALLOCATE(bsupsijh, stat=istat)
       
       nblock=SIZE(pmnch,1)*SIZE(pmnch,2)
       CALL to_full_mesh_par(pmnch, pmncf, nblock, nsmin, nsmax)
       CALL GradientFull_par(pmncf_ds, pmnch, nblock, nsmin, nsmax)

!
!      Calculate harmonics of net force (i.e., F = JxB - grad p) on full mesh
!      THIS IS THE TRUE NONLINEAR JXB-grad(p) FORCE when l_linearize=FALSE
!      AND THE PERTURBED FORCE WHEN l_linearize=TRUE

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
       IF (nsmin .EQ. 1) THEN
          pmncf(:,:,1) = 0
          pmncf(m0:m1,:,1) = pmnch(m0:m1,:,2)
 !         pmncf_ds(:,:,1) = 0
          pmncf_ds(m1,:,1) = 2*ohs*pmnch(m1,:,2)
       END IF

       nsmax=MIN(endglobrow,ns)
       jkxbsubsmnc = jkxbsubsmnc - pmncf_ds(:,:,nsmin:nsmax)

       CALL get_force_harmonics_par(pmncf, jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,   &
                                    fsubsmncf, fsubumnsf, fsubvmnsf)

       DEALLOCATE(jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,                &
                  pmncf, pmncf_ds, pmnch, stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL to_full_mesh_par(pmnsh, pmnsf, nblock, nsmin, nsmax)
          CALL GradientFull_par(pmnsf_ds, pmnch, nblock, nsmin, nsmax)

          IF (nsmin .EQ. 1) THEN
             pmnsf(:,:,1) = 0
             pmnsf(m0:m1,:,1) = pmnsh(m0:m1,:,2)
             pmnsf_ds(m1,:,1) = 2*ohs*pmnsh(m1,:,2)
          END IF

          jkxbsubsmns = jkxbsubsmns - pmnsf_ds(:,:,nsmin:nsmax)

          CALL get_force_harmonics_par_other(pmnsf, jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc,   &
                                             fsubsmnsf, fsubumncf, fsubvmncf)

          DEALLOCATE(jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc,                &
                     pmnsf, pmnsf_ds, pmnsh, stat=istat)
       END IF
#endif

       CALL second0(skstoff)
       time_update_force = time_update_force + (skstoff-skston)

      END SUBROUTINE update_force_par
#endif

!*******************************************************************************
!>  \brief Serial subroutine to update MHD force
!*******************************************************************************
      SUBROUTINE update_force
#if defined(SKS)
       USE descriptor_mod, ONLY: DIAGONALDONE, INHESSIAN
#endif
       IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, l, nblock
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m1=1, m0=0, pcos=0, psin=1
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    &
         pmncf_ds, pmncf, pmnch, jpmnch_i,                              &
         bsupsmnsh_i, bsupumnch_i, bsupvmnch_i,                         &
#if defined(MRC)
         pmnsf_ds, pmnsf, pmnsh, jpmnsh_i,                              &
         bsupsmnch_i, bsupumnsh_i, bsupvmnsh_i,                         &
#endif
         ksupsijf, ksupuijf, ksupvijf,                                  &
         bsupsijf1,bsupuijf1,bsupvijf1,                                 &
#if defined(MRC)
         jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,                         &
#endif
         jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc, tempij

       REAL(rprec), POINTER, DIMENSION(:,:,:) :: jkxbsubsij,            &
                                                 jkxbsubuij,            &
                                                 jkxbsubvij,            &
                                                 jpijh
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::            &
         bsupsijh, bsupuijh, bsupvijh
       REAL(rprec) :: ton, toff, xfact, skston, skstoff 
!-----------------------------------------------
       CALL second0(ton)
       skston=ton
!
!     Allocate space for TOTAL values of BSUPXMN?H_i and pmnch_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments are calculated in UPDATE_BFIELD and UPDATE_PRES
!     
       ALLOCATE(bsupsmnsh_i(0:mpol,-ntor:ntor,ns),                      &
                bsupumnch_i(0:mpol,-ntor:ntor,ns),                      &
                bsupvmnch_i(0:mpol,-ntor:ntor,ns),                      &  
                jpmnch_i   (0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbation here (add gc0 back in funct_island)
!
       xfact = 1
       IF (l_linearize) xfact = 0
       bsupsmnsh_i = xfact*jbsupsmnsh + djbsupsmnsh
       bsupumnch_i = xfact*jbsupumnch + djbsupumnch
       bsupvmnch_i = xfact*jbsupvmnch + djbsupvmnch
       jpmnch_i    = xfact*jpmnch     + djpmnch

!SPH:  Convert to real covariant components

       ALLOCATE(bsupsijh(ntheta,nzeta,ns),                             &
                bsupuijh(ntheta,nzeta,ns),                              &
                bsupvijh(ntheta,nzeta,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

       CALL toijsp(bsupsmnsh_i, bsupsijh, 0, 0, psin, ifull)             ! Calculate BSUP on half mesh
       CALL toijsp(bsupumnch_i, bsupuijh, 0, 0, pcos, ifull)
       CALL toijsp(bsupvmnch_i, bsupvijh, 0, 0, pcos, ifull)

       DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i)

#if defined(MRC)
       IF (l_force_asym) THEN

          ALLOCATE(bsupsmnch_i(0:mpol,-ntor:ntor,ns),                      &
                   bsupumnsh_i(0:mpol,-ntor:ntor,ns),                      &
                   bsupvmnsh_i(0:mpol,-ntor:ntor,ns),                      &
                   jpmnsh_i   (0:mpol,-ntor:ntor,ns), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbation here (add gc0 back in funct_island)
!
          bsupsmnch_i = xfact*jbsupsmnch + djbsupsmnch
          bsupumnsh_i = xfact*jbsupumnsh + djbsupumnsh
          bsupvmnsh_i = xfact*jbsupvmnsh + djbsupvmnsh
          jpmnsh_i    = xfact*jpmnsh     + djpmnsh

!SPH:  Convert to real covariant components

          ALLOCATE(tempij(ntheta,nzeta,ns))

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

          CALL toijsp(bsupsmnch_i, tempij, fourier_none, fourier_none, fourier_parity_even, fourier_full_grid)             ! Calculate BSUP on half mesh
          bsupsijh = bsupsijh + tempij

          CALL toijsp(bsupumnch_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid)
          bsupuijh = bsupuijh + tempij

          CALL toijsp(bsupvmnch_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid)
          bsupvijh = bsupvijh + tempij

          DEALLOCATE(bsupsmnch_i, bsupumnsh_i, bsupvmnsh_i)

       END IF
#endif

       bsupsijh = bsupsijh/jacobh
!       bsupuijh = bsupuijh/jacobh
       bsupvijh = bsupvijh/jacobh

!
!      Calculate contravariant components of J = curl B (Ksup = gsqrt*J)
!      and update magnetic energy if nonlinear force is being computed
!                             
       ALLOCATE(ksupsijf(ntheta,nzeta,ns),                             & ! Full mesh quantities (real space)
                ksupuijf(ntheta,nzeta,ns),                              & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijf(ntheta,nzeta,ns),                              &
                bsupsijf1(ntheta,nzeta,ns),                            &
                bsupuijf1(ntheta,nzeta,ns),                             &
                bsupvijf1(ntheta,nzeta,ns), stat=istat)
       IF (istat.ne.0) STOP 'Allocation failed in UPDATE_FORCE'

!
!      Calculate updated FULL-MESH components of BSUPX in real-space
!
       CALL bhalftobfull (bsupsijh,  bsupuijh,  bsupvijh,            &
                          bsupsijf1, bsupuijf1, bsupvijf1)
       IF (ANY(bsupsijf1(:,:,ns) .NE. zero)) PRINT *,'bsupsijf1(ns) != 0'

!Origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents (bsupsijh, bsupuijh, bsupvijh,               &
                         ksupsijf, ksupuijf, ksupvijf,               &
                         .NOT.l_linearize, .FALSE.)

!      alias work arrays for use as temp work space
       jkxbsubsij=>bsupsijh
       jkxbsubuij=>bsupuijh
       jkxbsubvij=>bsupvijh

       IF (l_linearize) THEN
!      dJ X B0 (note here, ksupXij?f is the perturbed current)
          jkxbsubsij =  bsupvijf0*ksupuijf - bsupuijf0*ksupvijf ! J(K^u*B^v - K^v*B^u)
          jkxbsubuij =  bsupsijf0*ksupvijf - bsupvijf0*ksupsijf ! J(K^v*B^s - K^s*B^v)
          jkxbsubvij =  bsupuijf0*ksupsijf - bsupsijf0*ksupuijf ! J(K^s*B^u - K^u*B^s)

!      J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!      Now, bsupXij?f1 = perturbed B-field
          xfact = 1
!!        IF (l_Compute_Hessian) xfact = xfact-xPosDef
          xfact = xfact-xPosDef

          jkxbsubsij = jkxbsubsij+ xfact*                                        &
      &          (bsupvijf1*ksupuijf0 - bsupuijf1*ksupvijf0)
          jkxbsubuij = jkxbsubuij + xfact*                                        &
      &          (bsupsijf1*ksupvijf0 - bsupvijf1*ksupsijf0)
          jkxbsubvij = jkxbsubvij + xfact*                                        &
      &          (bsupuijf1*ksupsijf0 - bsupsijf1*ksupuijf0)

       ELSE

          jkxbsubsij = bsupvijf1*ksupuijf - bsupuijf1*ksupvijf ! J(K^u*B^v - K^v*B^u)
          jkxbsubuij = bsupsijf1*ksupvijf - bsupvijf1*ksupsijf ! J(K^v*B^s - K^s*B^v)
          jkxbsubvij = bsupuijf1*ksupsijf - bsupsijf1*ksupuijf ! J(K^s*B^u - K^u*B^s)

       END IF
      
       DEALLOCATE(ksupsijf, ksupuijf, ksupvijf,                      &
                  bsupsijf1,bsupuijf1,bsupvijf1)

       ALLOCATE(jkxbsubsmnc(0:mpol,-ntor:ntor,ns),                            &
                jkxbsubumns(0:mpol,-ntor:ntor,ns),                            &
                jkxbsubvmns(0:mpol,-ntor:ntor,ns), stat = istat)
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       
       CALL tomnsp(jkxbsubsij, jkxbsubsmnc, pcos, ifull)                             ! Harmonics corresponding to..
       CALL tomnsp(jkxbsubuij, jkxbsubumns, psin, ifull)                             ! Lorentz Force on full grid
       CALL tomnsp(jkxbsubvij, jkxbsubvmns, psin, ifull)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(jkxbsubsmns(0:mpol,-ntor:ntor,ns),                            &
                   jkxbsubumnc(0:mpol,-ntor:ntor,ns),                            &
                   jkxbsubvmnc(0:mpol,-ntor:ntor,ns), stat = istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL tomnsp(jkxbsubsij, jkxbsubsmns, fourier_parity_odd, fourier_full_grid)                             ! Harmonics corresponding to..
          CALL tomnsp(jkxbsubuij, jkxbsubumnc, fourier_parity_even, fourier_full_grid)                             ! Lorentz Force on full grid
          CALL tomnsp(jkxbsubvij, jkxbsubvmnc, fourier_parity_even, fourier_full_grid)
       END IF
#endif

       DEALLOCATE(bsupuijh, bsupvijh)                                   ! aliased work2-3 deallocated now

!
!      Calculate radial derivative of pressure gradient 
!
       jpijh => bsupsijh

       ALLOCATE(pmncf_ds(0:mpol,-ntor:ntor,ns),                         &
                pmncf(0:mpol,-ntor:ntor,ns), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL toijsp(jpmnch_i, jpijh, 0, 0, pcos, ihalf)                    ! Calculate jac*pressure in real space (half mesh)
       DEALLOCATE (jpmnch_i)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsf_ds(0:mpol,-ntor:ntor,ns),                         &
                   pmnsf(0:mpol,-ntor:ntor,ns), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL toijsp(jpmnsh_i, tempij, fourier_none, fourier_none, fourier_parity_odd, fourier_half_grid)
          jpijh = jpijh + tempij

          DEALLOCATE (jpmnsh_i, stat=istat)
          DEALLOCATE (tempij)
       END IF
#endif

       IF (.not.l_linearize) THEN
!
!      UPDATE THERMAL (PRESSURE) ENERGY
!
          wp = 0
          DO l = 1, nu_i
             wp = wp + SUM(jpijh(l,:,2:))*cosmui(l,m0)
          END DO

          wp = 2*signjac*(twopi*pi * hs_i)*wp                            !Factor of 2 if 0-pi only in theta

       END IF

!SPH   NOW CONVERT THIS TO REAL PRESSURE
!SPH   Evolved pressure is jacobh * real-pressure
       jpijh = jpijh/jacobh
       ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL tomnsp(jpijh, pmnch, pcos, ihalf)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE(pmnsh(0:mpol,-ntor:ntor,ns), stat=istat)
          IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

          CALL tomnsp(jpijh, pmnsh, fourier_parity_odd, fourier_half_grid)
       END IF
#endif

       DEALLOCATE(bsupsijh)

       nblock=SIZE(pmnch,1)*SIZE(pmnch,2)
       CALL to_full_mesh(pmnch, pmncf, nblock)
       CALL GradientFull(pmncf_ds, pmnch, nblock)
!
!      Calculate harmonics of net force (i.e., F = JxB - grad p) on full mesh
!      THIS IS THE TRUE NONLINEAR JXB-grad(p) FORCE when l_linearize=FALSE
!      AND THE PERTURBED FORCE WHEN l_linearize=TRUE

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
       pmncf(:,:,1) = 0
       pmncf(m0:m1,:,1) = pmnch(m0:m1,:,2)
       pmncf_ds(:,:,1) = 0
       pmncf_ds(m1,:,1) = 2*ohs*pmnch(m1,:,2)

       CALL get_force_harmonics(pmncf, pmncf_ds,                             &
                                jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,       &
                                fsubsmncf, fsubumnsf, fsubvmnsf)

       DEALLOCATE(jkxbsubsmnc, jkxbsubumns, jkxbsubvmns,                     &
                  pmncf, pmncf_ds, pmnch)

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL to_full_mesh(pmnsh, pmnsf, nblock)
          CALL GradientFull(pmnsf_ds, pmnsh, nblock)

          pmnsf(:,:,1) = 0
          pmnsf(m0:m1,:,1) = pmnsh(m0:m1,:,2)
          pmnsf_ds(:,:,1) = 0
          pmnsf_ds(m1,:,1) = 2*ohs*pmnsh(m1,:,2)

          CALL get_force_harmonics_other(pmnsf, pmnsf_ds,                      &
                                   jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc,      &
                                   fsubsmnsf, fsubumncf, fsubvmncf)

          DEALLOCATE(jkxbsubsmns, jkxbsubumnc, jkxbsubvmnc,                    &
                     pmnsf, pmnsf_ds, pmnsh, stat=istat)
       END IF
#endif

       CALL second0(toff)
       time_update_force = time_update_force + (toff-ton)

       CALL second0(skstoff)
#if defined(SKS)
       IF(DIAGONALDONE.AND.INHESSIAN) time_update_force=time_update_force+(skstoff-skston)
#endif
      END SUBROUTINE update_force

      END MODULE siesta_force
