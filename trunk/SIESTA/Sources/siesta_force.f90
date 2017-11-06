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
      USE v3_utilities, ONLY: assert
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE island_params, ONLY: nu_i, hs_i, dnorm_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd, lverbose
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE quantities
      USE mpi_inc

      IMPLICIT NONE

      INTEGER, PRIVATE  :: nsmin, nsmax
      INTEGER, PRIVATE, PARAMETER :: isym=0, iasym=1, m0=0, m1=1, m2=2, &
                                     n0 = 0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE ::               &
         bsupsijh, bsupuijh, bsupvijh, pijh, KxBsij, KxBuij, KxBvij


      CONTAINS
!>  \brief Parallel/Serial update MHD force on full radial mesh
      SUBROUTINE update_force
!     
!     PURPOSE: UPDATES FORCE using the advanced values of B and p obtained
!              from calls to UPDATE_PRES and UPDATE_BFIELD
!              Updates values of FSUB?MNF (on the FULL radial mesh) and stores in QUANTITIES MODULE
!
!      Linearized force is Flin ~ delta_v has ONLY the linear part (not DC part and not non-linear)
!
       USE siesta_currents, ONLY: cv_currents
       USE bhtobf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, n1, n2, i, j, irefl, jrefl
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
         pmncf_ds, pmncf, pmnch,                                        &
         pmnsf_ds, pmnsf, pmnsh,                                        &
         ksupsijf, ksupuijf, ksupvijf,                                  &
         bsupsijf1,bsupuijf1,bsupvijf1, pijf1,                          &
         KxBsmnc, KxBumns, KxBvmns,                                     &
         KxBsmns, KxBumnc, KxBvmnc
       REAL(dp)    :: ton, toff, skston, skstoff, wp0
       LOGICAL     :: ltest=.TRUE.
!-----------------------------------------------
       CALL second0(skston)
       nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow+1)
       n1 = nsmin;  n2 = nsmax

       ALLOCATE(bsupsijh(ntheta,nzeta,n1:n2),                           &
                bsupuijh(ntheta,nzeta,n1:n2),                           &
                bsupvijh(ntheta,nzeta,n1:n2),                           &
                pijh(ntheta,nzeta,n1:n2), stat=istat)  
       CALL ASSERT(istat.EQ.0,'Allocation1 failed in UPDATE_FORCE')


!Convert to real-space contravariant components of B and p
       CALL IncFields ( jbsupsmnsh,  jbsupumnch,  jbsupvmnch,  jpmnch,  &
                       djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch,  &
                        isym)
       IF (lasym)                                                       &
       CALL IncFields ( jbsupsmnch,  jbsupumnsh,  jbsupvmnsh,  jpmnsh,  &
                       djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh,  &
                        iasym)

!UPDATE THERMAL (PRESSURE) ENERGY (pijh = jac * pres at this point)
       WPRES: IF (l_getwmhd) THEN
          CALL ASSERT(.NOT.INHESSIAN,                                   &
             'l_getwmhd must be set to FALSE in Hessian')
          n2=MIN(ns,endglobrow)
          IF (nsmin .EQ. 1) pijh(:,:,1) = 0
          wp = SUM(pijh(:,:,n1:n2)*wint(:,:,n1:n2))
          wp0 = signjac*(twopi*twopi * hs_i)*wp
          n2=nsmax
          IF (PARSOLVER) THEN
#if defined(MPI_OPT)
          CALL MPI_ALLREDUCE(wp0,wp,1,MPI_REAL8, MPI_SUM,               &
                             SIESTA_COMM,MPI_ERR)
#endif
          ELSE
          wp = wp0
          END IF
       END IF WPRES

       ALLOCATE(ksupsijf(ntheta,nzeta,n1:n2),                           & ! Full mesh quantities (real space)
                ksupuijf(ntheta,nzeta,n1:n2),                           & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijf(ntheta,nzeta,n1:n2),                           &
                bsupsijf1(ntheta,nzeta,n1:n2),                          &
                bsupuijf1(ntheta,nzeta,n1:n2),                          &
                bsupvijf1(ntheta,nzeta,n1:n2),                          &
                pijf1(ntheta,nzeta,n1:n2), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation2 failed in UPDATE_FORCE')

!
!Calculate updated FULL-MESH components of BSUPX and PRESSURE in real-space by averaging
!
       bsupsijh = bsupsijh/jacobh(:,:,n1:n2)
       bsupvijh = bsupvijh/jacobh(:,:,n1:n2)
       pijh = pijh/jacobh(:,:,n1:n2)
       CALL bhalftobfull (bsupsijh, bsupuijh, bsupvijh,                 &
                          bsupsijf1,bsupuijf1,bsupvijf1)
                                 
!
!Calculate contravariant components of J = curl B (Ksup = gsqrt*Jsup)
!and update the magnetic energy if nonlinear force is being computed
!Note: origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents (bsupsijh, bsupuijh, bsupvijh,                  &
                         ksupsijf, ksupuijf, ksupvijf,                  &
                         ( .NOT.l_linearize .OR. l_getwmhd), .FALSE.)

       DEALLOCATE(bsupsijh, bsupuijh, bsupvijh, pijf1, stat=istat)

       nsmax=MIN(ns,endglobrow); n2 = nsmax

       ALLOCATE(KxBsij(ntheta,nzeta,n1:n2),                             &
                KxBuij(ntheta,nzeta,n1:n2),                             &
                KxBvij(ntheta,nzeta,n1:n2), stat=istat)  
       CALL ASSERT(istat.EQ.0,'Allocation3 failed in UPDATE_FORCE')

       IF (l_linearize) THEN
!dJ X B0 (ksupXijf is the perturbed current)
          KxBsij = bsupvijf0(:,:,n1:n2)*ksupuijf(:,:,n1:n2)             &
                 - bsupuijf0(:,:,n1:n2)*ksupvijf(:,:,n1:n2)
          KxBuij = bsupsijf0(:,:,n1:n2)*ksupvijf(:,:,n1:n2)             &
                 - bsupvijf0(:,:,n1:n2)*ksupsijf(:,:,n1:n2)
          KxBvij = bsupuijf0(:,:,n1:n2)*ksupsijf(:,:,n1:n2)             &
                 - bsupsijf0(:,:,n1:n2)*ksupuijf(:,:,n1:n2)

!J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!bsupXijf1 is the perturbed B-field

          KxBsij = KxBsij                                               &
          +  (bsupvijf1(:,:,n1:n2)*ksupuijf0(:,:,n1:n2)                 &
          -   bsupuijf1(:,:,n1:n2)*ksupvijf0(:,:,n1:n2))
          KxBuij = KxBuij                                               &
          +  (bsupsijf1(:,:,n1:n2)*ksupvijf0(:,:,n1:n2)                 &
          -   bsupvijf1(:,:,n1:n2)*ksupsijf0(:,:,n1:n2))
          KxBvij = KxBvij                                               &
          +  (bsupuijf1(:,:,n1:n2)*ksupsijf0(:,:,n1:n2)                 &
          -   bsupsijf1(:,:,n1:n2)*ksupuijf0(:,:,n1:n2))
       ELSE
          KxBsij = bsupvijf1(:,:,n1:n2)*ksupuijf(:,:,n1:n2)             &
                 - bsupuijf1(:,:,n1:n2)*ksupvijf(:,:,n1:n2)
          KxBuij = bsupsijf1(:,:,n1:n2)*ksupvijf(:,:,n1:n2)             &
                 - bsupvijf1(:,:,n1:n2)*ksupsijf(:,:,n1:n2)
          KxBvij = bsupuijf1(:,:,n1:n2)*ksupsijf(:,:,n1:n2)             &
                 - bsupsijf1(:,:,n1:n2)*ksupuijf(:,:,n1:n2)
       END IF

       IF (nsmax .EQ. ns) ltest = ALL(bsupsijf1(:,:,ns).EQ.zero)
       CALL ASSERT(ltest, 'bsupsijf1(ns) != 0 in UPDATE_FORCE')

       CALL GetMHDForce(fsubsmncf, fsubumnsf, fsubvmnsf, isym)
       IF (lasym)                                                       &     
	   CALL GetMHDForce(fsubsmnsf, fsubumncf, fsubvmncf, iasym)

       DEALLOCATE(ksupsijf, ksupuijf, ksupvijf,                         &
	              bsupsijf1,bsupuijf1,bsupvijf1,                        &
                  KxBsij, KxBuij, KxBvij, pijh, stat=istat)

       CALL second0(skstoff)
       time_update_force = time_update_force + (skstoff-skston)

      END SUBROUTINE update_force


      SUBROUTINE IncFields( jbsupsmnh,  jbsupumnh,  jbsupvmnh,  jpmnh,  &
                           djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh,  &
                            ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN) ::            & !1:ns
          jbsupsmnh,  jbsupumnh,  jbsupvmnh,  jpmnh
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN) ::            & !nsmin:nsmax
         djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, fcomb, istat
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         bsupsmnh_i, bsupumnh_i, bsupvmnh_i, pmnh_i
      REAL(dp)           :: xfact
!----------------------------------------------- 
      IF (ipar .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
      END IF
!
!     Allocate space for TOTAL values of BSUPXMNH_i and pmnh_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments were calculated previously in UPDATE_BFIELD and UPDATE_PRES
!
      ALLOCATE(bsupsmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),               &
               bsupumnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),               &
               bsupvmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),               &  
               pmnh_i    (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
      CALL ASSERT(istat.EQ.0,'Allocation failed in INCFIELDS')
!
!     IF nonlinear force, add perturbation to DC (initial) value
!     IF linear force, compute JUST perturbed fields here (add unperturbed force gc0 in funct_island)
!
      xfact = 1
      IF (l_linearize) xfact = 0
      bsupsmnh_i = xfact*jbsupsmnh(:,:,nsmin:nsmax) + djbsupsmnh(:,:,nsmin:nsmax)
      bsupumnh_i = xfact*jbsupumnh(:,:,nsmin:nsmax) + djbsupumnh(:,:,nsmin:nsmax)
      bsupvmnh_i = xfact*jbsupvmnh(:,:,nsmin:nsmax) + djbsupvmnh(:,:,nsmin:nsmax)
      pmnh_i     = xfact*jpmnh(:,:,nsmin:nsmax)     + djpmnh(:,:,nsmin:nsmax)

! Calculate real-space BSUPX on half mesh
      CALL toijsp_par(bsupsmnh_i, bsupsijh, fcomb, fours)   
      CALL toijsp_par(bsupumnh_i, bsupuijh, fcomb, fouru)
      CALL toijsp_par(bsupvmnh_i, bsupvijh, fcomb, fourv)
      CALL toijsp_par(pmnh_i, pijh, fcomb, fourv)
 
      DEALLOCATE (bsupsmnh_i, bsupumnh_i, bsupvmnh_i, pmnh_i)
      
      END SUBROUTINE IncFields

!*******************************************************************************
!>  \brief Compute covariant (sub) components of the MHD force in parallel.
!>
!*******************************************************************************
     
      SUBROUTINE GetMHDForce(fsubsmnf, fsubumnf, fsubvmnf, ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), POINTER, DIMENSION(:,:,:) :: fsubsmnf, fsubumnf, fsubvmnf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  ::  fours, fouru, fourv, istat, k1, k2, k3, k4
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)                           &
               ::  KxBsmn, KxBumn, KxBvmn, pmnh, pmnf, pmnf_ds
      REAL(dp) :: p00
!----------------------------------------------- 
      IF (ipar .EQ. isym) THEN
         fours = f_cos; fouru = f_sin; fourv = f_sin
      ELSE
         fours = f_sin; fouru = f_cos; fourv = f_cos
      END IF

      k1 = nsmin; k2 = nsmax; k3 = k1; k4 = MIN(ns,endglobrow+1)

      ALLOCATE(KxBsmn(0:mpol,-ntor:ntor,k1:k2),                         &
               KxBumn(0:mpol,-ntor:ntor,k1:k2),                         &
               KxBvmn(0:mpol,-ntor:ntor,k1:k2),                         &
               pmnf_ds(0:mpol,-ntor:ntor,k3:k4),                        &
               pmnf(0:mpol,-ntor:ntor,k3:k4),                           &
               pmnh(0:mpol,-ntor:ntor,k3:k4), stat=istat)         
      CALL ASSERT(istat.EQ.0,'Allocation failed in GETMHDFORCE')

!Harmonics of Lorentz force on full grid
      CALL tomnsp_par(KxBsij, KxBsmn, fours)
      CALL tomnsp_par(KxBuij, KxBumn, fouru)
      CALL tomnsp_par(KxBvij, KxBvmn, fourv)

!Harmonics of "true" pressure (after dividing out jacobian factor in pijh)
      CALL tomnsp_par(pijh, pmnh, fours)

      CALL to_full_mesh(pmnh, pmnf)
      IF (nsmax .EQ. ns) THEN
         p00 = pmnf(m0, n0, ns)
         pmnf(:,:,ns) = 0; pmnf(m0, n0, ns) = p00
      END IF

!Radial derivative of pressure valid for [k3, k4-1]
      CALL GradientFull(pmnf_ds, pmnh)

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
      IF (k3 .EQ. 1) THEN
         pmnf(:,:,1) = 0; pmnf_ds(:,:,1) = 0
         pmnf(m0:m1,:,1) = pmnh(m0:m1,:,2)
         pmnf_ds(m1,:,1) = 2*ohs*pmnh(m1,:,2)
      END IF

!
!Calculate harmonics of MHD force (F = JxB - grad p) on full mesh
!For l_linearize=T:  PERTURBED Force (ignore F(xc=0) part)
!For l_linearize=F:  NONLINEAR Force
!
      CALL ASSERT(nsmax.eq.k2,'NSMAX != K2')
      KxBsmn = KxBsmn - pmnf_ds(:,:,k1:k2)
      CALL get_force_harmonics(pmnf, KxBsmn, KxBumn, KxBvmn,            &
                               fsubsmnf, fsubumnf, fsubvmnf, ipar)
      DEALLOCATE(KxBsmn, KxBumn, KxBvmn, pmnf, pmnf_ds, pmnh, stat=istat)

      END SUBROUTINE GetMHDForce

!*******************************************************************************
!>  \brief Compute harmonics of the MHD covariant force components
!>
!>  The MHD forces are pressure gradient and the magnetic pressure
!>
!>     F = J X B - Grad(p)                                                 (1)
!>
!>  where J is the plasma current density. In curvilinear coodinates, the forces
!>  components becomes
!>
!>     F_s = (K^u*B^v - K^v*B^u) - dpds                                    (2a)
!>
!>     F_u = (K^v*B^s - K^s*B^v) - dpdu                                    (2b)
!>
!>     F_v = (K^s*B^u - K^u*B^s) - dpdv                                    (2c)
!>
!>  where K^i = sqrt(g)*J^i . Note this routine only works on a
!>  single parity at a time.
!>
!>  \param[in]  pmnf        Fourier pressure.
!>  \param[in]  kxb_smn     Fourier magnetic pressure and radial pressure
!>                          gradient in the radial direction.
!>  \param[in]  kxb_umn     Fourier magnetic pressure in the poloidal direction.
!>  \param[in]  kxb_vmn     Fourier magnetic pressure in the toroidal direction.
!>  \param[out] f_smnf      F_s component of the MHD force. (cos, ipar=0)
!>  \param[out] f_umnf      F_u component of the MHD force. (sin, ipar=0)
!>  \param[out] f_vmnf      F_v component of the MHD force. (sin, ipar=0)
!>  \param[in]  ipar        =isym (default lasym=F), =iasym (lasym=T)
!*******************************************************************************
       SUBROUTINE get_force_harmonics(pmnf,                             &
                                      KxB_smn, KxB_umn, KxB_vmn,        &
                                      f_smnf, f_umnf, f_vmnf, ipar) 
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE shared_data, ONLY: nprecon, l_PrintOriginForces,             &
                              l_push_s, l_push_u, l_push_v, l_push_edge,&
                              col_scale
       USE nscalingtools, ONLY: FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN)  ::          &
           KxB_smn, KxB_umn, KxB_vmn, pmnf
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT) ::        &
           f_smnf, f_umnf, f_vmnf
       INTEGER, INTENT(IN) ::  ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp), PARAMETER :: p5 = 0.5_dp
       INTEGER             :: m, n, mp, np, mcol, ioff, joff
       REAL(dp)            :: skston, skstoff 
       REAL(dp)            :: f1, f2, f3, r0, tmps(6) 
       INTEGER             :: from, sparity
       LOGICAL             :: ltest
!-----------------------------------------------
!
!      NOTE: there is an m=1 part of pmnf retained at js=1, which is REALLY at js=hs/2
!      It contributes a spurious piece to f_vmnf, which is correctly zeroed
!
       CALL second0(skston)
       IF (ipar .EQ. isym) THEN
          mcol=1; sparity = 1
          r0 = hs_i/2
       ELSE
          mcol=4; sparity = -1
          r0 =-hs_i/2
       END IF
       
       nsmin = LBOUND(KXB_smn,3);  nsmax = UBOUND(KXB_smn,3)

       f_smnf(:,:,nsmin:nsmax) = KxB_smn
       DO m = 0, mpol
          mp = m*sparity
          DO n = -ntor, ntor       
             np = n*nfp*sparity
             f_umnf(m,n,nsmin:nsmax) = KxB_umn(m,n,nsmin:nsmax)         &
                               + mp*pmnf(m,n,nsmin:nsmax)                 ! mp ~ -(d/du)pres
             f_vmnf(m,n,nsmin:nsmax) = KxB_vmn(m,n,nsmin:nsmax)         &
                               + np*pmnf(m,n,nsmin:nsmax)                 ! np ~ -(d/dv)pres
          ENDDO
       ENDDO

       ltest = ALL(f_smnf(m0,-ntor:-1,nsmin:nsmax).EQ.zero)
       CALL ASSERT(ltest,'f_smnf != 0, m=0, n<0')
       ltest = ALL(f_umnf(m0,-ntor:-1,nsmin:nsmax).EQ.zero)
       CALL ASSERT(ltest,'f_umnf != 0, m=0, n<0')
       ltest = ALL(f_vmnf(m0,-ntor:-1,nsmin:nsmax).EQ.zero)
       CALL ASSERT(ltest,'f_vmnf != 0, m=0, n<0')

!      ADD SYMMETRIC DAMPING TERMS -B v dot B TO SUPPRESS NULL SPACE
!      WORKS BEST: comment out IF statement below, allowing -mu v|| term
!      to be added to BOTH linear and nonlinear forces. However, must
!      guarantee that v->0 at end for this to work well

!CALL funct_island WITH vsupXmn...=0 to get un-augmented MHD force
       IF (muPar.NE.zero .AND. l_Compute_Hessian)                       &
!       IF (muPar .NE. zero)                                             &
       CALL AddParDamping(f_smnf, f_umnf, f_vmnf, KxB_vmn, ipar)

!GATHER BDY FORCES AND PRINT OUT ON PROC=0
       PRINT_O: IF (l_PrintOriginForces) THEN
             fbdy(1) = SQRT(SUM(f_smnf(m1,:,1)**2))
             fbdy(2) = SQRT(SUM(f_smnf(m1,:,2)**2))
             fbdy(3) = SQRT(SUM(f_smnf(m0,:,2)**2) + SUM(f_smnf(m2:,:,2)**2))

             fbdy(4) = SQRT(SUM(f_umnf(m1,:,1)**2))
             fbdy(5) = SQRT(SUM(f_umnf(m1,:,2)**2))
             fbdy(6) = SQRT(SUM(f_umnf(m0,:,2)**2) + SUM(f_umnf(m2:,:,2)**2))

             fbdy(7) = SQRT(SUM(f_vmnf(m0,:,1)**2))
             fbdy(8) = SQRT(SUM(f_vmnf(m0,:,2)**2))
             fbdy(9) = SQRT(SUM(f_vmnf(m1:,:,2)**2))

             fbdy(11) = SQRT(SUM(f_umnf(:,:,ns)**2))
             fbdy(10) = SQRT(SUM(f_umnf(:,:,ns-1)**2))
             fbdy(13) = SQRT(SUM(f_vmnf(:,:,ns)**2))
             fbdy(12) = SQRT(SUM(f_vmnf(:,:,ns-1)**2))

#if defined(MPI_OPT)
           NPROCS_0: IF (nprocs.GT.1) THEN
             IF (.NOT.FIRSTPARTITION_GE_2) THEN
                tmps(1)=fbdy(2); tmps(2)=fbdy(3)
                tmps(3)=fbdy(5); tmps(4)=fbdy(6)
                tmps(5)=fbdy(8); tmps(6)=fbdy(9)
                from=2
                CALL MPI_BCast(tmps,6,MPI_REAL8,from-1,SIESTA_COMM,MPI_ERR)
                fbdy(2)=tmps(1); fbdy(3)=tmps(2)
                fbdy(5)=tmps(3); fbdy(6)=tmps(4)
                fbdy(8)=tmps(5); fbdy(9)=tmps(6)
             END IF

             IF (.NOT.LASTPARTITION_GE_2) THEN
                tmps(1)=fbdy(11); tmps(2)=fbdy(13)
                from=nprocs-1
                CALL MPI_BCast(tmps,2,MPI_REAL8,from-1,SIESTA_COMM,MPI_ERR)
                fbdy(11)=tmps(1); fbdy(13)=tmps(2)
             END IF

             tmps(1)=fbdy(10); tmps(2)=fbdy(11)
             tmps(3)=fbdy(12); tmps(4)=fbdy(13)
             from=nprocs
             CALL MPI_BCast(tmps,4,MPI_REAL8,from-1,SIESTA_COMM,MPI_ERR)
             fbdy(10)=tmps(1); fbdy(11)=tmps(2)
             fbdy(12)=tmps(3); fbdy(13)=tmps(4)
           END IF NPROCS_0
#endif
           IF (iam .EQ. 0 .AND. lverbose) THEN
              IF (ipar .EQ. 0) PRINT 120,'SYMMETRIC FORCES: '
              IF (ipar .EQ. 1) PRINT *,'ASYMMETRIC FORCES: '

              PRINT 121, fbdy(1),fbdy(2),fbdy(3)
              PRINT 122, fbdy(4),fbdy(5),fbdy(6)
              PRINT 123, fbdy(7),fbdy(8),fbdy(9)
              PRINT 124, fbdy(10),fbdy(11)
              PRINT 125, fbdy(12),fbdy(13)
           END IF

       END IF PRINT_O

 120   FORMAT(/,1x,a)
 121   FORMAT(' fs(1,m=1): ',1p,e12.3,' fs(2,m=1): ',1pe12.3,' fs(2,m!=1):',1pe12.3)
 122   FORMAT(' fu(1,m=1): ',1p,e12.3,' fu(2,m=1): ',1pe12.3,' fu(2,m!=1):',1pe12.3)
 123   FORMAT(' fv(1,m=0): ',1p,e12.3,' fv(2,m=0): ',1pe12.3,' fv(2,m>0): ',1pe12.3)
 124   FORMAT(' fu(ns-1) : ',1p,e12.3,' fu(ns)   : ',1pe12.3)
 125   FORMAT(' fv(ns-1) : ',1p,e12.3,' fv(ns)   : ',1pe12.3,/)

!      ORIGIN BOUNDARY CONDITIONS
       NSMIN_1: IF (nsmin .EQ. 1) THEN
       f_smnf(m0,:,1) = 0; f_smnf(m2:,:,1) = 0
!       IF (.NOT.l_natural) f_umnf(m1,:,1) = f_umnf(m1,:,1)-r0*f_smnf(m1,:,1) 
       IF (.NOT.l_push_s) THEN 
          f_smnf(m1,:,1) = 0
       ELSE
          f_smnf(m1,:,1) = p5*f_smnf(m1,:,1)
       END IF

       f_umnf(m0,:,1) = 0; f_umnf(m2:,:,1) = 0
       IF (.NOT.l_push_u) THEN
          f_umnf(m1,:,1) = 0
       ELSE
          f_umnf(m1,:,1) = p5*f_umnf(m1,:,1)
       END IF

       f_vmnf(m1:,:,1) = 0
       IF (.NOT.l_push_v) THEN
          f_vmnf(m0,:,1) = 0
       ELSE
          f_vmnf(m0,:,1) = p5*f_vmnf(m0,:,1)                              !MAINTAIN SYMMETRY OF HESSIAN
       END IF

       END IF NSMIN_1
 
!      EDGE BOUNDARY CONDITIONS
       NSMAX_NS: IF (nsmax .EQ. ns) THEN
          f_smnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             f_umnf(:,:,ns) = 0
             f_vmnf(:,:,ns) = 0
          ELSE                                                            !MAINTAIN SYMMETRY OF HESSIAN
             f_umnf(:,:,ns) = p5*f_umnf(:,:,ns)
             f_vmnf(:,:,ns) = p5*f_vmnf(:,:,ns)
          END IF
       END IF NSMAX_NS

       f_smnf(:,:,nsmin:nsmax) = signjac*f_smnf(:,:,nsmin:nsmax)        &
                               * col_scale(:,:,mcol,nsmin:nsmax)
       f_umnf(:,:,nsmin:nsmax) = signjac*f_umnf(:,:,nsmin:nsmax)        &
                               * col_scale(:,:,mcol+1,nsmin:nsmax)
       f_vmnf(:,:,nsmin:nsmax) = signjac*f_vmnf(:,:,nsmin:nsmax)        &
                               * col_scale(:,:,mcol+2,nsmin:nsmax)

       CALL second0(skstoff)
       get_force_harmonics_time=get_force_harmonics_time+(skstoff-skston)

       END SUBROUTINE get_force_harmonics


       SUBROUTINE addpardamping (f_smnf, f_umnf, f_vmnf, work, ipar)
       USE hessian, ONLY: muPar, mupar_norm, l_Compute_Hessian
       USE shared_data, ONLY:  col_scale
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
       REAL(dp)            :: work(0:mpol,-ntor:ntor,nsmin:nsmax)
       REAL(dp), INTENT(INOUT), DIMENSION(0:mpol,-ntor:ntor,ns) ::      &
                 f_smnf, f_umnf, f_vmnf          
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: js, istat, mcol, fours, fouru, fourv
       REAL(dp) :: av(ns)
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: vdotB, tmp
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: xsubsij, xsubuij, xsubvij
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!
       IF (ipar .EQ. isym) THEN
          mcol = 1; fours=f_cos; fouru=f_sin; fourv=f_sin
       ELSE
          mcol = 4; fours=f_sin; fouru=f_cos; fourv=f_cos
       END IF

       ALLOCATE (vdotB(ntheta,nzeta,nsmin:nsmax),                       &
                 tmp(ntheta,nzeta,nsmin:nsmax), stat=istat)

!SPH 05-6-16 UPDATE mupar_norm IN cv_currents and put in 1/jacobh factor
!SIMPLER (FASTER, APPROXIMATE) FORM
!
!       GOTO 200

       tmp = bsubsijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupsmncf(:,:,js)
          ELSE
             work(:,:,js) = jvsupsmnsf(:,:,js)
          END IF
          f_smnf(:,:,js) = f_smnf(:,:,js)                               &
                         + work(:,:,js)*av(js)*col_scale(:,:,mcol,js)
       END DO   

       tmp = bsubuijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupumnsf(:,:,js)
          ELSE
             work(:,:,js) = jvsupumncf(:,:,js)
          END IF
          f_umnf(:,:,js) = f_umnf(:,:,js)                               &
                         + work(:,:,js)*av(js)*col_scale(:,:,mcol+1,js) 
       END DO   

       tmp = bsubvijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
       CALL SurfAverage(av(nsmin:nsmax), tmp, nsmin, nsmax)

       DO js = nsmin,nsmax
          av(js) = av(js)*muPar*mupar_norm(js)
          IF (ipar .EQ. 0) THEN
             work(:,:,js) = jvsupvmnsf(:,:,js)
          ELSE
             work(:,:,js) = jvsupvmncf(:,:,js)
          END IF
          f_vmnf(:,:,js) = f_vmnf(:,:,js)                               &
                         + work(:,:,js)*av(js)*col_scale(:,:,mcol+2,js)
       END DO   

       GO TO 100       

 200   CONTINUE

!
!col_scale factor already in real-space jvsupXijf
       IF (nsmin .EQ. 1) bsubuijf(:,:,1) = 0
       vdotB(:,:,nsmin:nsmax) =                                         &
           jvsupsijf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)         &
         + jvsupuijf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)         &
         + jvsupvijf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

       DO js = nsmin,nsmax
          vdotB(:,:,js) = muPar*mupar_norm(js)*vdotB(:,:,js)/jacobf(:,:,js)
       END DO

       tmp = bsubsijf(:,:,nsmin:nsmax)*vdotB
       CALL tomnsp_par(tmp, work, fours)
       DO js=nsmin,nsmax
          f_smnf(:,:,js) = f_smnf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubuijf(:,:,nsmin:nsmax)*vdotB
       CALL tomnsp_par(tmp, work, fouru)
       DO js=nsmin,nsmax
          f_umnf(:,:,js) = f_umnf(:,:,js) + work(:,:,js)
       END DO

       tmp = bsubvijf(:,:,nsmin:nsmax)*vdotB
       CALL tomnsp_par(tmp, work, fourv)
       DO js=nsmin,nsmax
          f_vmnf(:,:,js) = f_vmnf(:,:,js) + work(:,:,js)
       END DO

 100   CONTINUE
       DEALLOCATE (vdotB, tmp)

       END SUBROUTINE addpardamping

       END MODULE siesta_force
