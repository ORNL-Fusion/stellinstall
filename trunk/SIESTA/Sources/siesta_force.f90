!*******************************************************************************
!>  \file siesta_force.f90
!>  \brief Contains the \ref siesta_force module.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  \brief Compute the JxB - Grad(p) covariant force components
!>  \author S. P. Hirshman and R. Sanchez
!>  \date August 21, 2006
!>
!>  The MHD is in equilibrium when the force of the magnetic pressure balances
!>  the plasma pressure.
!>
!>     F = J X B - Grad(p) = 0                                               (1)
!>
!>  where J is the current density. Updates for B (dB) and p (dp) were computed in
!>  \ref siesta_bfield and \ref siesta_pressure, respectively.
!>
!>     B_i+1 = B_i + dB                                                     (2a)
!>
!>     p_i+1 = p_i + dp                                                     (2b)
!>
!>  Note that the magnetic field and pressure terms contain a Jacobian factor. 
!>  Updated current density can be found from the updated magnetic field.
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
      USE v3_utilities, ONLY: assert
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE fourier 
      USE island_params, ONLY: nu_i, hs_i, dnorm_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd, lverbose, l_pedge
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE quantities
      USE mpi_inc

      IMPLICIT NONE

      PRIVATE
      
      INTEGER, PARAMETER :: isym=0, iasym=1, m0=0, m1=1, m2=2, n0=0
      INTEGER  :: nsmin, nsmax
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
         bsupsijh, bsupuijh, bsupvijh, pijh, KxBsij, KxBuij, KxBvij
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
         jbsupsmnh_i, jbsupumnh_i, jbsupvmnh_i, jpmnh_i
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
         pmnh, pmnf, pmnf_ds
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)    ::  pardamp
      
      PUBLIC :: update_force

      CONTAINS
!>  \brief Parallel/Serial update MHD forces on full radial mesh
      SUBROUTINE update_force
!     
!     PURPOSE: UPDATES FORCE using the advanced values of B and p obtained
!              from calls to UPDATE_PRES and UPDATE_BFIELD
!              Updates values of FSUB?MNF (on the FULL radial mesh) and stores in QUANTITIES MODULE
!
!      Linearized force is Flin ~ delta_v has ONLY the linear part (not DC part and not non-linear)
!
       USE shared_data, ONLY: gc, col_scale
       USE siesta_currents, ONLY: cv_currents
       USE hessian, ONLY: Apply_ColScale
       USE bhtobf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
                 pijf1
       INTEGER  :: istat, n1, n2, i, j, k
       REAL(dp)    :: ton, toff, skston, skstoff
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


!Compute real-space contravariant components of B and p
!
!     Allocate space for TOTAL values of BSUPXMNH_i and pmnh_i
!     (fields are updated in UPDATE_STATE, not here => UNPERTURBED fields 
!      are unchanged until UPDATE_STATE call)
!     perturbed fields were calculated in UPDATE_BFIELD and UPDATE_PRES
!
       ALLOCATE(jbsupsmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                jbsupumnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                jbsupvmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &  
                jpmnh_i    (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       CALL ASSERT(istat.EQ.0,'Allocation failed before INCFIELDS')
       CALL IncFields ( jbsupsmnsh,  jbsupumnch,  jbsupvmnch,  jpmnch,  &
                       djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch,  &
                        isym)
       IF (lasym)                                                       &
       CALL IncFields ( jbsupsmnch,  jbsupumnsh,  jbsupvmnsh,  jpmnsh,  &
                       djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh,  &
                        iasym)
       DEALLOCATE(jbsupsmnh_i,jbsupumnh_i,jbsupvmnh_i,jpmnh_i,stat=istat)

!UPDATE THERMAL (PRESSURE) ENERGY (pijh = jac * pres at this point)
       WPRES: IF (l_getwmhd) THEN
          CALL ASSERT(.NOT.INHESSIAN,                                   &
             'l_getwmhd must be set to FALSE in Hessian')
          n2=MIN(ns,endglobrow)
          IF (nsmin .EQ. 1) n1 = MAX(2,n1)
          wp = SUM(pijh(:,:,n1:n2)*wint(:,:,n1:n2))
          wp0 = signjac*(twopi*twopi * hs_i)*wp
          n1 = nsmin
          n2=nsmax
          IF (PARSOLVER) THEN
#if defined(MPI_OPT)
             CALL MPI_ALLREDUCE(wp0, wp, 1, MPI_REAL8, MPI_SUM,         & 
                                SIESTA_COMM, MPI_ERR)
#endif
          ELSE
             wp = wp0
          END IF
       END IF WPRES

       ALLOCATE(pijf1(ntheta,nzeta,n1:n2), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation2 failed in UPDATE_FORCE')
!
!Update FULL-MESH components of BSUPX and PRESSURE in real-space by averaging
       CALL bhalftobfull (bsupsijh, bsupuijh, bsupvijh,                 &
                          bsupsijf, bsupuijf, bsupvijf, pijh, pijf1)
                                 
!
!Calculate contravariant components of J = curl B (Ksup = gsqrt*Jsup)
!and update the magnetic energy if nonlinear force is being computed
       CALL cv_currents (bsupsijh, bsupuijh, bsupvijh,                  &
                         ksupsijf, ksupuijf, ksupvijf,                  &
                         (.NOT.l_linearize .OR. l_getwmhd), .FALSE.)

       DEALLOCATE(bsupsijh, bsupuijh, bsupvijh, stat=istat)

       nsmax=MIN(ns,endglobrow); n2 = nsmax

       ALLOCATE(KxBsij(ntheta,nzeta,n1:n2),                             &
                KxBuij(ntheta,nzeta,n1:n2),                             &
                KxBvij(ntheta,nzeta,n1:n2), stat=istat)  
       CALL ASSERT(istat.EQ.0,'Allocation3 failed in UPDATE_FORCE')

       CALL InitParDamping (KxBsij)                                       !KxBsij -> scratch array

!Compute curvilinear Lorentz force components
       IF (l_linearize) THEN
!dJ X B0
          CALL Lorentz(bsupsijf0, bsupuijf0, bsupvijf0,                 &
                       ksupsijf,  ksupuijf,  ksupvijf, f_none)
!J0 X dB
          CALL Lorentz(bsupsijf,  bsupuijf,  bsupvijf,                  &
                       ksupsijf0, ksupuijf0, ksupvijf0, f_sum)
       ELSE
!(J0+dJ) X (B0+dB)
          CALL Lorentz(bsupsijf, bsupuijf, bsupvijf,                    &
                       ksupsijf, ksupuijf, ksupvijf,  f_none)
       END IF

       IF (nsmax.EQ.ns .AND. .NOT.l_vessel) THEN
          ltest = ALL(bsupsijf(:,:,ns).EQ.zero)
          CALL ASSERT(ltest, 'bsupsijf(ns) != 0 in UPDATE_FORCE')
       END IF

       n2 = MIN(ns,endglobrow+1)

       ALLOCATE(pmnf_ds(0:mpol,-ntor:ntor,nsmin:n2),                    &
                pmnf(0:mpol,-ntor:ntor,nsmin:n2),                       &
                pmnh(0:mpol,-ntor:ntor,nsmin:n2), stat=istat)         
       CALL ASSERT(istat.EQ.0,'Allocation failed before GETMHDFORCE')

       CALL GetMHDForce(fsubsmncf, fsubumnsf, fsubvmnsf, isym)
       IF (lasym)                                                       &
	   CALL GetMHDForce(fsubsmnsf, fsubumncf, fsubvmncf, iasym)

       CALL Apply_ColScale (gc, col_scale, nsmin, nsmax)                  
	   
       DEALLOCATE(pmnf, pmnf_ds, pmnh, pijf1,                           &
                  KxBsij, KxBuij, KxBvij, pijh, stat=istat)
       IF (ALLOCATED(pardamp)) DEALLOCATE(pardamp)

       CALL second0(skstoff)
       time_update_force = time_update_force + (skstoff-skston)

      END SUBROUTINE update_force

!>  \brief Compute nonlinear or linearized contravariant magnetic field
!>   components and pressure in real space
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
      INTEGER  :: fours, fouru, fourv, fcomb, istat, i, j, k
!----------------------------------------------- 
      IF (ipar .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
      END IF
!
!     IF nonlinear force, add perturbation to unperturbed value
!     IF linear force, compute JUST perturbed fields here
!
      
      jbsupsmnh_i = djbsupsmnh
      jbsupumnh_i = djbsupumnh
      jbsupvmnh_i = djbsupvmnh
      jpmnh_i     = djpmnh
      IF (.NOT.l_linearize) THEN
         jbsupsmnh_i = jbsupsmnh_i + jbsupsmnh(:,:,nsmin:nsmax)
         jbsupumnh_i = jbsupumnh_i + jbsupumnh(:,:,nsmin:nsmax)
         jbsupvmnh_i = jbsupvmnh_i + jbsupvmnh(:,:,nsmin:nsmax)
         jpmnh_i     = jpmnh_i     + jpmnh(:,:,nsmin:nsmax)
      END IF

! js=1 origin values
      IF (nsmin .EQ. 1) THEN
         jbsupsmnh_i(:,:,1) = 0; jbsupsmnh_i(m1,:,1) = jbsupsmnh_i(m1,:,2)
         jbsupumnh_i(:,:,1) = 0; jbsupumnh_i(m0:m1,:,1) = jbsupumnh_i(m0:m1,:,2)
         jbsupvmnh_i(:,:,1) = 0; jbsupvmnh_i(m0,:,1) = jbsupvmnh_i(m0,:,2)
         jpmnh_i(:,:,1) = 0;     jpmnh_i(m0,:,1) = jpmnh_i(m0,:,2)
      END IF

! Calculate real-space JBSUPXIJH and JPIJH (half mesh)
      CALL toijsp_par(jbsupsmnh_i, bsupsijh, fcomb, fours)   
      CALL toijsp_par(jbsupumnh_i, bsupuijh, fcomb, fouru)
      CALL toijsp_par(jbsupvmnh_i, bsupvijh, fcomb, fourv)
      CALL toijsp_par(jpmnh_i, pijh, fcomb, fourv)
       
      END SUBROUTINE IncFields

!>  \brief Compute covariant (sub) eal-space components of the Lorentz K X B force
      SUBROUTINE Lorentz (bsupsijf, bsupuijf, bsupvijf,                 &
                          ksupsijf, ksupuijf, ksupvijf, fcomb)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: fcomb
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN) ::            &
             bsupsijf, bsupuijf, bsupvijf, ksupsijf, ksupuijf, ksupvijf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  ::  n1, n2
!-----------------------------------------------
      n1 = LBOUND(KxBsij,3);  n2 = UBOUND(KXBsij,3)

      
 !Initialize to avoid NaN errors!
      IF (fcomb .EQ. f_none) THEN
         KxBsij = 0; KxBuij = 0; KxBvij = 0
      END IF

      KxBsij = KxBsij                                                   &
             + ksupuijf(:,:,n1:n2)*bsupvijf(:,:,n1:n2)                  &
             - ksupvijf(:,:,n1:n2)*bsupuijf(:,:,n1:n2)
      KxBuij = KxBuij                                                   &
             + ksupvijf(:,:,n1:n2)*bsupsijf(:,:,n1:n2)                  &
             - ksupsijf(:,:,n1:n2)*bsupvijf(:,:,n1:n2)
      KxBvij = KxBvij                                                   &
             + ksupsijf(:,:,n1:n2)*bsupuijf(:,:,n1:n2)                  &
             - ksupuijf(:,:,n1:n2)*bsupsijf(:,:,n1:n2)
      
      END SUBROUTINE Lorentz

!>  \brief Compute covariant (sub) components of the MHD force J X B - grad p
!>  The MHD forces are sum of Lorentz and pressure gradient forces
!>
!>     F = J X B - Grad(p)                                                 (1)
!>
!>  where J is the plasma current density. In curvilinear coodinates, the forces
!>  components becomes
!>
!>     F_s = (K^u*B^v - K^v*B^u) - dp/ds                                   (2a)
!>
!>     F_u = (K^v*B^s - K^s*B^v) - dp/du                                   (2b)
!>
!>     F_v = (K^s*B^u - K^u*B^s) - dp/dv                                   (2c)
!>
!>  where K^i = sqrt(g)*J^i . 
!>  Note this routine works on one parity at a time.
      SUBROUTINE GetMHDForce(fsubsmnf, fsubumnf, fsubvmnf, ipar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), POINTER, DIMENSION(:,:,:) ::                           &
         fsubsmnf, fsubumnf, fsubvmnf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  ::  fours, fouru, fourv, n, i, j, k
      REAL(dp), POINTER, DIMENSION(:,:,:)  :: KxBsmn, KXBumn, KXBvmn
      LOGICAL  ::  ltest
!----------------------------------------------- 
!
!Calculate harmonics of MHD force (F = JxB - grad p) on full mesh
!For l_linearize=T:  PERTURBED Force (ignore F(xc=0) part)
!For l_linearize=F:  FULL NONLINEAR Force (F(xc=0) part only)
!
      IF (ipar .EQ. isym) THEN
         fours = f_cos; fouru = f_sin; fourv = f_sin
      ELSE
         fours = f_sin; fouru = f_cos; fourv = f_cos
      END IF
      
      KXBsmn => fsubsmnf(:,:,nsmin:nsmax)
      KXBumn => fsubumnf(:,:,nsmin:nsmax)
      KXBvmn => fsubvmnf(:,:,nsmin:nsmax)

!Harmonics of covariant Lorentz force components on full grid
      CALL tomnsp_par(KxBsij, KxBsmn, fours)
      CALL tomnsp_par(KxBuij, KxBumn, fouru)
      CALL tomnsp_par(KxBvij, KxBvmn, fourv)

!Harmonics of "true" pressure (jacobian factor has been divided out of pijh)
      CALL tomnsp_par(pijh, pmnh, fours)
      CALL to_full_mesh(pmnh, pmnf)

!Radial derivative of pressure
      CALL GradientFull(pmnf_ds, pmnh)

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmnf_ds(1)
!            NEED pmn(m1,:,1) so limit at r12=hs_i/2 for dp/du, dp/dv is numerically correct
      IF (nsmin .EQ. 1) THEN
         pmnf(:,:,1) = 0; pmnf_ds(:,:,1) = 0
         pmnf(m0:m1,:,1) = pmnh(m0:m1,:,2)
         pmnf_ds(m1,:,1) = 2*ohs*pmnh(m1,:,2)
      END IF
      
      IF (nsmax.EQ.ns .AND. l_pedge) THEN
         pmnf(:,:,ns) = 0
         pmnf(m0,0,ns) = pmnh(m0,0,ns)
      END IF

      CALL get_force_harmonics (KxBsmn, KxBumn, KxBvmn, ipar)

      END SUBROUTINE GetMHDForce
      
      
!>  \brief Compute scaling factors for parallel flow damping
!>
      SUBROUTINE InitParDamping(parscale)
      USE hessian, ONLY:  l_Compute_Hessian, mupar, mupar_norm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: parscale                 !Scratch array
!-----------------------------------------------
      IF (muPar.EQ.zero .OR. .NOT.l_Compute_Hessian) RETURN

      CALL ASSERT(LBOUND(parscale,3).EQ.nsmin,'PARSCALE LBOUND WRONG!')
      CALL ASSERT(UBOUND(parscale,3).EQ.nsmax,'PARSCALE UBOUND WRONG!')
      ALLOCATE(pardamp(nsmin:nsmax,3))

!SPH 05-6-16 UPDATE mupar_norm IN cv_currents and put in 1/jacobh factor
!SIMPLER, FASTER APPROXIMATION: KEEP ONLY DIAGONAL DISPLACEMENT TERMS FROM v dot B
      muPar = ABS(muPar)

      parscale = bsubsijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
      CALL SurfAverage(pardamp(nsmin:nsmax,1), parscale, nsmin, nsmax)
      pardamp(:,1) = pardamp(:,1)*muPar*mupar_norm(nsmin:nsmax)

      parscale = bsubuijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
      CALL SurfAverage(pardamp(nsmin:nsmax,2), parscale, nsmin, nsmax)
      pardamp(:,2) = pardamp(:,2)*muPar*mupar_norm(nsmin:nsmax)
      
      parscale = bsubvijf(:,:,nsmin:nsmax)**2/jacobf(:,:,nsmin:nsmax)
      CALL SurfAverage(pardamp(nsmin:nsmax,3), parscale, nsmin, nsmax)
      pardamp(:,3) = pardamp(:,3)*muPar*mupar_norm(nsmin:nsmax)
     
      END SUBROUTINE InitParDamping

!>  \brief Final computation of the MHD covariant Fourier force components.
!>  Add pressure gradients and parallel damping (to Hessian), and scale 
!>  and print boundary values.
!>
!>
!>  \param[inout] f_smnf      F_s component of the MHD force. (cos, ipar=0)
!>  \param[inout] f_umnf      F_u component of the MHD force. (sin, ipar=0)
!>  \param[inout] f_vmnf      F_v component of the MHD force. (sin, ipar=0)
!>  \param[in]  ipar          =isym (lasym=F/T), =iasym (only for lasym=T)
!*******************************************************************************
       SUBROUTINE get_force_harmonics(f_smnf, f_umnf, f_vmnf, ipar) 
       USE hessian, ONLY: muPar, l_Compute_Hessian
       USE shared_data, ONLY: nprecon, l_PrintOriginForces,             &
                              l_push_s, l_push_u, l_push_v,             &
                              l_push_edge, col_scale
       USE nscalingtools, ONLY: FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
       
!these were computed in siesta_displacment and are used for || damping      
       USE quantities, ONLY: gvsupsmncf => fsupsmncf, gvsupumnsf => fsupumnsf,  &
             gvsupvmnsf => fsupvmnsf, gvsupsmnsf => fsupsmnsf,          &
             gvsupumncf => fsupumncf, gvsupvmncf => fsupvmncf
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,nsmin:nsmax)  ::           &
           f_smnf, f_umnf, f_vmnf
       INTEGER, INTENT(IN) ::  ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp), PARAMETER :: p5 = 0.5_dp
       REAL(dp)            :: skston, skstoff 
       REAL(dp)            :: f1, f2, f3, tmps(6)
       INTEGER             :: from, sparity, m, mp, n, np
       LOGICAL             :: ltest
!-----------------------------------------------
       CALL second0(skston)
       
      IF (ipar .EQ. isym) THEN
         sparity=1
      ELSE
         sparity=-1
      END IF

!add pressure gradient to Lorentz force
      f_smnf = f_smnf - pmnf_ds(:,:,nsmin:nsmax)
      
      DO m = 0, mpol
         mp = m*sparity
         f_umnf(m,:,:) = f_umnf(m,:,:) + mp*pmnf(m,:,nsmin:nsmax)         ! mp ~ -(d/du)p
      END DO
       
      DO n = -ntor, ntor       
         np = n*nfp*sparity
         f_vmnf(:,n,:) = f_vmnf(:,n,:) + np*pmnf(:,n,nsmin:nsmax)         ! np ~ -(d/dv)p
      ENDDO
       
       
      ltest = ALL(f_smnf(m0,-ntor:-1,:).EQ.zero)
      CALL ASSERT(ltest,'f_smnf != 0, m=0, n<0')
      ltest = ALL(f_umnf(m0,-ntor:-1,:).EQ.zero)
      CALL ASSERT(ltest,'f_umnf != 0, m=0, n<0')
      ltest = ALL(f_vmnf(m0,-ntor:-1,:).EQ.zero)
      CALL ASSERT(ltest,'f_vmnf != 0, m=0, n<0')
      IF (ipar .EQ. isym) THEN
         CALL AddParDamping(f_smnf, f_umnf, f_vmnf,                     &
                            gvsupsmncf, gvsupumnsf, gvsupvmnsf, ipar)
         ltest = ALL(f_umnf(m0,0,:).EQ.zero)
         CALL ASSERT(ltest,'f_umnf(m0,0) != 0')
         ltest = ALL(f_vmnf(m0,0,:).EQ.zero)
         CALL ASSERT(ltest,'f_vmnf(m0,0) != 0')
      ELSE
         CALL AddParDamping(f_smnf, f_umnf, f_vmnf,                     &
                            gvsupsmnsf, gvsupumncf, gvsupvmncf, ipar)
         ltest = ALL(f_smnf(m0,0,:).EQ.zero)
         CALL ASSERT(ltest,'f_smnf(m0,0) != 0')
      END IF

!GATHER BDY FORCES AND PRINT OUT ON PROC=0
       PRINT_O: IF (l_PrintOriginForces) THEN
          IF (nsmin .EQ. 1) THEN
          fbdy(1) = SQRT(SUM(f_smnf(m1,:,1)**2))
          fbdy(2) = SQRT(SUM(f_smnf(m1,:,2)**2))
          fbdy(3) = SQRT(SUM(f_smnf(m0,:,2)**2) + SUM(f_smnf(m2:,:,2)**2))

          fbdy(4) = SQRT(SUM(f_umnf(m1,:,1)**2))
          fbdy(5) = SQRT(SUM(f_umnf(m1,:,2)**2))
          fbdy(6) = SQRT(SUM(f_umnf(m0,:,2)**2) + SUM(f_umnf(m2:,:,2)**2))

          fbdy(7) = SQRT(SUM(f_vmnf(m0,:,1)**2))
          fbdy(8) = SQRT(SUM(f_vmnf(m0,:,2)**2))
          fbdy(9) = SQRT(SUM(f_vmnf(m1:,:,2)**2))
          END IF

          IF (nsmax .EQ. ns) THEN
          fbdy(11) = SQRT(SUM(f_umnf(:,:,ns)**2))
          fbdy(10) = SQRT(SUM(f_umnf(:,:,ns-1)**2))
          fbdy(13) = SQRT(SUM(f_vmnf(:,:,ns)**2))
          fbdy(12) = SQRT(SUM(f_vmnf(:,:,ns-1)**2))
          END IF

#if defined(MPI_OPT)
           NPROCS_0: IF (nprocs .GT. 1) THEN
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
!      p5 factors below are needed to preserve hessian symmetry
       NSMIN_1: IF (nsmin .EQ. 1) THEN
       f_smnf(m0,:,1) = 0; f_smnf(m2:,:,1) = 0
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
          f_vmnf(m0,:,1) = p5*f_vmnf(m0,:,1)
       END IF

       END IF NSMIN_1
 
!      EDGE BOUNDARY CONDITIONS
       EDGE: IF (nsmax .EQ. ns) THEN
          f_smnf(:,:,ns) = 0
          IF (l_pedge) f_umnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             f_umnf(:,:,ns) = 0
             f_vmnf(:,:,ns) = 0
          ELSE
             f_umnf(:,:,ns) = p5*f_umnf(:,:,ns)
             f_vmnf(:,:,ns) = p5*f_vmnf(:,:,ns)
          END IF
       END IF EDGE

       CALL second0(skstoff)
       get_force_harmonics_time=get_force_harmonics_time+(skstoff-skston)

       END SUBROUTINE get_force_harmonics


!>  \brief Add parallel flow damping terms to forces (for Hessian calculation)
!>   to suppress null space at resonances
!>
       SUBROUTINE addpardamping (f_smnf, f_umnf, f_vmnf,                &
                                 jvsupsmnf, jvsupumnf, jvsupvmnf, ipar)
       USE shared_data, ONLY:  col_scale
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor,nsmin:nsmax) ::            &
                 f_smnf, f_umnf, f_vmnf          
       REAL(dp), POINTER, DIMENSION(:,:,:) ::                           &
                 jvsupsmnf, jvsupumnf, jvsupvmnf                          !scaled in siesta_displacement
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER    :: js  !, mcol, fours, fouru, fourv
!-----------------------------------------------
!
!      Compute parallel (to B) damping force and add it to (linearized) forces
!      to suppress null space in preconditioner
!
       IF (.NOT.ALLOCATED(pardamp)) RETURN
       CALL ASSERT(LBOUND(f_smnf,3).EQ.nsmin,'LBOUND(f_smnf) != nsmin')
       CALL ASSERT(UBOUND(f_smnf,3).EQ.nsmax,'UBOUND(f_smnf) != nsmax')

!SPH 05-6-16 UPDATE mupar_norm IN cv_currents and put in 1/jacobh factor
!SIMPLER (FASTER, APPROXIMATELY DIAGONAL) FORM

       DO js = nsmin,nsmax
          f_smnf(:,:,js) = f_smnf(:,:,js)                               &
                         + jvsupsmnf(:,:,js)*pardamp(js,1)

          f_umnf(:,:,js) = f_umnf(:,:,js)                               &
                         + jvsupumnf(:,:,js)*pardamp(js,2)

          f_vmnf(:,:,js) = f_vmnf(:,:,js)                               &
                         + jvsupvmnf(:,:,js)*pardamp(js,3)
       END DO   


       END SUBROUTINE addpardamping

       END MODULE siesta_force
