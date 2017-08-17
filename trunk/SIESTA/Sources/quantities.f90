      MODULE quantities
!     
!     WRITTEN 08-18-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: This module contains the dynamical variables evolved: 
!        Namely, pressure, velocity and magnetic field.
!
      USE stel_kinds
      USE stel_constants      
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,            &
          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,        &
          ohs=>ohs_i, nfp=>nfp_i    !, mpol32=>mpol32_i, ntor32=>ntor32_i
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: lasym, lverbose
      IMPLICIT NONE

!
!       Convention: The name of the variables may contain:
!           i.   The first one or two letters tells the quantity. So, "p"=pressure, "B" = magnetic field,
!                "jB" = magnetic field*jacobian, "v" = velocity, "jv" = velocity*jacobian 
!           ii.  The next three letters maybe "SUB" or "SUP", to distinguish covariant/contravariant 
!           ii.  The next letter maybe "s", "u" or "v" to distinguish radial/poloidal/toroidal components
!           iii. The next two letters maybe "MN" or "IJ", to distinguish Fourier harmonics or angular values 
!           iv.  The next letter maybe "C" or "S" to distinguish COSINE or SINE parity
!           v.   The last letter maybe "F" or "H" to distinguish FULL or HALF radial mesh
!
!       Radial meshes
!           velocity, force components       :: full mesh
!           p, contravariant components of B :: half mesh
!           jacobh (jacobian to curvilinear coordinates) :: half mesh (stores <jacobh(2)> in jacobh(1))
!     VARIABLE DECLARATIONS  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp) :: b_factor, w_factor, p_factor                             ! Scale internal energy <wb+wp/(gam-1)> to 1
       REAL(dp) :: signjac
       REAL(dp) :: wp, wb
!---- Fourier Space Arrays - lasym=F parity
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &  ! Fourier harmonics= BSUP?, VSUB?, P
        jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,                      & 
        ksupsmnsf,  ksupumncf,  ksupvmncf,                               &
        djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,                  &  ! Array names starting with "d" represent perturbations
        fsupsmncf, fsupumnsf, fsupvmnsf                                  

!----  Pointers are subsections of xc, gc
       REAL(dp), POINTER, DIMENSION(:,:,:) ::                            &  ! Fourier harmonics= JVSUP?mn?F
        jvsupsmncf, jvsupumnsf, jvsupvmnsf,                              &  ! Evolved quantities
        fsubsmncf, fsubumnsf, fsubvmnsf                                     ! Computed forces

!---- Fourier Space Arrays - lasym=T parity
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &  ! Fourier harmonics= BSUP?, VSUB?, P
        jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh,                      & 
        ksupsmncf,  ksupumnsf,  ksupvmnsf,                               &
        djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh,                  &  ! Array names starting with "d" represent perturbations
        fsupsmnsf, fsupumncf, fsupvmncf                                  

       REAL(dp), POINTER, DIMENSION(:,:,:) ::                            &  ! Fourier harmonics= JVSUP?mn?F
        jvsupsmnsf, jvsupumncf, jvsupvmncf,                              &  ! Evolved quantities
        fsubsmnsf, fsubumncf, fsubvmncf                                     ! Computed forces

!---- Fourier Initial Power Spectrum Storage
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) ::                      &
        pwr_spec_s, pwr_spec_a
                                
!---- Real Space Arrays
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupsijf,             &  ! all these are needed to compute convolutions
        jvsupuijf, jvsupvijf, jacobh, jacobf, wint                          ! in pressure and magn. field evolution equations.
  
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
        bsupsijf0, bsupuijf0, bsupvijf0,                                 &
        bsupsijh0, bsupuijh0, bsupvijh0,                                 &
        bsubsijf,  bsubuijf,  bsubvijf, bsq,                             &
        ksubsijf, ksubuijf, ksubvijf,                                    &
        ksupsijf0, ksupuijf0, ksupvijf0

       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
        pijh, pijh_du, pijh_dv,                                          &
        pijf, pijf_ds, pijf_du, pijf_dv

       REAL(dp) :: fbdy(13)                
!
!    Note that: DBSUP?MNF:              calculated in UPDATE_BFIELD;    used in UPDATE_FORCE and UPDATE_STATE
!               DPCMNF:                 calculated in UPDATE_PRES;      used in UPDATE_FORCE and UPDATE_STATE
!               BSUP?MN:                updated by UPDATE_STATE;        used in UPDATE_FORCE   
!               PMN:                    updated in UPDATE_STATE;        used in UPDATE_FORCE     
!               JVSUP?IJ:               updated in UPDATE_UPPERV;       used in UPDATE_BFIELD + UPDATE_PRES 
!               FSUB?MN:                updated in UPDATE_FORCE;        used in UPDATE_VEL
!               VSUB?MN:                updated by UPDATE_VEL;          used in UPDATE_UPPERV
!-----------------------------------------------

      CONTAINS

        SUBROUTINE init_quantities
        USE fourier, ONLY: getorigin, tomnsp, toijsp, f_none, f_cos,     &
                           f_sin, f_sum
        USE metrics, ONLY: sqrtg
        USE read_wout_mod, ONLY: lasym_vmec=>lasym
        USE island_params, ONLY: nu_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0
        INTEGER :: istat, l
        REAL(dp) :: sum1
!-----------------------------------------------
        w_factor = 1._dp/ABS(wb_i)
        p_factor = w_factor
        b_factor = SQRT(w_factor)
        gnorm_i  = gnorm_i/w_factor

!
!       Allocate and initialize dependent variables
!
        CALL alloc_quantities

        jacobh(:,:,:) = RESHAPE(sqrtg, SHAPE(jacobh))                       ! Reshapes sqrtg in local form

!AVOID DIVIDE BY ZERO AND STORE AVERAGED VALUE (OVER THETA) 
!AT FIRST 1/2 GRID PT AT ORIGIN IN JACOBH(1)
        CALL getorigin (jacobh, m0, 0)
        signjac = jacobh(1,1,2)/ABS(jacobh(1,1,2))

!TO PRESERVE UNPERTURBED (CONSTANT IN s) PART OF PRESSURE, NYQUIST-FILTER JACOBIAN
!ON: SPH 04-28-16
!        IF (.FALSE.) THEN                                                  !SET TO FALSE TO TURN THIS OFF
        jacobf = jacobh
        CALL tomnsp(jacobh, jbsupsmnsh, f_cos)                      !use jbsupsmnsh for dummy storage
        IF (lasym) THEN
           CALL tomnsp(jacobh, jbsupsmnch, f_sin)
        END IF
        CALL toijsp(jbsupsmnsh, jacobh, f_none, f_cos)
        jbsupsmnsh = 0
        IF (lasym .and. .not.lasym_vmec) THEN          ! MRC If the VMEC case is stellarator
        CALL toijsp(jbsupsmnch, jacobh, f_sum, f_sin)  ! symmetric zero out the jacobian terms.
        jbsupsmnch = 0
        END IF
        
        IF (ANY(jacobh*signjac .LE. 0)) STOP 'FILTERED JACOBIAN CHANGED SIGN!'
        sum1 = SUM((jacobh(:,:,2:)-jacobf(:,:,2:))**2/     &
                   (jacobh(:,:,2:)+jacobf(:,:,2:))**2)
        IF (iam .EQ. 0 .and. lverbose) PRINT '(a,1pe8.2)',                &
           ' JACOBIAN SPECTRAL TRUNCATION ERROR: ', SQRT(sum1/SIZE(jacobh))
!        END IF

        CALL to_full_mesh(jacobh, jacobf, nuv)
!Must be consistent with variational form of pressure evolution equation (DO NOT USE 2pt EXTRAPOLATIONS!)
        jacobf(:,:,1) = jacobh(:,:,1)
        jacobf(:,:,ns)= jacobh(:,:,ns)

        DO l = 1, nu_i
            wint(l,:,:) = cosmui(l,m0)
        END DO         

        ALLOCATE (vp_f(ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error in init_bcovar'        
        CALL SurfAverage(vp_f, jacobf, 1, ns)

!SANITY CHECK: COMPARE TO volume_i
        sum1 = 0.5*hs_i*(SUM(vp_f(2:ns)) + SUM(vp_f(1:ns-1)))*signjac
        sum1 = (twopi*twopi)*sum1
        IF (ABS((sum1-volume_i)/volume_i) .GT. 1.E-3_dp) THEN

            IF (IAM .EQ. 0 .and. lverbose) PRINT *,'VMEC VOLUME: ',volume_i,' THIS VOLUME: ', sum1
!           STOP 'Volume does not agree with VMEC input!'
        END IF
      
        END SUBROUTINE init_quantities


        SUBROUTINE Init_Fields
        USE fourier, ONLY: orthonorm, tomnsp_par, f_sin, f_cos
        USE metrics, ONLY: lmns_i, lmnc_i, sqrtg, iflipj
        USE nscalingtools, ONLY: startglobrow, endglobrow
        USE mpi_inc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0, m1=1
        INTEGER  :: istat, m, n, nsmin, nsmax, nloc, js
        REAL(dp)           :: r0
        REAL(dp), DIMENSION(:), ALLOCATABLE :: phiph, chiph
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pmnf
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE:: presif, presih,        &
                                                  tmps, tmpu, tmpv, tmpp
!-----------------------------------------------

!
!       Initialize half-mesh covariant components of jac*Bsuph-i and jac*presh
!       

#if defined(SKS)
        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
#else
        nsmin = 1;  nsmax = ns
#endif
        jbsupsmnsh(:,:,nsmin:nsmax) = 0

        ALLOCATE(phiph(ns), chiph(ns))

        phipf_i = signjac*iflipj*phipf_i
        chipf_i = signjac*iflipj*chipf_i

        phiph(2:ns) = (phipf_i(2:ns) + phipf_i(1:nsh))/2
        chiph(2:ns) = (chipf_i(2:ns) + chipf_i(1:nsh))/2
        phiph(1) = 0;  chiph(1) = 0

!SKIP THIS TO USE INTERPOLATED LAMBDA FROM VMEC
!        IF (lasym) THEN
!           IF (iam .EQ. 0) PRINT *,'--- IMPLEMENT COMPUTE LMNC! ---'
!        ELSE
           CALL ReCompute_Lambda(lmns_i, lmnc_i, jacobh, orthonorm,       &
                                phiph, chiph, nsmin, nsmax)
!        ENDIF

        DO m = 0,mpol
           DO n = -ntor, ntor
               jbsupumnch(m,n,nsmin:nsmax) = phiph(nsmin:nsmax)*(-n*nfp*lmns_i(m,n,nsmin:nsmax))
               jbsupvmnch(m,n,nsmin:nsmax) = phiph(nsmin:nsmax)*( m    *lmns_i(m,n,nsmin:nsmax))
           END DO
        END DO

        jbsupumnch(m0,0,nsmin:nsmax) = chiph(nsmin:nsmax) / orthonorm(0,0)
        jbsupvmnch(m0,0,nsmin:nsmax) = phiph(nsmin:nsmax) / orthonorm(0,0)

        r0 = hs_i/2
        jbsupumnch(m1,:,1) = 0 
        jbsupvmnch(m1,:,1) = 0 
!        jbsupsmnsh(m1,:,2) = r0*jbsupumnch(m1,:,2)                        ! Need to make jbsupsmnsh - r0*jbsupumnch = 0 at js=2
                                                                           ! MESSES UP DIV(B) = 0
!STORE m=0 component of jbsupvmnch at origin in half-grid js=1
        jbsupvmnch(m0,:,1) = jbsupvmnch(m0,:,2) 
        jbsupvmnch(m1:,:,1) = 0

        jbsupumnch(:,:,nsmin:nsmax) = b_factor*jbsupumnch(:,:,nsmin:nsmax) 
        jbsupvmnch(:,:,nsmin:nsmax) = b_factor*jbsupvmnch(:,:,nsmin:nsmax) 
        IF (lasym) THEN
           jbsupsmnch = 0

           DO m = 0, mpol
              DO n = -ntor, ntor
                 jbsupumnsh(m,n,:) = phiph(:)*( n*nfp*lmnc_i(m,n,1:ns))
                 jbsupvmnsh(m,n,:) = phiph(:)*(-m    *lmnc_i(m,n,1:ns))
              END DO
           END DO

           jbsupumnsh(m1,:,1) = 0   
           jbsupvmnsh(m1,:,1) = 0  
!           jbsupsmnch(m1,:,2) = -r0*jbsupumnsh(m1,:,2)                    ! Need to make jbsupsmnch + r0*jbsupumnsh = 0 at js=2
                                                                           ! MESSES UP DIV(B) = 0
           jbsupumnsh = b_factor*jbsupumnsh
           jbsupvmnsh = b_factor*jbsupvmnsh

!STORE m=0 component of jbsupvmnch at origin in half-grid js=1
           jbsupvmnsh(m0,:,1) = jbsupvmnsh(m0,:,2)
           jbsupvmnsh(m1:,:,1) = 0
        END IF

#if defined(MRC)
        DEALLOCATE (phiph, chiph, stat=istat)
        IF (istat .NE. 0) STOP "Deallocate error #1 in init_fields"
#else
        DEALLOCATE (phiph, chiph, sqrtg, stat=istat)
        IF (istat .NE. 0) STOP "Deallocate error #1 in init_fields"
#endif

!
!       Initialize half mesh pressure
!
        ALLOCATE(presih(ntheta,nzeta,ns),                                &
                 presif(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP "Allocate error #1 in init_field"        

        DO js = MAX(1,nsmin-1), nsmax
           presif(:,:,js) = p_factor*presf_i(js)
        END DO

        nloc = MAX(1,startglobrow-1)

        CALL to_half_mesh_par(presif(:,:,nloc:nsmax),                    &
                              presih(:,:,nloc:nsmax), nuv, nloc, nsmax) 
      
!       Evolve (internally) jacob*(real pressure)
        presih(:,:,nsmin:nsmax) = presih(:,:,nsmin:nsmax) * jacobh(:,:,nsmin:nsmax)

        CALL tomnsp_par(presih, jpmnch, f_cos, nsmin, nsmax)               ! Initialize pressure harmonics (half mesh)
        jpmnch(:,:,1) = 0
        jpmnch(m0,:,1) = jpmnch(m0,:,2)

        IF (lasym) THEN
        CALL tomnsp_par(presih, jpmnsh, f_sin, nsmin, nsmax)
        jpmnsh(:,:,1) = 0
        jpmnsh(m0,:,1) = jpmnsh(m0,:,2)
        END IF


        DEALLOCATE(presif, presih, stat=istat)
        IF (istat .NE. 0) STOP "Deallocate error #2 in init_fields"

!
!       Gather jbsupmnX's, jpmn onto all processors (same as in update_state_par but for initial values)
!
#if defined(SKS)
        CALL GATHER_FIELDS(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch)
        IF (lasym)                                                       &
        CALL GATHER_FIELDS(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh)
#endif
        END SUBROUTINE Init_Fields


        SUBROUTINE ReCompute_Lambda(lmns, lmnc, jacobh, orthonorm, phiph, &
                                   chiph, nsmin, nsmax)
!-----------------------------------------------
        USE island_params, ns=>ns_i, nuv=>nuv_i, mnmax=>mnmax_i,  &
                           nfp => nfp_i, ntheta => nu_i, nzeta => nv_i  
        USE metrics, ONLY: guu, guv, gvv

        INTEGER, INTENT(IN)   :: nsmin, nsmax
        REAL(dp), INTENT(OUT) :: lmns(mnmax,ns), lmnc(mnmax,ns)                        
        REAL(dp), INTENT(IN)  :: phiph(ns), chiph(ns), orthonorm(mnmax), &
                                 jacobh(nuv, ns)
        INTEGER               :: js, m, n, mp, np, mn, mnp, lk, lu, lv,  &
                                 info, isign, n2
        REAL(dp), ALLOCATABLE :: amat(:,:), brhs(:,:), atemp(:), btemp(:)
        REAL(dp), ALLOCATABLE :: cosmni(:,:), cosmnp(:,:),               &
                                 sinmni(:,:), sinmnp(:,:)
        REAL(dp) :: ton, toff
!-----------------------------------------------
        CALL second0(ton)

        lmns(:,1) = 0
        IF (lasym) lmnc(:,1) = 0
!
!       ALLOCATE cos(mu+nv) and sin(mu+nv)
! 
        n2 = 1
        IF (lasym) n2 = 2
        ALLOCATE (atemp(nuv), btemp(nuv), cosmni(nuv,mnmax),             &
                  cosmnp(nuv,mnmax),                                     &
                  brhs(mnmax,n2), amat(mnmax*n2, mnmax*n2), stat=info)
        IF (info .NE. 0) STOP 'Allocation error in Recompute_Lambda'
        IF (lasym) THEN
        ALLOCATE (sinmni(nuv,mnmax), sinmnp(nuv,mnmax),                  &
                  stat=info)
        IF (info .NE. 0) STOP 'Allocation error in Recompute_Lambda'
        END IF
!COMPUTE MATRIX ELEMENTS, SURFACE BY SURFACE (can parallelize this: Sudip)

        mn = 0
        DO n=-ntor,ntor
           isign = 1; IF (n .LT. 0) isign = -1
           np = ABS(n)
        DO m=0,mpol
	      mn = mn+1
	      lk = 0
	      DO lv = 1, nzeta
	         DO lu = 1, ntheta
     	         lk = lk+1
                 cosmni(lk,mn) = cosmui(lu,m)*cosnv(lv,np)               &
                               - isign*sinmui(lu,m)*sinnv(lv,np)
                 cosmnp(lk,mn) = cosmu(lu,m)*cosnv(lv,np)                &
                               - isign*sinmu(lu,m)*sinnv(lv,np)
                 IF (.NOT. lasym) CYCLE
                 sinmni(lk,mn) = sinmui(lu,m)*cosnv(lv,np)               &
                               + isign*cosmui(lu,m)*sinnv(lv,np)
                 sinmnp(lk,mn) = sinmu(lu,m)*cosnv(lv,np)                &
                               + isign*cosmu(lu,m)*sinnv(lv,np)
                 END DO
              END DO
              IF (m.EQ.0 .AND. n.LT.0) THEN
                 cosmni(:,mn) = 0
                 cosmnp(:,mn) = 0
                 IF (lasym) THEN
                 sinmni(:,mn) = 0
                 sinmnp(:,mn) = 0
                 END IF
              END IF
           END DO
        END DO

!
!       COMPUTE Lambda from sqrt(g)*J^s == mB_v - nB_u = 0
!
        DO js = MAX(nsmin,2), nsmax

        mn = 0
        DO n=-ntor,ntor
	     DO m=0,mpol
!
!       RHS (driving) TERMS [sinmi, -cosmni parts]
!
	     mn = mn+1
	     btemp = (-m*(chiph(js)*guv(:,js)+phiph(js)*gvv(:,js))       &
	            + n*nfp*(chiph(js)*guu(:,js)+phiph(js)*guv(:,js)))   &
                   / jacobh(:,js)
             brhs(mn,1) = SUM(cosmni(:,mn)*btemp)
             IF (lasym) THEN
             brhs(mn,2) =-SUM(sinmni(:,mn)*btemp)
             END IF
			  
!
!            MATRIX ELEMENTS - COEFFICIENTS OF Lambda(mn,js)
!	     
	     mnp = 0

             DO np=-ntor,ntor
	        DO mp=0,mpol
	           mnp=mnp+1
                   atemp = (m*mp*gvv(:,js)+(n*nfp)*(np*nfp)*guu(:,js)    &
  	                 - nfp*(m*np+mp*n)*guv(:,js))/jacobh(:,js)
                   amat(mn,mnp) = SUM(cosmni(:,mn)*cosmnp(:,mnp)*atemp)
	           IF (mnp.EQ.mn .AND. amat(mn,mnp).EQ.zero)             &
                      amat(mn,mnp) = signjac
                   IF (.NOT.lasym) CYCLE
                   amat(mn,mnp+mnmax) =-SUM(cosmni(:,mn)*sinmnp(:,mnp)*atemp)
                   amat(mn+mnmax,mnp) =-SUM(sinmni(:,mn)*cosmnp(:,mnp)*atemp)
                   amat(mn+mnmax,mnp+mnmax)                              &
                                      = SUM(sinmni(:,mn)*sinmnp(:,mnp)*atemp)
	           IF (mnp.EQ.mn .AND. amat(mn+mnmax,mnp+mnmax).EQ.zero) &
                      amat(mn+mnmax,mnp+mnmax) = signjac
             	END DO
	     END DO
	     END DO
	  END DO

	  CALL solver(amat, brhs, mnmax*n2, 1, info)
	  IF (info .NE. 0) STOP 'INFO != 0 IN RECOMPUTE_LAMBDA'
	  lmns(:,js) = brhs(:,1)/orthonorm/phiph(js)                     !Divide by phiph -> old style
          IF (lasym)                                                    &
	  lmnc(:,js) = brhs(:,2)/orthonorm/phiph(js)                     

      END DO

      IF (mn.NE.mnmax .OR. mnp.NE.mnmax) STOP 'mn or mnp != mnmax in RECOMPUTE_LAMBDA'

      DEALLOCATE(atemp, btemp, cosmni, cosmnp, brhs, amat)
      IF (lasym) DEALLOCATE(sinmni, sinmnp)

      CALL second0(toff)
      IF (iam .EQ. 0 .and. lverbose) PRINT 100,' RECOMPUTE_LAMBDA - TIME: ', toff-ton,' S'

 100  FORMAT(a,f10.2,a)
      END SUBROUTINE Recompute_Lambda



        SUBROUTINE alloc_quantities
!
!       Allocates variables that will be needed through the advancing subroutines. 
!       More may be added in the future. 
!       All variables allocated here are then deallocated in DEALLOC_QUANTITIES, 
!       which should be called at the end of the run.
! 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: istat
!-----------------------------------------------
!       Half mesh quantities (harmonics) evolved: 3 magnetic field components and pressure                
        ALLOCATE(jbsupsmnsh(0:mpol,-ntor:ntor,ns),                       &  ! Half mesh quantities (harmonics)
                 jbsupumnch(0:mpol,-ntor:ntor,ns),                       &  ! B^s (sine), B^u (cosine), B^v (cosine)
                 jbsupvmnch(0:mpol,-ntor:ntor,ns),                       & 
                 jpmnch    (0:mpol,-ntor:ntor,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'      ! CONTRAVARIANT MAGNETIC FIELDS
        jbsupsmnsh = zero; jbsupumnch = zero; jbsupvmnch = zero        
        jpmnch = zero

        ALLOCATE (pwr_spec_s(0:mpol,-ntor:ntor,ns,4),                    &
                  pwr_spec_a(0:mpol,-ntor:ntor,ns,4), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'      ! CONTRAVARIANT MAGNETIC FIELDS


        ALLOCATE(fsupsmncf(0:mpol,-ntor:ntor,ns),                        &  ! Full mesh quantities (harmonics)
                 fsupumnsf(0:mpol,-ntor:ntor,ns),                        &  ! F_s (cosine), F_u (sine), F_v (sine)
                 fsupvmnsf(0:mpol,-ntor:ntor,ns), stat=istat) 
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'      ! CONTRAVARIANT MHD FORCES

        IF (lasym) THEN

        ALLOCATE(jbsupsmnch(0:mpol,-ntor:ntor,ns),                       &  ! Half mesh quantities (harmonics)
                 jbsupumnsh(0:mpol,-ntor:ntor,ns),                       &  ! B^s (sine), B^u (cosine), B^v (cosine)
                 jbsupvmnsh(0:mpol,-ntor:ntor,ns),                       & 
                 jpmnsh    (0:mpol,-ntor:ntor,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'      ! CONTRAVARIANT MAGNETIC FIELDS
        jbsupsmnch = zero; jbsupumnsh = zero; jbsupvmnsh = zero        
        jpmnsh = zero

        ALLOCATE(fsupsmnsf(0:mpol,-ntor:ntor,ns),                        &  ! Full mesh quantities (harmonics)
                 fsupumncf(0:mpol,-ntor:ntor,ns),                        &  ! F_s (cosine), F_u (sine), F_v (sine)
                 fsupvmncf(0:mpol,-ntor:ntor,ns), stat=istat) 
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'      ! CONTRAVARIANT MHD FORCES

        ENDIF
                
        ALLOCATE (bsq(ntheta,nzeta,ns),                                  &
                  jacobh(ntheta,nzeta,ns), jacobf(ntheta,nzeta,ns),      &
                  wint(ntheta,nzeta,ns), stat=istat)

        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'

        END SUBROUTINE alloc_quantities
        
        
        SUBROUTINE dealloc_quantities
        USE metrics, ONLY: lmns_i, lmnc_i
!
!       Deallocates all variables allocated in ALLOC_QUANTITIES
!        
        INTEGER :: istat
        
        DEALLOCATE(jpmnch, jbsupsmnsh, jbsupumnch, jbsupvmnch, djpmnch,  &
          djbsupsmnsh, djbsupumnch, djbsupvmnch,                         &
          fsupsmncf, fsupumnsf, fsupvmnsf,                               &
          ksupsmnsf, ksupumncf, ksupvmncf, pwr_spec_s, pwr_spec_a,       &
          stat=istat)
        IF(ALLOCATED(lmns_i)) DEALLOCATE(lmns_i)

        IF (lasym) THEN
        DEALLOCATE(jpmnsh, jbsupsmnch, jbsupumnsh, jbsupvmnsh, djpmnsh,  &
          djbsupsmnch, djbsupumnsh, djbsupvmnsh,                         &
          fsupsmnsf, fsupumncf, fsupvmncf,                               &
          ksupsmncf, ksupumnsf, ksupvmnsf, stat=istat)
        IF(ALLOCATED(lmnc_i)) DEALLOCATE(lmnc_i)
        END IF

        DEALLOCATE(jacobh, jacobf, jvsupsijf, jvsupuijf, jvsupvijf,      &
          ksubsijf, ksubuijf, ksubvijf, pijh, pijf,                      &
          pijf_ds, bsupuijf0, vp_f, bsq, wint,                           &
          bsupvijf0, bsupsijf0, pijf_du, pijf_dv,                        &
          pijh_du, pijh_dv, stat=istat)                      
        IF(ALLOCATED(bsupsijf0)) DEALLOCATE(bsupsijf0, bsupuijf0, bsupvijf0)
        IF(ALLOCATED(bsupsijh0)) DEALLOCATE(bsupsijh0, bsupuijh0, bsupvijh0)
        IF(ALLOCATED(ksupsijf0)) DEALLOCATE(ksupsijf0, ksupuijf0, ksupvijf0)
        IF(ALLOCATED(bsubsijf))  DEALLOCATE(bsubsijf,  bsubuijf, bsubvijf)
 
        IF (istat .NE. 0) STOP 'Deallocation failed in Quantities!'
        

        END SUBROUTINE dealloc_quantities
        

      SUBROUTINE GATHER_FIELDS(djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh)
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE hessian, ONLY: gather_array
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   ::                       &
                djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tmps, tmpu, tmpv, tmpp
      INTEGER  :: istat, nsmin, nsmax, k1, k2
!-----------------------------------------------
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      ALLOCATE (tmps(0:mpol,-ntor:ntor,ns), tmpu(0:mpol,-ntor:ntor,ns),  &
                tmpv(0:mpol,-ntor:ntor,ns), tmpp(0:mpol,-ntor:ntor,ns),  &
                stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error #1 in update_state_par'

      CALL gather_array(djbsupsmnh(:,:,nsmin:nsmax), tmps)
      CALL gather_array(djbsupumnh(:,:,nsmin:nsmax), tmpu)
      CALL gather_array(djbsupvmnh(:,:,nsmin:nsmax), tmpv)
      CALL gather_array(djpmnh(:,:,nsmin:nsmax), tmpp)

      k1 = LBOUND(djpmnh,3);  k2 = UBOUND(djpmnh,3)
      IF (k1.NE.1 .OR. k2.NE.ns) THEN
         DEALLOCATE(djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh)
         ALLOCATE (djbsupsmnh(0:mpol,-ntor:ntor,ns),                     &
                   djbsupumnh(0:mpol,-ntor:ntor,ns),                     &
                   djbsupvmnh(0:mpol,-ntor:ntor,ns),                     &
                   djpmnh(0:mpol,-ntor:ntor,ns), stat=istat)
         IF (istat .NE. 0) STOP 'Allocation error #2 in update_state_par'
      END IF

      djbsupsmnh = tmps
      djbsupumnh = tmpu
      djbsupvmnh = tmpv
      djpmnh     = tmpp
      DEALLOCATE (tmps, tmpu, tmpv, tmpp)

      END SUBROUTINE GATHER_FIELDS
          
     
      END MODULE quantities
          
