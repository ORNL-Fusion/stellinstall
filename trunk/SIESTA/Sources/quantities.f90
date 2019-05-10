!>  \brief Module contained subroutines for allocating and initializing
!>  curvilinear magnetic covariant and pressure Fourier harmonics on the
!>  half-radial mesh
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
      MODULE quantities
!     
!     WRITTEN 08-18-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: This module contains the dynamical variables evolved: 
!        Namely, pressure, velocity and magnetic field.
!
      USE stel_kinds
      USE stel_constants      
      USE v3_utilities, ONLY: assert
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,           &
          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,       &
          ohs=>ohs_i, nfp=>nfp_i    !, mpol32=>mpol32_i, ntor32=>ntor32_i
      USE descriptor_mod, ONLY: iam, nprocs, SIESTA_COMM
      USE shared_data, ONLY: lasym, lverbose
      USE nscalingtools, ONLY: startglobrow, endglobrow, PARSOLVER, MPI_ERR
      USE utilities, ONLY: to_half_mesh, to_full_mesh
      USE siesta_namelist, ONLY: l_vessel
      IMPLICIT NONE

!
!       Convention: The name of the variables may contain:
!           i.   The first one or two letters identifies the quantity. So, "p"=pressure, "B" = magnetic field,
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
!           jacobh (jacobian to curvilinear coordinates) :: half mesh
!                  store m0 (poloidal average) component of jacobh(2) in jacobh(1)
!           jacobf : full mesh average of jacobh, with jacobf(1) = jacobh(1)
!                    so values of p and B at origin can be extrapolated
!     VARIABLE DECLARATIONS  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PRIVATE :: isym=0, iasym=1
       REAL(dp) :: b_factor, w_factor, p_factor                            ! Scale internal energy <wb+wp/(gam-1)> to 1
       REAL(dp) :: signjac
       REAL(dp) :: wp, wb, wp0, wb0                                        ! summed (over procs) energies/local energies
!---- Fourier Space Arrays - lasym=F parity
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &  ! Fourier harmonics= BSUP?, VSUB?, P
        jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,                     & 
        ksupsmnsf,  ksupumncf,  ksupvmncf,                              &
        djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch                     ! Array names starting with "d" represent perturbations

!----  Pointers are subsections of xc, gc
       REAL(dp), POINTER, DIMENSION(:,:,:) ::                           &  ! Fourier harmonics= JVSUP?mn?F
        jvsupsmncf, jvsupumnsf, jvsupvmnsf,                             &  ! Evolved quantities
        fsubsmncf, fsubumnsf, fsubvmnsf,                                &  ! Computed forces
        fsupsmncf, fsupumnsf, fsupvmnsf                                  

!---- Fourier Space Arrays - lasym=T parity
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &  ! Fourier harmonics= BSUP?, VSUB?, P
        jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh,                     & 
        ksupsmncf,  ksupumnsf,  ksupvmnsf,                              &
        djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh                     ! Array names starting with "d" represent perturbations

       REAL(dp), POINTER, DIMENSION(:,:,:) ::                           &  ! Fourier harmonics= JVSUP?mn?F
        jvsupsmnsf, jvsupumncf, jvsupvmncf,                             &  ! Evolved quantities
        fsubsmnsf, fsubumncf, fsubvmncf,                                &  ! Computed forces
        fsupsmnsf, fsupumncf, fsupvmncf                                  

!---- Fourier Initial Power Spectrum Storage
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) ::                     &
        pwr_spec_s, pwr_spec_a
                                
!---- Real Space Arrays
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupsijf,            &  ! all these are needed to compute convolutions
        jvsupuijf, jvsupvijf, jacobh, jacobf, wint                         ! in pressure and magn. field evolution equations.
  
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
        bsupsijf0, bsupuijf0, bsupvijf0,                                &
        bsupsijh0, bsupuijh0, bsupvijh0,                                &
        bsupsijf,  bsupuijf,  bsupvijf,                                 &
        bsubsijf,  bsubuijf,  bsubvijf, bsq,                            &
        ksubsijf, ksubuijf, ksubvijf,                                   &
        ksupsijf0, ksupuijf0, ksupvijf0,                                &
        ksupsijf,  ksupuijf,  ksupvijf

       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
        pijh0, pijh0_du, pijh0_dv,                                      &
        pijf0, pijf0_ds

       REAL(dp) :: fbdy(13)                
!
!    Note that: DBSUP?MNF:              calculated in UPDATE_BFIELD;    used in UPDATE_FORCE and UPDATE_STATE
!               DPCMNF:                 calculated in UPDATE_PRES;      used in UPDATE_FORCE and UPDATE_STATE
!               BSUP?MN:                updated by UPDATE_STATE;        used in UPDATE_FORCE   
!               PMN:                    updated in UPDATE_STATE;        used in UPDATE_FORCE     
!               JVSUP?IJ:               updated in UPDATE_UPPERV;       used in UPDATE_BFIELD + UPDATE_PRES 
!               FSUB?MN:                updated in UPDATE_FORCE;        used in UPDATE_STATE
!               VSUB?MN:                updated by UPDATE_VEL;          used in UPDATE_DISPLACEMENT
!-----------------------------------------------

      CONTAINS
!*******************************************************************************
!>  \brief Initializes mpi and computes jacobian
!*******************************************************************************
        SUBROUTINE init_quantities
        USE descriptor_mod, ONLY: LSCALAPACK
        USE fourier, ONLY: getorigin, tomnsp, toijsp, f_none, f_cos,    &
                           f_sin, f_sum, Test_Fourier
        USE metrics, ONLY: sqrtg
        USE shared_data, ONLY: mblk_size, ndims, lasym
#if defined(MPI_OPT)        
        USE prof_mod, ONLY: SetUpScalingAllGather
#endif        
#if defined(SKS)
        USE blocktridiagonalsolver_s, ONLY: Initialize
#endif
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
        
!Store dimensions
        IF (lasym) THEN
           ndims = 6
        ELSE
           ndims = 3
        END IF      
                
        mblk_size=ndims*(mpol+1)*(2*ntor+1)
        
        CALL ASSERT(mblk_size.NE.0, 'mblk_size = 0 in init_quantities')
        
        IF (PARSOLVER) THEN
#if defined(SKS)
           CALL Initialize(ns, mblk_size) 
#endif
        ELSE IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
        CALL SetUpScalingAllGather(mblk_size)
#endif
        ELSE
           startglobrow = 1; endglobrow = ns
        END IF
        
!Test orthogonality AFTER mpi is initialized        
        CALL Test_Fourier

!       Allocate and initialize dependent variables
        CALL alloc_quantities

        jacobh(:,:,:) = RESHAPE(sqrtg, SHAPE(jacobh))  ! Reshapes sqrtg in local form

!AVOID DIVIDE BY ZERO AND STORE AVERAGED VALUE (OVER THETA) 
!AT FIRST 1/2 GRID PT AT ORIGIN IN JACOBH(1)
        CALL getorigin (jacobh, m0)
        signjac = jacobh(1,1,2)/ABS(jacobh(1,1,2))

!TO PRESERVE UNPERTURBED (CONSTANT IN s) PART OF PRESSURE, NYQUIST-FILTER JACOBIAN
!ON: SPH 04-28-16
!        IF (.FALSE.) THEN
        CALL tomnsp(jacobh, jacobf, f_cos)                            !use jacobf for dummy storage
        CALL toijsp(jacobf, jacobh, f_none, f_cos)
        IF (lasym) THEN
           CALL tomnsp(jacobh, jacobf, f_sin)
           CALL toijsp(jacobf, jacobh, f_sum, f_sin)
        END IF
!       END IF        
        
        CALL ASSERT(ALL(jacobh*signjac.GT.0),'FILTERED JACOBIAN CHANGED SIGN!')
        sum1 = SUM((jacobh(:,:,2:)-jacobf(:,:,2:))**2/                  &
                   (jacobh(:,:,2:)+jacobf(:,:,2:))**2)
        IF (iam .EQ. 0 .AND. lverbose) PRINT '(a,1pe8.2)',              &
           ' JACOBIAN SPECTRAL TRUNCATION ERROR: ', SQRT(sum1/SIZE(jacobh))
!        END IF

        CALL to_full_mesh(jacobh, jacobf)
!Must be consistent with variational form of pressure evolution equation (DO NOT USE 2pt EXTRAPOLATIONS!)
        jacobf(:,:,1) = jacobh(:,:,1)
        jacobf(:,:,ns)= jacobh(:,:,ns)

        DO l = 1, ntheta
            wint(l,:,:) = cosmui(l,m0)
        END DO         

        ALLOCATE (vp_f(ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation error in init_quantities')
        CALL SurfAverage(vp_f, jacobf, 1, ns)

!SANITY CHECK: COMPARE TO volume_i
        sum1 = 0.5*hs_i*(SUM(vp_f(2:ns)) + SUM(vp_f(1:ns-1)))*signjac
        sum1 = (twopi*twopi)*sum1
        IF (ABS((sum1-volume_i)/volume_i) .GT. 1.E-3_dp) THEN
            IF (IAM .EQ. 0 .and. lverbose) PRINT *,'VMEC VOLUME: ',     &
               volume_i,' SIESTA VOLUME: ', sum1
        END IF
        END SUBROUTINE init_quantities

!*******************************************************************************
!>  \brief Initializes half-mesh covariant components of jac*Bsuph-i and jac*presh
!>   Fields are gathered to all processors at the end of this routine
!>   Only called if lrestart=F (otherwise data read in from restart file)    
!*******************************************************************************
        SUBROUTINE Init_Fields
        USE fourier, ONLY: orthonorm, tomnsp_par, f_sin, f_cos
        USE metrics, ONLY: lmns_i, lmnc_i, sqrtg, iflipj
        USE mpi_inc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0, m1=1
        INTEGER  :: istat, nsmin, nsmax, nloc, js
        REAL(dp), DIMENSION(:), ALLOCATABLE :: phiph, chiph
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pmnf
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE:: presif, presih
!-----------------------------------------------
       
        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

        ALLOCATE(phiph(ns), chiph(ns))

        phipf_i = signjac*iflipj*phipf_i
        chipf_i = signjac*iflipj*chipf_i

        phiph(2:ns) = (phipf_i(2:ns) + phipf_i(1:nsh))/2
        chiph(2:ns) = (chipf_i(2:ns) + chipf_i(1:nsh))/2
        phiph(1) = 0;  chiph(1) = 0

        IF (.NOT.l_vessel) THEN
          CALL ReCompute_Lambda(lmns_i, lmnc_i, jacobh, orthonorm,      &
                                phiph, chiph, nsmin, nsmax)
        END IF
         
  		CALL Init_Bfield(jbsupsmnsh, jbsupumnch, jbsupvmnch,            &
                         lmns_i, phiph, chiph, nsmin, nsmax, isym)
		IF (lasym)                                                      &
		    CALL Init_Bfield(jbsupsmnch, jbsupumnsh, jbsupvmnsh,        &
                             lmnc_i, phiph, chiph, nsmin, nsmax, iasym)

        DEALLOCATE (phiph, chiph, stat=istat)
        CALL ASSERT(istat.EQ.0,'Deallocate error #1 in init_fields')

!
!       Initialize half mesh pressure
!
        nloc = MAX(1,startglobrow-1)
        ALLOCATE(presih(ntheta,nzeta,nloc:nsmax),                       &
                 presif(ntheta,nzeta,nloc:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0, 'Allocate error #1 in init_fields')

        DO js = nloc, nsmax
           presif(:,:,js) = p_factor*presf_i(js)
        END DO

        CALL to_half_mesh(presif, presih)
      
! Initialize (internally) jacob*(real pressure), which is evolved
        presih(:,:,nsmin:nsmax) = presih(:,:,nsmin:nsmax)*jacobh(:,:,nsmin:nsmax)

! Initialize internal pressure harmonics (half mesh)        
        CALL tomnsp_par(presih(:,:,nsmin:nsmax), jpmnch(:,:,nsmin:nsmax), &
                        f_cos)
        jpmnch(:,:,1) = 0
        jpmnch(m0,:,1) = jpmnch(m0,:,2)

        IF (lasym) THEN
           CALL tomnsp_par(presih(:,:,nsmin:nsmax), jpmnsh(:,:,nsmin:nsmax), &
                           f_sin)
           jpmnsh(:,:,1) = 0
           jpmnsh(m0,:,1) = jpmnsh(m0,:,2)
        END IF

        DEALLOCATE(presif, presih, stat=istat)
        CALL ASSERT(istat.EQ.0, 'Deallocate error #2 in init_fields')

!
!       Gather jbsupmnX's, jpmn onto all processors (same as in update_state but for initial values)
!
        CALL GatherFields

        END SUBROUTINE Init_Fields


        SUBROUTINE Init_Bfield(jbsupsmnh, jbsupumnh, jbsupvmnh,         &
                               lmn, phiph, chiph, nsmin, nsmax, ipar)
        USE fourier, ONLY: orthonorm
        USE siesta_namelist, ONLY: l_vessel
#if defined(MPI_OPT)        
        USE mpi_inc
#endif        
!-----------------------------------------------
	    INTEGER, INTENT(IN)  :: ipar, nsmin, nsmax
		REAL(dp), INTENT(IN) :: phiph(ns), chiph(ns)
		REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(IN) :: lmn
		REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) ::       &
		                        jbsupsmnh, jbsupumnh, jbsupvmnh
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0, m1=1
        INTEGER  :: m, n, sparity, mp, np
!-----------------------------------------------

		IF (ipar .EQ. isym) THEN
		   sparity = 1
		ELSE
		   sparity = -1
		END IF
		
        IF (.NOT.l_vessel) THEN
      
        DO m = 0,mpol
	  	   mp = m*sparity
           DO n = -ntor, ntor
		      np = n*nfp*sparity
              jbsupumnh(m,n,nsmin:nsmax) = phiph(nsmin:nsmax)*(-np*lmn(m,n,nsmin:nsmax))
              jbsupvmnh(m,n,nsmin:nsmax) = phiph(nsmin:nsmax)*( mp*lmn(m,n,nsmin:nsmax))
           END DO
        END DO

        IF (ipar .EQ. isym) THEN
           jbsupumnh(m0,0,nsmin:nsmax) = chiph(nsmin:nsmax) / orthonorm(0,0)
           jbsupvmnh(m0,0,nsmin:nsmax) = phiph(nsmin:nsmax) / orthonorm(0,0)
        END IF
        
        END IF
        
        jbsupsmnh(:,:,nsmin:nsmax) = b_factor*jbsupsmnh(:,:,nsmin:nsmax)
        jbsupumnh(:,:,nsmin:nsmax) = b_factor*jbsupumnh(:,:,nsmin:nsmax)
        jbsupvmnh(:,:,nsmin:nsmax) = b_factor*jbsupvmnh(:,:,nsmin:nsmax)


!STORE m=0 component of jbsupvmnh (and m=0,1 of jbsupumnh)  at origin of half-grid js=1
        IF (nsmin .EQ. 1) THEN
           jbsupsmnh(:,:,1) = 0;  jbsupsmnh(m1,:,1) = jbsupsmnh(m1,:,2)
           jbsupumnh(:,:,1) = 0;  jbsupumnh(m0:m1,:,1) = jbsupumnh(m0:m1,:,2)
           jbsupvmnh(m0,:,1) = jbsupvmnh(m0,:,2) 
           jbsupvmnh(m1:,:,1) = 0
        END IF
        
		END SUBROUTINE Init_Bfield

!*******************************************************************************
!>  \brief recomputes "lambda" for VMEC equilibrium on sqrt(flux) mesh
!>   by solving J^s = 0 in flux coordinates
!*******************************************************************************
        SUBROUTINE ReCompute_Lambda(lmns, lmnc, jacobh, orthonorm,      &
                                    phiph, chiph, nsmin, nsmax)
!-----------------------------------------------
        USE island_params, ns=>ns_i, nuv=>nuv_i, mnmax=>mnmax_i,  &
                           nfp => nfp_i, ntheta => nu_i, nzeta => nv_i  
        USE metrics, ONLY: guu, guv, gvv

        INTEGER, INTENT(IN)   :: nsmin, nsmax
        REAL(dp), INTENT(OUT) :: lmns(mnmax,ns), lmnc(mnmax,ns)
        REAL(dp), INTENT(IN)  :: phiph(ns), chiph(ns), orthonorm(mnmax),&
                                 jacobh(nuv, ns)
        INTEGER               :: js, m, n, mp, np, mn, mnp, lk, lu, lv, &
                                 info, isign, n2
        REAL(dp), ALLOCATABLE :: amat(:,:), brhs(:,:), atemp(:), btemp(:)
        REAL(dp), ALLOCATABLE :: cosmni(:,:), cosmnp(:,:),              &
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
        ALLOCATE (atemp(nuv), btemp(nuv), cosmni(nuv,mnmax),            &
                  cosmnp(nuv,mnmax),                                    &
                  brhs(mnmax,n2), amat(mnmax*n2, mnmax*n2), stat=info)
        CALL ASSERT(info.EQ.0,'Allocation error in Recompute_Lambda')
        IF (lasym) THEN
        ALLOCATE (sinmni(nuv,mnmax), sinmnp(nuv,mnmax), stat=info)
        CALL ASSERT(info.EQ.0,'Allocation error in Recompute_Lambda')
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
                 cosmni(lk,mn) = cosmui(lu,m)*cosnv(lv,np)              &
                               - isign*sinmui(lu,m)*sinnv(lv,np)
                 cosmnp(lk,mn) = cosmu(lu,m)*cosnv(lv,np)               &
                               - isign*sinmu(lu,m)*sinnv(lv,np)
                 IF (.NOT. lasym) CYCLE
                 sinmni(lk,mn) = sinmui(lu,m)*cosnv(lv,np)              &
                               + isign*cosmui(lu,m)*sinnv(lv,np)
                 sinmnp(lk,mn) = sinmu(lu,m)*cosnv(lv,np)               &
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
	          btemp = (-m*(chiph(js)*guv(:,js)+phiph(js)*gvv(:,js))     &
	            + n*nfp*(chiph(js)*guu(:,js)+phiph(js)*guv(:,js)))      &
                   / jacobh(:,js)
              brhs(mn,1) = SUM(cosmni(:,mn)*btemp)
              IF (lasym) brhs(mn,2) =-SUM(sinmni(:,mn)*btemp)
			  
!
!            MATRIX ELEMENTS - COEFFICIENTS OF Lambda(mn,js)
!	     
              mnp = 0

              DO np=-ntor,ntor
                 DO mp=0,mpol
                    mnp=mnp+1
                    atemp = (m*mp*gvv(:,js)+(n*nfp)*(np*nfp)*guu(:,js)  &
  	                      - nfp*(m*np+mp*n)*guv(:,js))/jacobh(:,js)
                    amat(mn,mnp) = SUM(cosmni(:,mn)*cosmnp(:,mnp)*atemp)
	                IF (mnp.EQ.mn .AND. amat(mn,mnp).EQ.zero)           &
                       amat(mn,mnp) = signjac
                    IF (.NOT.lasym) CYCLE
                    amat(mn,mnp+mnmax) =-SUM(cosmni(:,mn)*sinmnp(:,mnp)*atemp)
                    amat(mn+mnmax,mnp) =-SUM(sinmni(:,mn)*cosmnp(:,mnp)*atemp)
                    amat(mn+mnmax,mnp+mnmax)                            &
                                       = SUM(sinmni(:,mn)*sinmnp(:,mnp)*atemp)
	                IF (mnp.EQ.mn .AND. amat(mn+mnmax,mnp+mnmax).EQ.zero) &
                       amat(mn+mnmax,mnp+mnmax) = signjac
                 END DO
	          END DO
	       END DO
        END DO

        CALL solver(amat, brhs, mnmax*n2, 1, info)
        CALL ASSERT(info.EQ.0,'INFO != 0 IN RECOMPUTE_LAMBDA')
        lmns(:,js) = brhs(:,1)/orthonorm/phiph(js)                        !Divide by phiph -> old style
        IF (lasym) lmnc(:,js) = brhs(:,2)/orthonorm/phiph(js)

        END DO

        CALL ASSERT(mn.EQ.mnmax, mnp.EQ.mnmax, 'mn or mnp != mnmax in RECOMPUTE_LAMBDA')

        DEALLOCATE(atemp, btemp, cosmni, cosmnp, brhs, amat)
        IF (lasym) DEALLOCATE(sinmni, sinmnp)

        CALL second0(toff)
        IF (iam .EQ. 0 .and. lverbose) PRINT 100,' RECOMPUTE_LAMBDA - TIME: ', toff-ton,' S'

100     FORMAT(a,f10.2,a)
        END SUBROUTINE Recompute_Lambda

!*******************************************************************************
!>  \brief Allocates variables that will be updated in the time-advancing subroutines. 
!>       All variables allocated here are then deallocated in DEALLOC_QUANTITIES, 
!>       which should be called at the end of the run.
!*******************************************************************************
        SUBROUTINE alloc_quantities
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: istat
!-----------------------------------------------
! Half mesh quantities (harmonics) evolved: 3 magnetic field components and pressure                
! CONTRAVARIANT MAGNETIC FIELDS
        ALLOCATE(jbsupsmnsh(0:mpol,-ntor:ntor,ns),                      &  ! Half mesh quantities (harmonics)
                 jbsupumnch(0:mpol,-ntor:ntor,ns),                      &  ! B^s (sine), B^u (cosine), B^v (cosine)
                 jbsupvmnch(0:mpol,-ntor:ntor,ns),                      & 
                 jpmnch    (0:mpol,-ntor:ntor,ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #1 failed in ALLOC_QUANTITIES')
        jbsupsmnsh = zero; jbsupumnch = zero; jbsupvmnch = zero        
        jpmnch = zero

        ALLOCATE (pwr_spec_s(0:mpol,-ntor:ntor,ns,4),                   &
                  pwr_spec_a(0:mpol,-ntor:ntor,ns,4), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #2 failed in ALLOC_QUANTITIES')


        IF (lasym) THEN

        ALLOCATE(jbsupsmnch(0:mpol,-ntor:ntor,ns),                      &  ! Half mesh quantities (harmonics)
                 jbsupumnsh(0:mpol,-ntor:ntor,ns),                      &  ! B^s (sine), B^u (cosine), B^v (cosine)
                 jbsupvmnsh(0:mpol,-ntor:ntor,ns),                      & 
                 jpmnsh    (0:mpol,-ntor:ntor,ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #4 failed in ALLOC_QUANTITIES')
        jbsupsmnch = zero; jbsupumnsh = zero; jbsupvmnsh = zero        
        jpmnsh = zero

        ENDIF
                
        ALLOCATE (bsq(ntheta,nzeta,ns),                                 &
                  jacobh(ntheta,nzeta,ns), jacobf(ntheta,nzeta,ns),     &
                  wint(ntheta,nzeta,ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #6 failed in ALLOC_QUANTITIES')

        END SUBROUTINE alloc_quantities
        
        
        SUBROUTINE dealloc_quantities
        USE metrics, ONLY: lmns_i, lmnc_i
!
!       Deallocates all variables allocated in ALLOC_QUANTITIES
!        
        INTEGER :: istat
        
        DEALLOCATE(jpmnch, jbsupsmnsh, jbsupumnch, jbsupvmnch, djpmnch, &
          djbsupsmnsh, djbsupumnch, djbsupvmnch,                        &
          ksupsmnsf, ksupumncf, ksupvmncf, pwr_spec_s, pwr_spec_a,      &
          stat=istat)

        IF (lasym) THEN
        DEALLOCATE(jpmnsh, jbsupsmnch, jbsupumnsh, jbsupvmnsh, djpmnsh, &
          djbsupsmnch, djbsupumnsh, djbsupvmnsh,                        &
          ksupsmncf, ksupumnsf, ksupvmnsf, stat=istat)
        END IF

        DEALLOCATE(jacobh, jacobf, jvsupsijf, jvsupuijf, jvsupvijf,     &
          ksubsijf, ksubuijf, ksubvijf, pijh0, pijf0,                   &
          pijf0_ds, bsupuijf0, vp_f, bsq, wint,                         &
          bsupvijf0, bsupsijf0, pijh0_du, pijh0_dv, stat=istat)                      
        IF(ALLOCATED(bsupsijf0)) DEALLOCATE(bsupsijf0, bsupuijf0, bsupvijf0)
        IF(ALLOCATED(bsupsijh0)) DEALLOCATE(bsupsijh0, bsupuijh0, bsupvijh0)
        IF(ALLOCATED(ksupsijf0)) DEALLOCATE(ksupsijf0, ksupuijf0, ksupvijf0)
        IF(ALLOCATED(bsubsijf))  DEALLOCATE(bsubsijf,  bsubuijf, bsubvijf)
 
        CALL ASSERT(istat.EQ.0,'Deallocation failed in dealloc_quantities')
        

      END SUBROUTINE dealloc_quantities

!>  \brief Parallel/Serial routine to compute contravariant force components
!!  \author S. P. Hirshman
      SUBROUTINE toupper_forces
      USE island_params, ONLY: ntheta=>nu_i, nzeta=>nv_i,               &
                         mpol=>mpol_i, ntor=>ntor_i
      USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_cos, f_sin,  &
                         f_sum
      USE metrics, ONLY: toupper
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                                    fsubsijf, fsubuijf, fsubvijf,       &
                                    fsupsijf, fsupuijf, fsupvijf
!-----------------------------------------------
      INTEGER :: nsmin, nsmax, istat
!-----------------------------------------------
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)

      ALLOCATE(fsubsijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsubuijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsubvijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupsijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupuijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,' Allocation error in toupper_forces')

      CALL toijsp_par(fsubsmncf(:,:,nsmin:nsmax), fsubsijf, f_none, f_cos)
      CALL toijsp_par(fsubumnsf(:,:,nsmin:nsmax), fsubuijf, f_none, f_sin)
      CALL toijsp_par(fsubvmnsf(:,:,nsmin:nsmax), fsubvijf, f_none, f_sin)
 
      IF (lasym) THEN

      CALL toijsp_par(fsubsmnsf(:,:,nsmin:nsmax), fsubsijf, f_sum, f_sin)
      CALL toijsp_par(fsubumncf(:,:,nsmin:nsmax), fsubuijf, f_sum, f_cos)
      CALL toijsp_par(fsubvmncf(:,:,nsmin:nsmax), fsubvijf, f_sum, f_cos)
 
      END IF

!      
!     COMPUTE CONTRAVARIANT FORCES NEEDED TO COMPUTE <||F||^2>
!

      CALL toupper(fsubsijf, fsubuijf, fsubvijf,                        &
                   fsupsijf, fsupuijf, fsupvijf, nsmin, nsmax) 
      DEALLOCATE(fsubsijf, fsubuijf, fsubvijf)

      CALL tomnsp_par(fsupsijf, fsupsmncf(:,:,nsmin:nsmax), f_cos)
      CALL tomnsp_par(fsupuijf, fsupumnsf(:,:,nsmin:nsmax), f_sin)
      CALL tomnsp_par(fsupvijf, fsupvmnsf(:,:,nsmin:nsmax), f_sin)

      IF (lasym) THEN

      CALL tomnsp_par(fsupsijf, fsupsmnsf(:,:,nsmin:nsmax), f_sin)
      CALL tomnsp_par(fsupuijf, fsupumncf(:,:,nsmin:nsmax), f_cos)
      CALL tomnsp_par(fsupvijf, fsupvmncf(:,:,nsmin:nsmax), f_cos)

      END IF

      DEALLOCATE(fsupsijf, fsupuijf, fsupvijf, stat=istat)

      END SUBROUTINE toupper_forces


      SUBROUTINE GATHERFIELDS
      
      IF (PARSOLVER) THEN
      CALL GATHER_FIELDS(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch)
      IF (lasym)                                                        &
      CALL GATHER_FIELDS(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh)
      END IF
      
      END SUBROUTINE GATHERFIELDS  

#if 0
      SUBROUTINE PADFIELDS
      
      IF (PARSOLVER) THEN
      CALL PAD_FIELDS(djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch)
      IF (lasym)                                                        &
      CALL PAD_FIELDS(djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh)
      END IF
      
      END SUBROUTINE PADFIELDS      
#endif

      SUBROUTINE GATHER_FIELDS(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh)
      USE hessian, ONLY: gather_array
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   ::                      &
                jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tmps, tmpu, tmpv, tmpp
      INTEGER  :: istat, nsmin, nsmax, k1, k2
!-----------------------------------------------
!
!     gathers (total) field arrays for [1,ns] onto all processors
!
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      CALL ASSERT(ALLOCATED(jbsupsmnh), 'input arrays not allocated in gather_fields')

      ALLOCATE (tmps(0:mpol,-ntor:ntor,ns), tmpu(0:mpol,-ntor:ntor,ns), &
                tmpv(0:mpol,-ntor:ntor,ns), tmpp(0:mpol,-ntor:ntor,ns), &
                stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation error #1 in gather_fields')
      k1 = SIZE(tmps,1)*SIZE(tmps,2)

      CALL gather_array(jbsupsmnh(:,:,nsmin:nsmax), tmps)
      CALL gather_array(jbsupumnh(:,:,nsmin:nsmax), tmpu)
      CALL gather_array(jbsupvmnh(:,:,nsmin:nsmax), tmpv)
      CALL gather_array(jpmnh(:,:,nsmin:nsmax), tmpp)

      k1 = LBOUND(jpmnh,3);  k2 = UBOUND(jpmnh,3)
      IF (k1.NE.1 .OR. k2.NE.ns) THEN
         DEALLOCATE(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh, stat=istat)
         IF (istat .EQ. 0)                                              &
         ALLOCATE (jbsupsmnh(0:mpol,-ntor:ntor,ns),                     &
                   jbsupumnh(0:mpol,-ntor:ntor,ns),                     &
                   jbsupvmnh(0:mpol,-ntor:ntor,ns),                     &
                   jpmnh(0:mpol,-ntor:ntor,ns), stat=istat)
         CALL ASSERT(istat.EQ.0,'Allocation error #2 in gather_fields')
      END IF

      jbsupsmnh = tmps
      jbsupumnh = tmpu
      jbsupvmnh = tmpv
      jpmnh     = tmpp
      DEALLOCATE (tmps, tmpu, tmpv, tmpp)

      END SUBROUTINE GATHER_FIELDS

#if 0          
      SUBROUTINE PAD_FIELDS(djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh)
      USE blocktridiagonalsolver_s, ONLY: PadSides
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   ::                      &
                djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tmps, tmpu, tmpv, tmpp
      INTEGER  :: istat, nsmin, nsmax, k1, k2
!-----------------------------------------------
!
!     pads (perturbed) field arrays to [nsmin-1, nsmax+2] as required in siesta_init
!     From siesta_bfield, siesta_pressure, djbsupXmnh, djpmnh are on
!     [nsmin,nsmax+1] as required for Hessian computation in siesta_forces
!
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)

      CALL ASSERT(ALLOCATED(djbsupsmnh),'input arrays not allocated in pad_fields')

      ALLOCATE (tmps(0:mpol,-ntor:ntor,ns), tmpu(0:mpol,-ntor:ntor,ns),  &
                tmpv(0:mpol,-ntor:ntor,ns), tmpp(0:mpol,-ntor:ntor,ns),  &
                stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation error #1 in pad_fields')

      k1 = SIZE(tmps,1)*SIZE(tmps,2)
      tmps(:,:,nsmin:nsmax) = djbsupsmnh(:,:,nsmin:nsmax)
      CALL PadSides(tmps, k1, 1, 2)
      tmpu(:,:,nsmin:nsmax) = djbsupumnh(:,:,nsmin:nsmax)
      CALL PadSides(tmpu, k1, 1, 2)
      tmpv(:,:,nsmin:nsmax) = djbsupvmnh(:,:,nsmin:nsmax)
      CALL PadSides(tmpv, k1, 1, 2)
      tmpp(:,:,nsmin:nsmax) = djpmnh(:,:,nsmin:nsmax)
      CALL PadSides(tmpp, k1, 1, 2)

      k1 = LBOUND(djpmnh,3);  k2 = UBOUND(djpmnh,3)
 	  CALL ASSERT(k1.EQ.MAX(1,startglobrow-1),'k1 wrong in pad_fields')
	  CALL ASSERT(k2.EQ.MIN(ns,endglobrow+2),'k2 wrong in pad_fields')

      djbsupsmnh(:,:,k1:k2) = tmps(:,:,k1:k2)
      djbsupumnh(:,:,k1:k2) = tmpu(:,:,k1:k2)
      djbsupvmnh(:,:,k1:k2) = tmpv(:,:,k1:k2)
      djpmnh(:,:,k1:k2)     = tmpp(:,:,k1:k2)
      DEALLOCATE (tmps, tmpu, tmpv, tmpp)

      END SUBROUTINE PAD_FIELDS
#endif


      SUBROUTINE SurfAverage(average, q3d, nsmin, nsmax)
      USE stel_kinds, ONLY: dp
      USE island_params, ONLY: ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: nsmin, nsmax
      REAL(dp), INTENT(OUT) :: average(nsmin:nsmax)
      REAL(dp), INTENT(IN)  :: q3d(ntheta,nzeta,nsmin:nsmax)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: js
!-----------------------------------------------

      DO js = nsmin, nsmax
		 average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
      END DO

      END SUBROUTINE SurfAverage
     
      END MODULE quantities
          
