      MODULE siesta_currents
      USE stel_kinds
      USE stel_constants
      USE shared_data, ONLY: l_update_state, l_natural, lasym,           &
                             jsju_ratio_s, jsju_ratio_a, ste 
      USE fourier 
      IMPLICIT NONE
      PRIVATE
      INTEGER  :: nsmin, nsmax, ns_span, n1, n2
      INTEGER, PARAMETER :: isym=0, iasym=1, m0=0, m1=1, m2=2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                        &
                bsubsijh, bsubuijh, bsubvijh, bs_filter, bu_filter,      &
                bv_filter
      REAL(dp) :: edge_cur
      PUBLIC   :: cv_currents

      CONTAINS

      SUBROUTINE cv_currents (bsupsijh, bsupuijh, bsupvijh,              &
                              ksupsijf, ksupuijf, ksupvijf,              & 
                              lmagen, lcurr)
!*******************************************************************************
!>  \file siesta_currents.f90
!>  \brief Contains routines to compute the current density from Curl(B).
!>  Takes the Curl(B) to return the contravariant current density. Contravariant
!>  components of the magnetic field must be converted to covariant components
!>  first.
!>
!>  \param[in]  bsupsijh B^s on the half grid.
!>  \param[in]  bsupuijh B^u on the half grid.
!>  \param[in]  bsupvijh B^v on the half grid.
!>  \param[out] ksupsijf K^s on the full rgid.
!>  \param[out] ksupuijf K^u on the full rgid.
!>  \param[out] ksupvijf K^v on the full rgid.
!>  \param[in]  lmagen   logical controlling calc of magnetic energy wb
!>  \param[in]  lcurr    UNKNOWN

!
!  Note separating the Doxygen comment block here so detailed description is
!  found in the Module not the file.
!
!>  \brief Computes contravariant currents on full-grid from covariant B-fields
!>  on half-grid.
!>  \author S. P. Hirshman
!>  \date February 26 2007
!>
!>  Current density K is defined as
!>
!>     K = Curl(B)/mu0                                                       (1)
!>
!>  In curvilinear coordinates this becomes
!>
!>     K^s = 1/sqrt(g)(dB_v/du - dB_u/dv)                                   (2a)
!>
!>     K^u = 1/sqrt(g)(dB_s/dv - dB_v/ds)                                   (2b)
!>
!>     K^v = 1/sqrt(g)(dB_u/ds - dB_s/du)                                   (2c)
!>
!>  where sqrt(g) is the coordinate system Jacobian. Since there are derivatives 
!>  with respect to s, these routines take the half grid values as input. The
!>  contravariant current computed will be on the full grid. These routines also
!>  output currents in real space. Note the computed Kmn are [sqrt(g)*K]mn with the 
!>  embedded jacobian factor.
!*******************************************************************************
      USE island_params, ONLY: nu_i, rmajor_i
      USE Hessian, ONLY: mupar_norm
      USE nscalingtools, ONLY: nranks, startglobrow, endglobrow, MPI_ERR
      USE descriptor_mod, ONLY: iam
      USE quantities
      USE timer_mod
      USE mpi_inc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iwrite
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)  ::            &  !FOR NOW, ASSUME THESE ARE ARRAYS FOR ALL NS!
           bsupsijh, bsupuijh, bsupvijh,                                 &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(IN)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
         bsubsmnsh, bsubumnch, bsubvmnch,                                &
         bsubsmnsf, bsubumncf, bsubvmncf,                                &
         bsubsmnch, bsubumnsh, bsubvmnsh,                                &
         bsubsmncf, bsubumnsf, bsubvmnsf,                                &
         bsq_loc
      REAL(dp) :: beta, r0, ton, toff, stes(6), rtes(6), wb0
!-----------------------------------------------
      CALL second0(ton)

      n1=MAX(1, startglobrow); n2=MIN(ns,endglobrow)

      IF (lcurr) THEN
         nsmin=MAX(1, startglobrow-1); nsmax=MIN(ns,endglobrow+2)
      ELSE
         nsmin=MAX(1, startglobrow);   nsmax=MIN(ns,endglobrow+1)
      END IF

      ns_span = nsmax-nsmin+1

      IF (LBOUND(bsupsijh,3) .NE. nsmin) STOP 'BSUPSIJH HAS WRONG LB'
      IF (UBOUND(bsupsijh,3) .NE. nsmax) STOP 'BSUPSIJH HAS WRONG UB'
      IF (LBOUND(ksupsijf,3) .NE. nsmin) STOP 'KSUPSIJF HAS WRONG LB'
      IF (UBOUND(ksupsijf,3) .NE. nsmax) STOP 'KSUPSIJF HAS WRONG UB'

! Calculate covariant BsubXh on half mesh 
!       
      ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                       &
               bsubuijh(ntheta,nzeta,nsmin:nsmax),                       &
               bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh_par(bsupsijh, bsupuijh, bsupvijh,                    &
                        bsubsijh, bsubuijh, bsubvijh,                    &
                        nsmin, nsmax)
!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
      IF (.NOT. ALLOCATED(mupar_norm)) THEN
         ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
         IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
         beta = MAX(1.E-5_dp, wp_i/wb_i)
         mupar_norm = signjac*beta / rmajor_i**2  !/ vp_f(nsmin:nsmax)     !put 1/jac correction factor in AddDamping
         IF (iam .EQ. 0) PRINT '(a,1p,e12.3)',' Rmajor_i = ',rmajor_i
      END IF

      IF (lmagen) THEN
         ALLOCATE (bsq_loc(ntheta,nzeta,nsmin:nsmax), stat=istat)
         bsq_loc = bsupsijh*bsubsijh  + bsupuijh*bsubuijh                &
                 + bsupvijh*bsubvijh
         IF (n1 .EQ. 1) bsq_loc(:,:,1) = 0
         IF (n1 .GT. 1) THEN
            IF(ANY(bsq_loc(:,:,n1:n2) .LE. zero)) THEN
               IF (iam .EQ. 0) PRINT *,'BSQ <= 0 in cv_currents_par!'
#if defined(MPI_OPT)
               CALL MPI_ABORT(MPI_COMM_WORLD, 2, MPI_ERR)
#endif
            END IF
         END IF

         wb = SUM(bsq_loc(:,:,n1:n2)*jacobh(:,:,n1:n2)*wint(:,:,n1:n2))
         wb0 = signjac*(twopi*pi*hs_i)*wb
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(wb0,wb,1,MPI_REAL8, MPI_SUM,                 &
                            MPI_COMM_WORLD,MPI_ERR)
#else
         wb = wb0
#endif
         DEALLOCATE (bsq_loc)
      END IF

!
! Calculate BsubXmnh and KsupXmnf harmonics
!       
      CALL CURLB (ksupsmnsf, ksupumncf, ksupvmncf, jsju_ratio_s,         &
                  isym, lcurr)
      IF (lasym)                                                         &
      CALL CURLB (ksupsmncf, ksupumnsf, ksupvmnsf, jsju_ratio_a,         &
                  iasym, lcurr)

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF

!
!     Calculate full mesh, contravariant current real-space current components
!     KsupXF = sqrt(g)*JsupXF 
!        
      ns_span = nsmax-nsmin+1
      CALL GETCV (ksupsmnsf, ksupumncf, ksupvmncf, ksupsijf, ksupuijf,   &
                  ksupvijf, isym)                       
      IF (lasym)                                                         &
      CALL GETCV (ksupsmncf, ksupumnsf, ksupvmnsf, ksupsijf, ksupuijf,   &
                  ksupvijf, iasym)                       

!Co-variant K in real space needed for current diffusion
      IF (lcurr) THEN
         CALL tolowerf_par(ksupsijf, ksupuijf, ksupvijf,                 & 
                           ksubsijf, ksubuijf, ksubvijf, nsmin, nsmax)
         IF (nsmin .EQ. 1) ksubuijf(:,:,1) = 0 
      ENDIF

!Diagnostic output (only compute when state is updated)
      DIAGNO: IF (l_update_state) THEN
!
!     Compute spectral truncation error measure
!     (Use ksupsijf for temporary FILTERED bsubXijh values)
!
         IF (nsmax .EQ. ns) THEN 
            r0 = SUM(wint(:,:,ns)*(bsubvijh(:,:,ns)**2                   &
               +                   bsubuijh(:,:,ns)**2))
            r0 = SQRT(edge_cur/r0)
         END IF
#if defined(MPI_OPT)
         CALL MPI_BCAST(r0,1,MPI_REAL8,nranks-1,MPI_COMM_WORLD,MPI_ERR)
#endif
         IF (iam .EQ. 0) THEN
         DO iwrite = 6, 33, 27
            WRITE(iwrite, 1008)' RMS EDGE CURRENT ||KSUPS(NS)|| = ', r0
         END DO
         END IF
 1008 FORMAT(a,1p,e10.3)

         ns_span = n2-n1+1
         stes(1) = SUM(wint(:,:,n1:n2)*(bs_filter - bsubsijh(:,:,n1:n2))**2)
         stes(2) = SUM(wint(:,:,n1:n2)*(bs_filter + bsubsijh(:,:,n1:n2))**2)         
         stes(3) = SUM(wint(:,:,n1:n2)*(bu_filter - bsubuijh(:,:,n1:n2))**2)
         stes(4) = SUM(wint(:,:,n1:n2)*(bu_filter + bsubuijh(:,:,n1:n2))**2)         
         stes(5) = SUM(wint(:,:,n1:n2)*(bv_filter - bsubvijh(:,:,n1:n2))**2)
         stes(6) = SUM(wint(:,:,n1:n2)*(bv_filter + bsubvijh(:,:,n1:n2))**2)         
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(stes,rtes,6,MPI_REAL8,MPI_SUM,               &
                            MPI_COMM_WORLD,MPI_ERR)
#else
         rtes = stes
#endif
         IF (rtes(2) .NE. zero) ste(2) = SQRT(rtes(1)/rtes(2))
         IF (rtes(4) .NE. zero) ste(3) = SQRT(rtes(3)/rtes(4))
         IF (rtes(6) .NE. zero) ste(4) = SQRT(rtes(5)/rtes(6))
         DEALLOCATE (bs_filter, bu_filter, bv_filter)
      END IF DIAGNO

      DEALLOCATE (bsubsijh, bsubuijh, bsubvijh) 


      CALL second0(toff)
      cv_current_time=cv_current_time+(toff-ton)
      time_current = time_current+(toff-ton)

      END SUBROUTINE cv_currents

  
      SUBROUTINE CURLB (ksupsmnf, ksupumnf, ksupvmnf, jsju_ratio,        &
                        iparity, lcurr)
      
      USE nscalingtools, ONLY: endglobrow
      USE quantities, ONLY: ohs
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: iparity
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE   ::                       &
           ksupsmnf, ksupumnf, ksupvmnf
      REAL(dp)             :: jsju_ratio
      LOGICAL, INTENT(IN)  :: lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, fcomb, nblock, nmin, nmax,        &
                  m, n, mp, np, sparity, istat, ns_span1
!----------------------------------------------- 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
         bsubsmnf, bsubumnf, bsubvmnf, bsubsmnh, bsubumnh, bsubvmnh,     &
         db_uds, db_vds
      REAL(dp) :: r0, temp(-ntor:ntor)
!-----------------------------------------------
      IF (iparity .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
         r0 = hs_i/2; sparity=1; edge_cur = 0
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
         r0 = -hs_i/2; sparity=-1
      END IF

      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                  & ! Get BSUBMN on half mesh
               bsubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
               bsubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .NE. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp_par(bsubsijh, bsubsmnh, fours, 1, ns_span)                ! Calculate BSUBXMNH
      CALL tomnsp_par(bsubuijh, bsubumnh, fouru, 1, ns_span)
      CALL tomnsp_par(bsubvijh, bsubvmnh, fourv, 1, ns_span)

      ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                  &  ! Next, get BSUBMN on full mesh
               bsubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
               bsubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .NE. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
      nblock=SIZE(bsubsmnh,1)*SIZE(bsubsmnh,2)
      CALL to_full_mesh_par(bsubsmnh, bsubsmnf, nblock, nsmin, nsmax)                       
      CALL to_full_mesh_par(bsubumnh, bsubumnf, nblock, nsmin, nsmax)
      CALL to_full_mesh_par(bsubvmnh, bsubvmnf, nblock, nsmin, nsmax)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      IF (nsmin .EQ. 1) THEN
         bsubsmnf(:,:,1) = 0
         bsubumnf(:,:,1) = 0
         bsubsmnf(m1,:,1) = bsubsmnh(m1,:,1)
         bsubumnf(m1,:,1) = bsubumnh(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
         bsubvmnf(:,:,1) = 0
         bsubvmnf(m0:m1,:,1) = bsubvmnh(m0:m1,:,2)
      END IF
!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                    &
               db_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation4 failed in CV_CURRENTS'

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull_par(db_uds, bsubumnh, nblock, nsmin, nsmax)
      CALL GradientFull_par(db_vds, bsubvmnh, nblock, nsmin, nsmax)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      IF (nsmin .EQ. 1) THEN
         db_uds(m0:m1,:,1) = 2*ohs*bsubumnh(m0:m1,:,2)
         bsubvmnf(m1,:,1) = bsubvmnh(m1,:,2)
         db_vds(m1,:,1)    = 2*ohs*bsubvmnh(m1,:,2)                      !--> 0
         bsubumnf(m1,:,1) = r0*bsubsmnf(m1,:,1)                          !--> 0
      END IF

!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1

      nmin = nsmin
      IF (lcurr) THEN
         nmax=MIN(ns,endglobrow+1)
      ELSE
         nmax=MIN(ns,endglobrow)
      END IF
      DO m = 0, mpol
         mp = sparity*m
         DO n = -ntor, ntor
           np = sparity*n*nfp
           ksupsmnf(m,n,nmin:nmax) = -mp*bsubvmnf(m,n,nmin:nmax)         &
                                   +  np*bsubumnf(m,n,nmin:nmax)
           ksupumnf(m,n,nmin:nmax) =  np*bsubsmnf(m,n,nmin:nmax)         &
                                   -  db_vds(m,n,nmin:nmax)                     
           ksupvmnf(m,n,nmin:nmax) = -mp*bsubsmnf(m,n,nmin:nmax)         &
                                   +  db_uds(m,n,nmin:nmax)
         ENDDO
      ENDDO

!Diagnostics
      DIAGNO: IF (l_update_state) THEN
         IF (nsmax .EQ. ns) edge_cur = edge_cur                          & 
                                     + SUM(ksupsmnf(:,:,ns)**2)

!Deallocated in calling routine
         ns_span1 = n2-n1+1
         IF (.NOT. ALLOCATED(bs_filter))                                 &
         ALLOCATE(bs_filter(ntheta,nzeta,ns_span1),                      &
                  bu_filter(ntheta,nzeta,ns_span1),                      &
                  bv_filter(ntheta,nzeta,ns_span1), stat=istat)
         IF (istat .NE. 0) STOP 'Allocation5 failed in CV_CURRENTS'
         CALL toijsp_par(bsubsmnh(:,:,n1:n2), bs_filter, fcomb, fours, 1, ns_span1)
         CALL toijsp_par(bsubumnh(:,:,n1:n2), bu_filter, fcomb, fouru, 1, ns_span1)
         CALL toijsp_par(bsubvmnh(:,:,n1:n2), bv_filter, fcomb, fourv, 1, ns_span1)
	 IF (n1 .EQ. 1) bv_filter(:,:,1) = bv_filter(:,:,2)
	 IF (n1 .EQ. 1)	bsubvijh(:,:,1) = bsubvijh(:,:,2)

!K^s-r0*K^u at origin BEFORE bc applied
         IF (nsmin .EQ. 1) THEN
            temp = ksupsmnf(m1,:,1) + r0*ksupumnf(m1,:,1)
            jsju_ratio = SUM(temp**2)
            IF (jsju_ratio .GT. zero) THEN
               temp = ksupsmnf(m1,:,1) - r0*ksupumnf(m1,:,1)
               jsju_ratio = SUM(temp**2) / jsju_ratio
               jsju_ratio = SQRT(jsju_ratio)
            END IF
         END IF

      END IF DIAGNO

      DEALLOCATE (bsubsmnh, bsubumnh, bsubvmnh)
      DEALLOCATE (bsubsmnf, bsubumnf, bsubvmnf)   

!Apply exact bc K^s = r0*K^u for m=1 at origin
      IF (nsmin .EQ. 1) THEN
         IF (l_natural) ksupsmnf(m1,:,1) = r0*ksupumnf(m1,:,1)
         ksupsmnf(m0,:,1) = 0;  ksupsmnf(m2:,:,1) = 0
         ksupumnf(m0,:,1) = 0;  ksupumnf(m2:,:,1) = 0
         ksupvmnf(m1:,:,1)= 0
      END IF

      END SUBROUTINE CURLB


      SUBROUTINE GETCV (ksupsmnf, ksupumnf, ksupvmnf,                    &
                        ksupsijf, ksupuijf, ksupvijf, iparity)                       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: iparity
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                         &
           ksupsmnf, ksupumnf, ksupvmnf,                                 &
           ksupsijf, ksupuijf, ksupvijf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, fcomb
!-----------------------------------------------
      IF (iparity .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
      END IF

      CALL toijsp_par(ksupsmnf(:,:,nsmin:nsmax), ksupsijf, fcomb, fours, &
                      1, ns_span)
      CALL toijsp_par(ksupumnf(:,:,nsmin:nsmax), ksupuijf, fcomb, fouru, &
                      1, ns_span)
      CALL toijsp_par(ksupvmnf(:,:,nsmin:nsmax), ksupvijf, fcomb, fourv, &
                      1, ns_span)

      END SUBROUTINE GETCV


      END MODULE siesta_currents
