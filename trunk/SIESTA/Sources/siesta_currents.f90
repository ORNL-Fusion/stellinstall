      MODULE siesta_currents
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE descriptor_mod, ONLY: SIESTA_COMM
      USE shared_data, ONLY: l_update_state, l_natural, lasym,          &
                             jsju_ratio_s, jsju_ratio_a, ste 
      USE fourier
      USE nscalingtools, ONLY: nranks, startglobrow, endglobrow, MPI_ERR
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE metrics, ONLY: tolowerf, tolowerh
      USE quantities

      IMPLICIT NONE
      PRIVATE
      INTEGER  :: nsmin, nsmax, n1, n2
      INTEGER, PARAMETER :: isym=0, iasym=1, m0=0, m1=1, m2=2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
                bsubsijh, bsubuijh, bsubvijh, bs_filter, bu_filter,     &
                bv_filter
      REAL(dp) :: edge_cur
      PUBLIC   :: cv_currents

      CONTAINS

      SUBROUTINE cv_currents (bsupsijh, bsupuijh, bsupvijh,             &
                              ksupsijf, ksupuijf, ksupvijf,             & 
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
!>  \param[out] ksupsijf K^s on the full grid.
!>  \param[out] ksupuijf K^u on the full grid.
!>  \param[out] ksupvijf K^v on the full grid.
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
      USE descriptor_mod, ONLY: iam
      USE timer_mod
      USE mpi_inc
      USE descriptor_mod, ONLY: SIESTA_COMM
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iwrite
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)  ::           &
           bsupsijh, bsupuijh, bsupvijh,                                &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(IN)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsq_loc
      REAL(dp) :: beta, ton, toff, stes(6), rtes(6), tmp
!-----------------------------------------------
      CALL second0(ton)

      n1=MAX(1, startglobrow); n2=MIN(ns,endglobrow)
      nsmin = LBOUND(bsupsijh,3); nsmax = UBOUND(bsupsijh,3)

!      CALL ASSERT(LBOUND(ksupsijf,3).EQ.nsmin,                          &
!                 'KSUPSIJF HAS WRONG LB IN CV_CURRENTS')
!      CALL ASSERT(UBOUND(ksupsijf,3).EQ.nsmax,                          &
!                 'KSUPSIJF HAS WRONG UB IN CV_CURRENTS')

! Calculate covariant BsubXh on half mesh
      ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                      &
               bsubuijh(ntheta,nzeta,nsmin:nsmax),                      &
               bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.EQ.0, 'Allocation1 failed in CV_CURRENTS')
     
      CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                       &
                    bsubsijh, bsubuijh, bsubvijh, nsmin, nsmax)
                    
!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
      IF (.NOT. ALLOCATED(mupar_norm)) THEN
         ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
         CALL ASSERT(istat.EQ.0,'MUPAR_NORM ALLOCATION FAILED!')
         beta = MAX(1.E-5_dp, wp_i/wb_i)
         mupar_norm = signjac*beta / rmajor_i**2  !/ vp_f(nsmin:nsmax)     !add 1/jac correction factor in AddDamping
         IF (iam.EQ.0 .AND. lverbose)                                   &
            PRINT '(a,1p,e12.3)',' Rmajor_i = ',rmajor_i
      END IF

      IF (lmagen) THEN
         ALLOCATE (bsq_loc(ntheta,nzeta,nsmin:nsmax), stat=istat)
         bsq_loc = bsupsijh*bsubsijh  + bsupuijh*bsubuijh               &
                 + bsupvijh*bsubvijh
         IF (n1 .EQ. 1) bsq_loc(:,:,1) = 0
         CALL ASSERT(n1.EQ.1 .OR. ALL(bsq_loc(:,:,n1:n2) .GT. zero),    &
                     'BSQ <= 0 in cv_currents_par!')

         wb = SUM(bsq_loc(:,:,n1:n2)*jacobh(:,:,n1:n2)*wint(:,:,n1:n2))
         wb0 = signjac*(twopi*pi*hs_i)*wb
         
         IF (PARSOLVER) THEN
#if defined(MPI_OPT)
            CALL MPI_ALLREDUCE(wb0, wb, 1, MPI_REAL8, MPI_SUM,          &
                               SIESTA_COMM, MPI_ERR)
#endif
         ELSE
            wb = wb0
         END IF
         DEALLOCATE (bsq_loc)
      END IF

!
! Calculate BsubXmnh and KsupXmnf harmonics
!       
      CALL CURLB (ksupsmnsf, ksupumncf, ksupvmncf, jsju_ratio_s,        &
                  isym, lcurr)
      IF (lasym)                                                        &
      CALL CURLB (ksupsmncf, ksupumnsf, ksupvmnsf, jsju_ratio_a,        &
                  iasym, lcurr)

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF

!
!     Calculate full mesh, contravariant real-space current components
!     KsupXF = sqrt(g)*JsupXF 
!        
      CALL GETCV (ksupsmnsf, ksupumncf, ksupvmncf, ksupsijf, ksupuijf,  &
                  ksupvijf, isym)                       
      IF (lasym)                                                        &
      CALL GETCV (ksupsmncf, ksupumnsf, ksupvmnsf, ksupsijf, ksupuijf,  &
                  ksupvijf, iasym)                       

!Covariant K components (with jacob factor still there) in real space
!==> needed for current diffusion
      IF (lcurr) THEN
         CALL tolowerf(ksupsijf, ksupuijf, ksupvijf,                    & 
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
            tmp = SUM(wint(:,:,ns)*(bsubvijh(:,:,ns)**2                  &
                +                   bsubuijh(:,:,ns)**2))
            IF (tmp .NE. zero) tmp = SQRT(edge_cur/tmp)
         END IF
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
         CALL MPI_BCAST(tmp,1,MPI_REAL8,nranks-1,SIESTA_COMM,MPI_ERR)
         END IF
#endif
         IF (iam .EQ. 0) THEN
         DO iwrite = 6, unit_out, unit_out-6
            IF (.NOT.lverbose .AND. iwrite.EQ.6) CYCLE
            WRITE(iwrite, 1008)' RMS EDGE CURRENT ||KSUPS(NS)|| = ', tmp
         END DO
         CALL FLUSH(unit_out)
         END IF
 1008 FORMAT(a,1p,e10.3)

         stes(1) = SUM(wint(:,:,n1:n2)*(bs_filter - bsubsijh(:,:,n1:n2))**2)
         stes(2) = SUM(wint(:,:,n1:n2)*(bs_filter + bsubsijh(:,:,n1:n2))**2)         
         stes(3) = SUM(wint(:,:,n1:n2)*(bu_filter - bsubuijh(:,:,n1:n2))**2)
         stes(4) = SUM(wint(:,:,n1:n2)*(bu_filter + bsubuijh(:,:,n1:n2))**2)         
         stes(5) = SUM(wint(:,:,n1:n2)*(bv_filter - bsubvijh(:,:,n1:n2))**2)
         stes(6) = SUM(wint(:,:,n1:n2)*(bv_filter + bsubvijh(:,:,n1:n2))**2)         
         IF (PARSOLVER) THEN
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(stes, rtes, 6, MPI_REAL8, MPI_SUM,          &
                            SIESTA_COMM, MPI_ERR)
#endif
         ELSE
         rtes = stes
         END IF
         IF (rtes(2) .NE. zero) ste(2) = SQRT(rtes(1)/rtes(2))
         IF (rtes(4) .NE. zero) ste(3) = SQRT(rtes(3)/rtes(4))
         IF (rtes(6) .NE. zero) ste(4) = SQRT(rtes(5)/rtes(6))
         DEALLOCATE (bs_filter, bu_filter, bv_filter)
      END IF DIAGNO

      DEALLOCATE (bsubsijh, bsubuijh, bsubvijh) 

      CALL second0(toff)
#if defined(SKS)      
      cv_current_time=cv_current_time+(toff-ton)
#endif      
      time_current = time_current+(toff-ton)

      END SUBROUTINE cv_currents

  
!>     \brief Compute fourier components of contravariant currents on 
!>      full radial grid
      SUBROUTINE CURLB (ksupsmnf, ksupumnf, ksupvmnf, jsju_ratio,       &
                        ipar, lcurr)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE   ::                      &
           ksupsmnf, ksupumnf, ksupvmnf
      REAL(dp)             :: jsju_ratio
      LOGICAL, INTENT(IN)  :: lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, fcomb, nmin, nmax,               &
                  m, n, mp, np, sparity, istat, ns_span
!----------------------------------------------- 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         bsubsmnf, bsubumnf, bsubvmnf, bsubsmnh, bsubumnh, bsubvmnh,    &
         db_uds, db_vds
      REAL(dp) :: r12, temp(-ntor:ntor)
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
         sparity=1; edge_cur = 0
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
         sparity=-1
      END IF
      r12 = sparity*hs_i/2

      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      CALL ASSERT(istat.EQ.0,'Allocation1 failed in CURLB')

      CALL tomnsp_par(bsubsijh, bsubsmnh, fours)                          ! Compute BSUBXMNH    
      CALL tomnsp_par(bsubuijh, bsubumnh, fouru)
      CALL tomnsp_par(bsubvijh, bsubvmnh, fourv)

      ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                 & ! Compute BSUBXMNF
               bsubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      CALL ASSERT(istat.EQ.0,'Allocation2 failed in CURLB')
 
      CALL to_full_mesh(bsubsmnh, bsubsmnf)                       
      CALL to_full_mesh(bsubumnh, bsubumnf)
      CALL to_full_mesh(bsubvmnh, bsubvmnf)

!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               db_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation3 failed in CURLB')

!     Radial gradients of bsubmnu,v on full grid (excluding origin and edge: treat separately)
      CALL GradientFull(db_uds, bsubumnh)
      CALL GradientFull(db_vds, bsubvmnh)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
!     NOTE: the bc B^s = 0 at js=ns means there are no contributions from 
!           these derivative terms there
      IF (nsmin .EQ. 1) THEN
         bsubsmnf(m2:,:,1) = 0
         bsubumnf(:,:,1) = 0;  bsubumnf(m1,:,1) = bsubumnh(m1,:,2)
!     to satisfy k^s = r12*k^u at s=r12, for m=1, must have b_u = r12*b_s for m=1
!     this is not completely consistent with the delta-W jogs in CheckForces!
         IF (.NOT.l_natural) THEN
            temp = (r12*bsubsmnf(m1,:,1) - bsubumnf(m1,:,1))/2
            bsubsmnf(m1,:,1) = bsubsmnf(m1,:,1) - temp/r12
            bsubumnf(m1,:,1) = bsubumnf(m1,:,1) + temp
         END IF

         db_uds(:,:,1) = 0;  db_uds(m0:m1,:,1) = 2*ohs*bsubumnh(m0:m1,:,2)
         db_vds(:,:,1) = 0;  db_vds(m1,:,1) = 2*ohs*bsubvmnh(m1,:,2)
         bsubvmnf(m1+1:,:,1) = 0  
         CALL ASSERT(ALL(bsubvmnf(m1,:,1).EQ.bsubvmnh(m1,:,2)),         &
                     'b_vmnf(1) ! b_vmnh(2) in CURLB')
      END IF

!     NOTE: ksupXmnf(ns) is approximately ksupXmnh(ns) a half-grid point in from edge

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
           ksupsmnf(m,n,nmin:nmax) = -mp*bsubvmnf(m,n,nmin:nmax)        &
                                   +  np*bsubumnf(m,n,nmin:nmax)
           ksupumnf(m,n,nmin:nmax) =  np*bsubsmnf(m,n,nmin:nmax)        &
                                   -  db_vds(m,n,nmin:nmax)                     
           ksupvmnf(m,n,nmin:nmax) = -mp*bsubsmnf(m,n,nmin:nmax)        &
                                   +  db_uds(m,n,nmin:nmax)
         ENDDO
      ENDDO

!Diagnostics
      DIAGNO: IF (l_update_state) THEN
         IF (nsmax .EQ. ns) edge_cur = edge_cur                         & 
                                     + SUM(ksupsmnf(:,:,ns)**2)

!Deallocated in calling routine
         ns_span = n2-n1+1
         IF (.NOT. ALLOCATED(bs_filter))                                &
         ALLOCATE(bs_filter(ntheta,nzeta,ns_span),                      &
                  bu_filter(ntheta,nzeta,ns_span),                      &
                  bv_filter(ntheta,nzeta,ns_span), stat=istat)
         CALL ASSERT(istat.EQ.0,'Allocation4 failed in CURLB')
         CALL toijsp_par(bsubsmnh(:,:,n1:n2), bs_filter, fcomb, fours)
         CALL toijsp_par(bsubumnh(:,:,n1:n2), bu_filter, fcomb, fouru)
         CALL toijsp_par(bsubvmnh(:,:,n1:n2), bv_filter, fcomb, fourv)
	     IF (n1 .EQ. 1) bv_filter(:,:,1) = bv_filter(:,:,2)
	     IF (n1 .EQ. 1)	bsubvijh(:,:,1) = bsubvijh(:,:,2)

!K^s-r12*K^u at origin BEFORE bc applied
         IF (nsmin .EQ. 1) THEN
            temp = ksupsmnf(m1,:,1) + r12*ksupumnf(m1,:,1)
            jsju_ratio = SUM(temp**2)
            IF (jsju_ratio .GT. zero) THEN
               temp = ksupsmnf(m1,:,1) - r12*ksupumnf(m1,:,1)
               jsju_ratio = SUM(temp**2) / jsju_ratio
               jsju_ratio = SQRT(jsju_ratio)
            END IF
         END IF

      END IF DIAGNO

      DEALLOCATE (bsubsmnh, bsubumnh, bsubvmnh,                         &
                  bsubsmnf, bsubumnf, bsubvmnf, stat=istat)   
      CALL ASSERT(istat.EQ.0,'Deallocation failed in CURLB')

      IF (nsmin .EQ. 1) THEN
          ksupsmnf(m0,:,1) = 0; ksupsmnf(m2:,:,1) = 0
          ksupumnf(m2:,:,1) = 0
          ksupvmnf(m1:,:,1) = 0
      END IF

      END SUBROUTINE CURLB


!>     \brief Compute full-grid contravariant (sup) currents in real space
      SUBROUTINE GETCV (ksupsmnf, ksupumnf, ksupvmnf,                   &
                        ksupsijf, ksupuijf, ksupvijf, ipar)                       
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: ipar
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                        &
           ksupsmnf, ksupumnf, ksupvmnf,                                &
           ksupsijf, ksupuijf, ksupvijf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, fcomb
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
      END IF

      CALL toijsp_par(ksupsmnf(:,:,nsmin:nsmax), ksupsijf, fcomb, fours)
      CALL toijsp_par(ksupumnf(:,:,nsmin:nsmax), ksupuijf, fcomb, fouru)
      CALL toijsp_par(ksupvmnf(:,:,nsmin:nsmax), ksupvijf, fcomb, fourv)

      END SUBROUTINE GETCV


      END MODULE siesta_currents
