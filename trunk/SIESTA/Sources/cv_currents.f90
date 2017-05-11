!#if defined(SKS)
      SUBROUTINE cv_currents_par (bsupsijh, bsupuijh, bsupvijh,          &
                                  ksupsijf, ksupuijf, ksupvijf,          &
                                  lmagen, lcurr)
!*******************************************************************************
!>  \file cv_currents.f90
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
!  Note separating the Doxygen comment block here so detailed decription is
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
!>  In curvalinear coordinates this becomes
!>
!>     K^s = 1/J(dB_v/du - dB_u/dv)                                         (2a)
!>
!>     K^u = 1/J(dB_s/dv - dB_v/ds)                                         (3b)
!>
!>     K^v = 1/J(dB_u/ds - dB_s/du)                                         (4c)
!>
!>  where J is the coordinate system Jacobian. Since there are derivatives with
!>  respect to s, these routines take the half grid values as input. The
!>  contravariant current computed will be on the full grid. These routines also
!>  operate in real space.
!*******************************************************************************
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i, rmajor_i
      USE evolution, ONLY: l_natural, jsju_ratio_s, jsju_ratio_a, ste 
      USE shared_data, ONLY: l_update_state
      USE timer_mod
      USE Hessian, ONLY: mupar_norm
      USE nscalingtools, ONLY: nranks, startglobrow, endglobrow, MPI_ERR
      USE descriptor_mod, ONLY: iam
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iwrite
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(IN)  ::               &  !FOR NOW, ASSUME THESE ARE ARRAYS FOR ALL NS!
           bsupsijh, bsupuijh, bsupvijh
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(OUT) ::               &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(IN)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, m2=2
      INTEGER :: m, n, istat, l, js, nblock
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
         bsubsijh, bsubuijh, bsubvijh,                                   & 
         bsubsmnsh, bsubumnch, bsubvmnch,                                &
         bsubsmnsf, bsubumncf, bsubvmncf,                                &
         bsubsmnch, bsubumnsh, bsubvmnsh,                                &
         bsubsmncf, bsubumnsf, bsubvmnsf,                                &
         db_udss, db_vdss, db_udsc, db_vdsc, bsq_loc
      REAL(dp) :: beta, r0, temp(-ntor:ntor), skston, skstoff,           &
                     ton, toff, stes(6), rtes(6), wb0
      INTEGER     :: k, nsmin, nsmax, ns_span, n1, n2
!-----------------------------------------------
      CALL second0(ton)
      skston=ton

      n1=MAX(1, startglobrow); n2=MIN(ns,endglobrow)

      IF (lcurr) THEN
         nsmin=MAX(1, startglobrow-1); nsmax=MIN(ns,endglobrow+2)
      ELSE
         nsmin=MAX(1, startglobrow);   nsmax=MIN(ns,endglobrow+1)
      END IF

      ns_span = nsmax-nsmin+1
!
! Calculate covariant BsubXh on half mesh 
!       
      ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                       &
               bsubuijh(ntheta,nzeta,nsmin:nsmax),                       &
               bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh_par(bsupsijh, bsupuijh, bsupvijh,                    &
                        bsubsijh, bsubuijh, bsubvijh,                    &
                        nsmin, nsmax)

      IF (lmagen) THEN
         ALLOCATE (bsq_loc(ntheta,nzeta,nsmin:nsmax), stat=istat)
         bsq_loc = bsupsijh(:,:,1:ns_span)*bsubsijh                      & 
                 + bsupuijh(:,:,1:ns_span)*bsubuijh                      &
                 + bsupvijh(:,:,1:ns_span)*bsubvijh
         IF (n1 .EQ. 1) bsq_loc(:,:,1) = 0
         IF (n1 .GT. 1) THEN
            IF(ANY(bsq_loc(:,:,n1:n2) .LE. zero)) THEN
               IF (iam .EQ. 0) PRINT *,'BSQ <= 0 in cv_currents_par!'
#if defined(MPI_OPT)
               CALL MPI_ABORT(MPI_COMM_WORLD, 2, MPI_ERR)
#endif
            END IF
         END IF

         wb = 0
         DO l = 1, nu_i
            wb = wb + SUM(bsq_loc(l,:,n1:n2)*jacobh(l,:,n1:n2))*cosmui(l,m0)
         END DO         
         wb0 = signjac*(twopi*twopi*hs_i)*wb/2                              !2X : theta 0-pi
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(wb0,wb,1,MPI_REAL8, MPI_SUM,                 &
                            MPI_COMM_WORLD,MPI_ERR)
#else
         wb = wb0
#endif
!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
         IF (.NOT. ALLOCATED(mupar_norm)) THEN
            ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
            IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
            beta = MAX(1.E-5_dp, wp_i/wb_i)
            mupar_norm = signjac*beta / rmajor_i**2  !/ vp_f(nsmin:nsmax)               !!put 1/jac correction factor in AddDamping
            IF (iam .EQ. 0) PRINT '(a,1p,e12.3)',' Rmajor_i = ',rmajor_i
         END IF
         DEALLOCATE (bsq_loc)
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!       
      ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                 & ! First, get BSUBMN on half mesh
               bsubumnch(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .NE. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp_par(bsubsijh, bsubsmnsh, f_sin, f_full, 1, ns_span)      ! Calculate BSUB?MN?H
      CALL tomnsp_par(bsubuijh, bsubumnch, f_cos, f_full, 1, ns_span)
      CALL tomnsp_par(bsubvijh, bsubvmnch, f_cos, f_full, 1, ns_span)
      CALL second0(toff)

!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF 
!        
                    
      ALLOCATE(bsubsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &  ! Next, get BSUBMN on full mesh
               bsubumncf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .NE. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
      nblock=SIZE(bsubsmnsh,1)*SIZE(bsubsmnsh,2)
      CALL to_full_mesh_par(bsubsmnsh, bsubsmnsf, nblock, nsmin, nsmax)                       
      CALL to_full_mesh_par(bsubumnch, bsubumncf, nblock, nsmin, nsmax)
      CALL to_full_mesh_par(bsubvmnch, bsubvmncf, nblock, nsmin, nsmax)
      CALL second0(toff)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      IF (nsmin .EQ. 1) THEN
         bsubsmnsf(:,:,1) = 0
         bsubumncf(:,:,1) = 0
         bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,1)
         bsubumncf(m1,:,1) = bsubumnch(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
         bsubvmncf(:,:,1) = 0
         bsubvmncf(m0:m1,:,1) = bsubvmnch(m0:m1,:,2)
      END IF
!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(db_udsc(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               db_vdsc(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation4 failed in CV_CURRENTS'
      CALL second0(toff)

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull_par(db_udsc, bsubumnch, nblock, nsmin, nsmax)
      CALL GradientFull_par(db_vdsc, bsubvmnch, nblock, nsmin, nsmax)

      CALL second0(toff)
      r0 = hs_i/2

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      IF (nsmin .EQ. 1) THEN
         db_udsc(m0:m1,:,1) = 2*ohs*bsubumnch(m0:m1,:,2)
         bsubvmncf(m1,:,1) = bsubvmnch(m1,:,2)
         db_vdsc(m1,:,1)    = 2*ohs*bsubvmnch(m1,:,2)                      !--> 0
         bsubumncf(m1,:,1) = r0*bsubsmnsf(m1,:,1)                          !--> 0
      END IF


!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
      CALL second0(toff)

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF
      DO n = -ntor, ntor
         DO m = 0, mpol
           ksupsmnsf(m,n,nsmin:nsmax) = -m*bsubvmncf(m,n,nsmin:nsmax) + n*nfp*bsubumncf(m,n,nsmin:nsmax)
           ksupumncf(m,n,nsmin:nsmax) =  n*nfp*bsubsmnsf(m,n,nsmin:nsmax) - db_vdsc(m,n,nsmin:nsmax)                     
           ksupvmncf(m,n,nsmin:nsmax) = -m*bsubsmnsf(m,n,nsmin:nsmax)     + db_udsc(m,n,nsmin:nsmax)
         ENDDO
      ENDDO
      CALL second0(toff)

      IF (.NOT. l_asym) GOTO 2000

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+2)
      ELSE
         nsmax=MIN(ns,endglobrow+1)
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!
      ALLOCATE(bsubsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),                 & ! First, get BSUBMN on half mesh
               bsubumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp_par(bsubsijh, bsubsmnch, f_cos, f_full, 1, ns_span)      ! Calculate BSUB?MN?H
      CALL tomnsp_par(bsubuijh, bsubumnsh, f_sin, f_full, 1, ns_span)
      CALL tomnsp_par(bsubvijh, bsubvmnsh, f_sin, f_full, 1, ns_span)

!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF
!

      ALLOCATE(bsubsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),                 & ! Next, get BSUBMN on full mesh
               bsubumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                 &
               bsubvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation3 failed in CV_CURRENTS'

      CALL to_full_mesh_par(bsubsmnch, bsubsmncf, nblock, nsmin, nsmax)
      CALL to_full_mesh_par(bsubumnsh, bsubumnsf, nblock, nsmin, nsmax)
      CALL to_full_mesh_par(bsubvmnsh, bsubvmnsf, nblock, nsmin, nsmax)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      IF (nsmin .EQ. 1) THEN
         bsubsmncf(:,:,1) = 0
         bsubumnsf(:,:,1) = 0
         bsubsmncf(m1,:,1) = bsubsmnch(m1,:,1)
         bsubumnsf(m1,:,1) = bsubumnsh(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
         bsubvmnsf(:,:,1) = 0
         bsubvmnsf(m0:m1,:,1) = bsubvmnsh(m0:m1,:,2)
      END IF

!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!
!
      ALLOCATE(db_udss(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
               db_vdss(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation4 failed in CV_CURRENTS'
      CALL second0(toff)

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull_par(db_udss, bsubumnsh, nblock, nsmin, nsmax)
      CALL GradientFull_par(db_vdss, bsubvmnsh, nblock, nsmin, nsmax)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      IF (nsmin .EQ. 1) THEN
         db_udss(m0:m1,:,1) = 2*ohs*bsubumnsh(m0:m1,:,2)
         bsubvmnsf(m1,:,1) = bsubvmnsh(m1,:,2)
         db_vdss(m1,:,1)    = 2*ohs*bsubvmnsh(m1,:,2)                      !--> 0
         bsubumnsf(m1,:,1) = -r0*bsubsmncf(m1,:,1)                         !--> 0
      END IF

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF

!     NOTE: ksupsmncf(ns) is REALLY ksupsmnch(ns) a half-grid point in from edge
!           to satisfy k^s = -r0*k^u at s=r0, for m=1, must have b_u = -r0*b_s for m=1

      DO n = -ntor, ntor
         DO m = 0, mpol
            ksupsmncf(m,n,nsmin:nsmax) =  m*bsubvmnsf(m,n,nsmin:nsmax) - n*nfp*bsubumnsf(m,n,nsmin:nsmax) ! dB_v/du - dB_u/dv  ! dSin ->  Cos
            ksupumnsf(m,n,nsmin:nsmax) = -n*nfp*bsubsmncf(m,n,nsmin:nsmax) - db_vdss(m,n,nsmin:nsmax)     ! dB_s/dv - dB_v/ds  ! dCos -> -Sin
            ksupvmnsf(m,n,nsmin:nsmax) =  m*bsubsmncf(m,n,nsmin:nsmax)     + db_udss(m,n,nsmin:nsmax)     ! dB_u/ds - dB_s/du  ! dCos -> -Sin
         END DO
      END DO

 2000 CONTINUE

      IF (nsmin .EQ. 1) THEN
!Diagnostic output: K^s-r0*K^u at origin BEFORE bc applied
         IF (l_update_state) THEN
            temp = (ksupsmnsf(m1,:,1) + r0*ksupumncf(m1,:,1))/2
            jsju_ratio_s = SUM((2*temp)**2)
            IF (jsju_ratio_s .GT. zero) THEN
               jsju_ratio_s = SUM((ksupsmnsf(m1,:,1) - r0*ksupumncf(m1,:,1))**2) / jsju_ratio_s
               jsju_ratio_s = SQRT(jsju_ratio_s)
            END IF

            IF (l_asym) THEN
            temp = (ksupsmncf(m1,:,1) - r0*ksupumnsf(m1,:,1))/2
            jsju_ratio_a = SUM((2*temp)**2)
            IF (jsju_ratio_a .GT. zero) THEN
               jsju_ratio_a = SUM((ksupsmncf(m1,:,1) + r0*ksupumnsf(m1,:,1))**2) / jsju_ratio_a
               jsju_ratio_a = SQRT(jsju_ratio_a)
            END IF
            END IF
         END IF

!Apply exact bc K^s = r0*K^u for m=1 at origin
         IF (l_natural) ksupsmnsf(m1,:,1) = r0*ksupumncf(m1,:,1)
         ksupsmnsf(m0,:,1) = 0;  ksupsmnsf(m2:,:,1) = 0
         ksupumncf(m0,:,1) = 0;  ksupumncf(m2:,:,1) = 0
         ksupvmncf(m1:,:,1)= 0
 
         IF (l_asym) THEN
!Apply exact bc K^s = -r0*K^u for m=1 at origin
         IF (l_natural) ksupsmncf(m1,:,1) = -r0*ksupumnsf(m1,:,1)
         ksupsmncf(m0,:,1) = 0;  ksupsmncf(m2:,:,1) = 0
         ksupumnsf(m0,:,1) = 0;  ksupumnsf(m2:,:,1) = 0
         ksupvmnsf(m1:,:,1)= 0
         END IF

      END IF

!
!     Compute spectral truncation error measure
!     (Use ksupsijf for temporary FILTERED bsubXijh values)
!
      IF (l_update_state) THEN
         IF (nsmax .EQ. ns) THEN 
         r0 = SQRT(SUM(ksupsmnsf(:,:,ns)**2))
         IF (l_asym) r0 = r0 + SQRT(SUM(ksupsmncf(:,:,ns)**2))
         r0 = r0/SUM(bsubvmncf(:,:,ns)**2 + bsubumncf(:,:,ns)**2)
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
         CALL toijsp_par(bsubsmnsh(:,:,n1:n2), ksupsijf, f_none, f_sin, f_full, 1, ns_span)
         CALL toijsp_par(bsubumnch(:,:,n1:n2), ksupuijf, f_none, f_cos, f_full, 1, ns_span)
         CALL toijsp_par(bsubvmnch(:,:,n1:n2), ksupvijf, f_none, f_cos, f_full, 1, ns_span)
         IF (l_asym) THEN
         CALL toijsp_par(bsubsmnch(:,:,n1:n2), ksupsijf, f_sum, f_cos, f_full, 1, ns_span)
         CALL toijsp_par(bsubumnsh(:,:,n1:n2), ksupuijf, f_sum, f_sin, f_full, 1, ns_span)
         CALL toijsp_par(bsubvmnsh(:,:,n1:n2), ksupvijf, f_sum, f_sin, f_full, 1, ns_span)
         END IF
         stes(1) = SUM((ksupsijf(:,:,1:ns_span)-bsubsijh(:,:,n1:n2))**2)
         stes(2) = SUM((ksupsijf(:,:,1:ns_span)+bsubsijh(:,:,n1:n2))**2)         
         stes(3) = SUM((ksupuijf(:,:,1:ns_span)-bsubuijh(:,:,n1:n2))**2)
         stes(4) = SUM((ksupuijf(:,:,1:ns_span)+bsubuijh(:,:,n1:n2))**2)         
         stes(5) = SUM((ksupvijf(:,:,1:ns_span)-bsubvijh(:,:,n1:n2))**2)
         stes(6) = SUM((ksupvijf(:,:,1:ns_span)+bsubvijh(:,:,n1:n2))**2)         
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(stes,rtes,6,MPI_REAL8,MPI_SUM,               &
                            MPI_COMM_WORLD,MPI_ERR)
#else
         rtes = stes
#endif
         IF (rtes(2) .NE. zero) ste(2) = SQRT(rtes(1)/rtes(2))
         IF (rtes(4) .NE. zero) ste(3) = SQRT(rtes(3)/rtes(4))
         IF (rtes(6) .NE. zero) ste(4) = SQRT(rtes(5)/rtes(6))
      END IF

      DEALLOCATE (bsubsmnsh, bsubumnch, bsubvmnch,                       &
                  bsubsmnsf, bsubumncf, bsubvmncf,                       &
                  db_vdsc,   db_udsc)
      IF (l_asym) THEN
      DEALLOCATE (bsubsmnch, bsubumnsh, bsubvmnsh,                       &
                  bsubsmncf, bsubumnsf, bsubvmnsf,                       &
                  db_vdss, db_udss)
      END IF


      DEALLOCATE (bsubsijh, bsubuijh, bsubvijh) 
 
      CALL second0(toff)

      ns_span = nsmax-nsmin+1
      CALL toijsp_par(ksupsmnsf(:,:,nsmin:), ksupsijf, f_none, f_sin, f_full, 1, ns_span)
      CALL toijsp_par(ksupumncf(:,:,nsmin:), ksupuijf, f_none, f_cos, f_full, 1, ns_span)
      CALL toijsp_par(ksupvmncf(:,:,nsmin:), ksupvijf, f_none, f_cos, f_full, 1, ns_span)
      IF (l_asym) THEN
      CALL toijsp_par(ksupsmncf(:,:,nsmin:), ksupsijf, f_sum, f_cos, f_full, 1, ns_span)
      CALL toijsp_par(ksupumnsf(:,:,nsmin:), ksupuijf, f_sum, f_sin, f_full, 1, ns_span)
      CALL toijsp_par(ksupvmnsf(:,:,nsmin:), ksupvijf, f_sum, f_sin, f_full, 1, ns_span)
      END IF


      IF (lcurr) THEN
         CALL tolowerf_par(ksupsijf, ksupuijf, ksupvijf,                 & 
                           ksubsijf, ksubuijf, ksubvijf, nsmin, nsmax)
         IF (nsmin .EQ. 1) ksubuijf(:,:,1) = 0 
      ENDIF


      CALL second0(toff)
      skstoff = toff
      cv_current_time=cv_current_time+(skstoff-skston)

      time_current = time_current+(skstoff-skston)

      END SUBROUTINE cv_currents_par
!#endif

      SUBROUTINE cv_currents (bsupsijh, bsupuijh, bsupvijh,              &
                              ksupsijf, ksupuijf, ksupvijf,              &
                              lmagen, lcurr)
!     
!     WRITTEN 02-26-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes contravariant currents on full-grid from 
!              contravariant B-fields on half-grid
!
      USE stel_kinds
      USE fourier 
      USE quantities
      USE nscalingtools, ONLY: startglobrow, endglobrow

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(IN)  ::               &
           bsupsijh, bsupuijh, bsupvijh
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(OUT) ::               &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(IN)  :: lmagen, lcurr

      INTEGER     :: n1, n2

      n1 = startglobrow
      n2 = endglobrow

      startglobrow = 1
      endglobrow = ns

      CALL cv_currents_par (bsupsijh, bsupuijh, bsupvijh,                &
                            ksupsijf, ksupuijf, ksupvijf,                &
                            lmagen, lcurr)

      startglobrow = n1
      endglobrow = n2

      END SUBROUTINE cv_currents
