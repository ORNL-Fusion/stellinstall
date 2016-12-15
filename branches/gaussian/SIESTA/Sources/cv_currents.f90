!*******************************************************************************
!>  \file cv_currents.f90
!>  \brief Contains routines to compute the current density from Curl(B).
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

#if defined(SKS)
!*******************************************************************************
!>  \brief Compute current density in parallel.
!>
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
!>  \param[in]  lmagen   UNKNOWN
!>  \param[in]  lcurr    UNKNOWN
!*******************************************************************************
      SUBROUTINE cv_currents_par (bsupsijh, bsupuijh, bsupvijh,      &
                                  ksupsijf, ksupuijf, ksupvijf,      &
                                  lmagen, lcurr)
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i
      USE evolution, ONLY: l_natural, jsju_ratio, ste 
      USE shared_data, ONLY: l_update_state
      USE timer_mod
      USE Hessian, ONLY: mupar_norm
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(IN)  ::           &  !FOR NOW, ASSUME THESE ARE ARRAYS FOR ALL NS!
           bsupsijh, bsupuijh, bsupvijh
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(OUT) ::           &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(in)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1, m2=2
      INTEGER :: m, n, istat, l, js, nblock
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                     &
         bsubsijh, bsubuijh, bsubvijh,                                  &
#if defined(MRC)
         bsubsmnch, bsubumnsh, bsubvmnsh,                               &
         bsubsmncf, bsubumnsf, bsubvmnsf,                               &
#endif
         bsubsmnsh, bsubumnch, bsubvmnch,                               &
         bsubsmnsf, bsubumncf, bsubvmncf,                               &
         db_uds, db_vds, bsq_loc, tempij
      REAL(rprec) :: beta, r0, temp(-ntor:ntor), skston, skstoff,       &
                     ton, toff, stes(6), rtes(6), wb0
      INTEGER     :: k, nsmin, nsmax, ns_span, n1, n2
!-----------------------------------------------
      CALL second0(ton)
      skston=ton

      nsmin=MAX(1, startglobrow); nsmax=MIN(ns,endglobrow+1)
      n1=MAX(1, startglobrow); n2=MIN(ns,endglobrow)

      IF (lcurr) THEN
         nsmin=MAX(1, startglobrow-1); nsmax=MIN(ns,endglobrow+2)
      END IF

      ns_span = nsmax-nsmin+1
!
! Calculate covariant BsubXh on half mesh 
!       
      ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                     &
               bsubuijh(ntheta,nzeta,nsmin:nsmax),                     &
               bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh_par(bsupsijh, bsupuijh, bsupvijh,                &
                        bsubsijh, bsubuijh, bsubvijh,                &
                        nsmin, nsmax)  

      IF (lmagen) THEN
         ALLOCATE (bsq_loc(ntheta,nzeta,nsmin:nsmax), stat=istat)
         bsq_loc = bsupsijh(:,:,1:ns_span)*bsubsijh                   & 
                 + bsupuijh(:,:,1:ns_span)*bsubuijh                   &
                 + bsupvijh(:,:,1:ns_span)*bsubvijh
         IF (n1 .EQ. 1) bsq_loc(:,:,1) = 0
         IF (n1 .GT. 1) THEN
            IF(ANY(bsq_loc(:,:,n1:n2) .LE. zero)) STOP 'BSQ <= 0 in cv_currents_par!'
         END IF

         wb = 0
         DO l = 1, nu_i
            wb = wb + SUM(bsq_loc(l,:,n1:n2)*jacobh(l,:,n1:n2))*cosmui(l,m0)
         END DO         
         wb0 = 2*signjac*(twopi*pi*hs_i)*wb/2                              !2X : theta 0-pi
         CALL MPI_ALLREDUCE(wb0,wb,1,MPI_REAL8, MPI_SUM,                &
                            MPI_COMM_WORLD,MPI_ERR)

!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
         IF (.NOT. ALLOCATED(mupar_norm)) THEN
            ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
            IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
            mupar_norm = 0; k=0
            DO js = nsmin, nsmax
               k = k+1
               IF (js .EQ. 1) CYCLE
               DO l = 1, nu_i
                  mupar_norm(js) = mupar_norm(js) +                     &
                  SUM(bsupvijh(l,:,k)**2/bsq_loc(l,:,js))*cosmui(l,m0)
               END DO
            END DO
            beta = MAX(1.E-5_dp, wp_i/wb_i)
            IF (nsmin .EQ. 1) mupar_norm(1) = mupar_norm(2)
            mupar_norm = -signjac*beta*mupar_norm
         END IF

         DEALLOCATE (bsq_loc)
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!       
      ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                & ! First, get BSUBMN on half mesh
               bsubumnch(0:mpol,-ntor:ntor,nsmin:nsmax),                &
               bsubvmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp_par(bsubsijh, bsubsmnsh, 1, ifull, 1, ns_span)    ! Calculate BSUB?MN?H
      CALL tomnsp_par(bsubuijh, bsubumnch, 0, ifull, 1, ns_span)
      CALL tomnsp_par(bsubvijh, bsubvmnch, 0, ifull, 1, ns_span)
      CALL second0(toff)

!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF 
!        
                    
      ALLOCATE(bsubsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &  ! Next, get BSUBMN on full mesh
               bsubumncf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
               bsubvmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
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
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                    &
               db_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation4 failed in CV_CURRENTS'
      CALL second0(toff)

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull_par(db_uds, bsubumnch, nblock, nsmin, nsmax)
      CALL GradientFull_par(db_vds, bsubvmnch, nblock, nsmin, nsmax)

      CALL second0(toff)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      IF (nsmin .EQ. 1) THEN
         r0 = hs_i/2
         db_uds(m0:m1,:,1) = 2*ohs*bsubumnch(m0:m1,:,2)
         bsubvmncf(m1,:,1) = bsubvmnch(m1,:,2)
         db_vds(m1,:,1)    = 2*ohs*bsubvmnch(m1,:,2)         !--> 0
         bsubumncf(m1,:,1) = r0*bsubsmnsf(m1,:,1)           !--> 0
      END IF

!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
      CALL second0(toff)

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF
      ns_span = nsmax-nsmin+1
      DO n = -ntor, ntor
         DO m = 0, mpol
           ksupsmnsf(m,n,nsmin:nsmax) = -m*bsubvmncf(m,n,nsmin:nsmax) + n*nfp*bsubumncf(m,n,nsmin:nsmax) ! dB_v/du - dB_u/dv  ! dCos -> -Sin
           ksupumncf(m,n,nsmin:nsmax) =  n*nfp*bsubsmnsf(m,n,nsmin:nsmax) - db_vds(m,n,nsmin:nsmax)      ! dB_s/dv - dB_v/ds  ! dSin ->  Cos
           ksupvmncf(m,n,nsmin:nsmax) = -m*bsubsmnsf(m,n,nsmin:nsmax)     + db_uds(m,n,nsmin:nsmax)      ! dB_u/ds - dB_s/du  ! dSin ->  Cos
         ENDDO
      ENDDO
      CALL second0(toff)

      IF (nsmin .EQ. 1) THEN
!Diagnostic output: K^s-r0*K^u at origin BEFORE bc applied
         temp = (ksupsmnsf(m1,:,1) + r0*ksupumncf(m1,:,1))/2
         IF (l_update_state) THEN
            jsju_ratio = SUM((2*temp)**2)
            IF (jsju_ratio .GT. zero) THEN
               jsju_ratio = SUM((ksupsmnsf(m1,:,1) - r0*ksupumncf(m1,:,1))**2) / jsju_ratio
               jsju_ratio = SQRT(jsju_ratio)
            END IF
         END IF

!Apply exact bc K^s = r0*K^u for m=1 at origin
         IF (l_natural) ksupsmnsf(m1,:,1) = r0*ksupumncf(m1,:,1)
         
         ksupsmnsf(m0,:,1) = 0;  ksupsmnsf(m2:,:,1) = 0
         ksupumncf(m0,:,1) = 0;  ksupumncf(m2:,:,1) = 0
         ksupvmncf(m1:,:,1)= 0
      END IF

      DEALLOCATE (bsubsmnsf, bsubumncf, bsubvmncf)

#if defined(MRC)
      IF (l_force_asym) THEN
!
! Calculate BsubXmnh Fourier coefficients
!
         ALLOCATE(bsubsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),                & ! First, get BSUBMN on half mesh
                  bsubumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                  bsubvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

         CALL tomnsp_par(bsubsijh, bsubsmnch, fourier_parity_even, fourier_full_grid, 1, ns_span)    ! Calculate BSUB?MN?H
         CALL tomnsp_par(bsubuijh, bsubumnsh, fourier_parity_odd, fourier_full_grid, 1, ns_span)
         CALL tomnsp_par(bsubvijh, bsubvmnsh, fourier_parity_odd, fourier_full_grid, 1, ns_span)

!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF
!

         ALLOCATE(bsubsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),                &  ! Next, get BSUBMN on full mesh
                  bsubumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                  bsubvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'

         nblock=SIZE(bsubsmnch,1)*SIZE(bsubsmnch,2)
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

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
         CALL GradientFull_par(db_uds, bsubumnsh, nblock, nsmin, nsmax)
         CALL GradientFull_par(db_vds, bsubvmnsh, nblock, nsmin, nsmax)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
         IF (nsmin .EQ. 1) THEN
            r0 = hs_i/2
            db_uds(m0:m1,:,1) = 2*ohs*bsubumnsh(m0:m1,:,2)
            bsubvmnsf(m1,:,1) = bsubvmnsh(m1,:,2)
            db_vds(m1,:,1)    = 2*ohs*bsubvmnsh(m1,:,2)         !--> 0
            bsubumnsf(m1,:,1) = r0*bsubsmncf(m1,:,1)           !--> 0
         END IF

! MRC:  Check Signs !!!!!!!!!!
!     NOTE: ksupsmncf(ns) is REALLY ksupsmnch(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1

         DO n = -ntor, ntor
            DO m = 0, mpol
               ksupsmncf(m,n,nsmin:nsmax) =  m*bsubvmnsf(m,n,nsmin:nsmax) - n*nfp*bsubumnsf(m,n,nsmin:nsmax) ! dB_v/du - dB_u/dv  ! dSin ->  Cos
               ksupumnsf(m,n,nsmin:nsmax) = -n*nfp*bsubsmncf(m,n,nsmin:nsmax) - db_vds(m,n,nsmin:nsmax)      ! dB_s/dv - dB_v/ds  ! dCos -> -Sin
               ksupvmnsf(m,n,nsmin:nsmax) =  m*bsubsmncf(m,n,nsmin:nsmax)     + db_uds(m,n,nsmin:nsmax)      ! dB_u/ds - dB_s/du  ! dCos -> -Sin
            END DO
         END DO

         IF (nsmin .EQ. 1) THEN
!Apply exact bc K^s = r0*K^u for m=1 at origin
            IF (l_natural) ksupsmncf(m1,:,1) = r0*ksupumnsf(m1,:,1)

            ksupsmncf(m0,:,1) = 0;  ksupsmncf(m2:,:,1) = 0
            ksupumnsf(m0,:,1) = 0;  ksupumnsf(m2:,:,1) = 0
            ksupvmnsf(m1:,:,1)= 0
         END IF

         DEALLOCATE (bsubsmncf, bsubumnsf, bsubvmnsf)
      END IF
#endif

      DEALLOCATE (db_vds, db_uds)

!
!     Compute spectral truncation error measure
!     (Use ksupsijf for temporary FILTERED bsubXijh values)
!
      IF (l_update_state) THEN
         ns_span = n2-n1+1
         CALL toijsp_par(bsubsmnsh(:,:,n1:n2), ksupsijf, 0, 0, 1, 0, 1, ns_span)
         CALL toijsp_par(bsubumnch(:,:,n1:n2), ksupuijf, 0, 0, 0, 0, 1, ns_span)
         CALL toijsp_par(bsubvmnch(:,:,n1:n2), ksupvijf, 0, 0, 0, 0, 1, ns_span)

#if defined(MRC)
         IF (l_force_asym) THEN
            ALLOCATE(tempij(LBOUND(ksupsijf, 1):UBOUND(ksupsijf, 1),           &
                            LBOUND(ksupsijf, 2):UBOUND(ksupsijf, 2),           &
                            LBOUND(ksupsijf, 3):UBOUND(ksupsijf, 3)))

            CALL toijsp_par(bsubsmnch(:,:,n1:n2), tempij, fourier_none, fourier_none, &
                            fourier_parity_even, fourier_full_grid, 1, ns_span)
            ksupsijf = ksupsijf + tempij

            CALL toijsp_par(bsubumnsh(:,:,n1:n2), tempij, fourier_none, fourier_none, &
                            fourier_parity_odd, fourier_full_grid, 1, ns_span)
            ksupuijf = ksupuijf + tempij

            CALL toijsp_par(bsubvmnsh(:,:,n1:n2), tempij, fourier_none, fourier_none, &
                            fourier_parity_odd, fourier_full_grid, 1, ns_span)
            ksupvijf = ksupvijf + tempij

            DEALLOCATE(tempij)
         END IF
#endif

         stes(1) = SUM((ksupsijf(:,:,1:ns_span)-bsubsijh(:,:,n1:n2))**2)
         stes(2) = SUM((ksupsijf(:,:,1:ns_span)+bsubsijh(:,:,n1:n2))**2)
         stes(3) = SUM((ksupuijf(:,:,1:ns_span)-bsubuijh(:,:,n1:n2))**2)
         stes(4) = SUM((ksupuijf(:,:,1:ns_span)+bsubuijh(:,:,n1:n2))**2)         
         stes(5) = SUM((ksupvijf(:,:,1:ns_span)-bsubvijh(:,:,n1:n2))**2)
         stes(6) = SUM((ksupvijf(:,:,1:ns_span)+bsubvijh(:,:,n1:n2))**2)         
         CALL MPI_ALLREDUCE(stes,rtes,6,MPI_REAL8,MPI_SUM,              &
                            MPI_COMM_WORLD,MPI_ERR)
         IF (rtes(2) .NE. zero) ste(2) = SQRT(rtes(1)/rtes(2))
         IF (rtes(4) .NE. zero) ste(3) = SQRT(rtes(3)/rtes(4))
         IF (rtes(6) .NE. zero) ste(4) = SQRT(rtes(5)/rtes(6))
         ns_span = nsmax-nsmin+1
      END IF

      DEALLOCATE (bsubsmnsh, bsubumnch, bsubvmnch,                      &
                  bsubsijh, bsubuijh, bsubvijh) 

#if defined(MRC)
      IF (l_force_asym) THEN
         DEALLOCATE (bsubsmnch, bsubumnsh, bsubvmnsh)
      END IF
#endif

      CALL second0(toff)

      CALL toijsp_par(ksupsmnsf(:,:,nsmin:), ksupsijf, 0, 0, 1, 0, 1, ns_span) 
      CALL toijsp_par(ksupumncf(:,:,nsmin:), ksupuijf, 0, 0, 0, 0, 1, ns_span)
      CALL toijsp_par(ksupvmncf(:,:,nsmin:), ksupvijf, 0, 0, 0, 0, 1, ns_span)

#if defined(MRC)
      IF (l_force_asym) THEN
         ALLOCATE(tempij(LBOUND(ksupsijf, 1):UBOUND(ksupsijf, 1),              &
                         LBOUND(ksupsijf, 2):UBOUND(ksupsijf, 2),              &
                         LBOUND(ksupsijf, 3):UBOUND(ksupsijf, 3)))

         CALL toijsp_par(ksupsmncf(:,:,nsmin:), tempij, fourier_none, fourier_none, &
                         fourier_parity_even, fourier_full_grid, 1, ns_span)
         ksupsijf = ksupsijf + tempij

         CALL toijsp_par(ksupumnsf(:,:,nsmin:), tempij, fourier_none, fourier_none, &
                         fourier_parity_odd, fourier_full_grid, 1, ns_span)
         ksupuijf = ksupuijf + tempij

         CALL toijsp_par(ksupvmnsf(:,:,nsmin:), tempij, fourier_none, fourier_none, &
                         fourier_parity_odd, fourier_full_grid, 1, ns_span)
         ksupvijf = ksupvijf + tempij

         DEALLOCATE(tempij)
      END IF
#endif

      IF (lcurr) THEN
         CALL tolowerf_par(ksupsijf, ksupuijf, ksupvijf,             & 
                           ksubsijf, ksubuijf, ksubvijf, nsmin, nsmax)
!DON'T KNOW WHY THIS WOULD BE NEEDED?
!         IF (nsmax .EQ. ns) THEN
!            ksubsijf(:,:,ns) = 0
!            ksubuijf(:,:,ns) = 0
!            ksubvijf(:,:,ns) = 0
!         END IF
         IF (nsmin .EQ. 1) ksubuijf(:,:,1) = 0 
      ENDIF


      CALL second0(toff)
      skstoff = toff
      cv_current_time=cv_current_time+(skstoff-skston)

      time_current = time_current+(skstoff-skston)

      END SUBROUTINE cv_currents_par
#endif

!*******************************************************************************
!>  \brief Compute current density in parallel.
!>
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
!>  \param[in]  lmagen   UNKNOWN
!>  \param[in]  lcurr    UNKNOWN
!*******************************************************************************
      SUBROUTINE cv_currents (bsupsijh, bsupuijh, bsupvijh,          &
                              ksupsijf, ksupuijf, ksupvijf,          &
                              lmagen, lcurr)
!     
!     WRITTEN 02-26-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes contravariant currents on full-grid from 
!              contravariant B-fields on half-grid
!
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i
      USE evolution, ONLY: l_natural, jsju_ratio, ste
      USE shared_data, ONLY: l_update_state
      USE Hessian, ONLY: mupar_norm
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(in)  ::           &
           bsupsijh, bsupuijh, bsupvijh
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(out) ::           &
           ksupsijf, ksupuijf, ksupvijf
      LOGICAL, INTENT(in)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1, m2=2
      INTEGER :: m, n, istat, l, js, i, j, k, nblock
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                     &
         bsubsijh, bsubuijh, bsubvijh,                                  &
         bsubsmnsh, bsubumnch, bsubvmnch,                               &
         bsubsmnsf, bsubumncf, bsubvmncf,                               &
#if defined(MRC)
         bsubsmnch, bsubumnsh, bsubvmnsh,                               &
         bsubsmncf, bsubumnsf, bsubvmnsf,                               &
#endif
         db_uds, db_vds, tempij
      REAL(rprec) :: beta, r0, temp(-ntor:ntor)
      REAL(rprec) :: skston, skstoff
      INTEGER :: nsmin, nsmax
!-----------------------------------------------
!
! Calculate BsubXh on half mesh 
!       
      CALL second0(skston)

      ALLOCATE(bsubsijh(ntheta,nzeta,ns),                              &
               bsubuijh(ntheta,nzeta,ns),                              &
               bsubvijh(ntheta,nzeta,ns), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                    &
                    bsubsijh, bsubuijh, bsubvijh)  

      IF (lmagen) THEN
         bsq = bsupsijh*bsubsijh + bsupuijh*bsubuijh                &
             + bsupvijh*bsubvijh
         IF (ANY(bsq(:,:,2) .le. zero)) STOP 'BSQ <= 0 in cv_currents!'

         wb = zero
         DO l = 1, nu_i
            wb = wb + SUM(bsq(l,:,2:ns)*jacobh(l,:,2:ns))*cosmui(l,m0)
         END DO         
         wb = (signjac*twopi*twopi*hs_i)*wb/2 

!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
         IF (.NOT. ALLOCATED(mupar_norm)) THEN
            ALLOCATE(mupar_norm(ns), stat=istat)
            IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
            mupar_norm = 0
            DO js = 2, ns
               DO l = 1, nu_i
                  mupar_norm(js) = mupar_norm(js) +                     &
                  SUM(bsupvijh(l,:,js)**2/bsq(l,:,js))*cosmui(l,m0)
               END DO
            END DO
            beta = MAX(1.E-5_dp, wp_i/wb_i)
            mupar_norm(1) = mupar_norm(2)
            mupar_norm = -signjac*beta*mupar_norm
         END IF
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!       
      ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,ns),                         & ! First, get BSUBMN on half mesh
               bsubumnch(0:mpol,-ntor:ntor,ns),                         &
               bsubvmnch(0:mpol,-ntor:ntor,ns), stat=istat)   
      IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp(bsubsijh, bsubsmnsh, 1, ifull)                         ! Calculate BSUB?MN?H
      CALL tomnsp(bsubuijh, bsubumnch, 0, ifull)
      CALL tomnsp(bsubvijh, bsubvmnch, 0, ifull)
!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF 
!        
      ALLOCATE(bsubsmnsf(0:mpol,-ntor:ntor,ns),                         &  ! Next, get BSUBMN on full mesh
               bsubumncf(0:mpol,-ntor:ntor,ns),                         &
               bsubvmncf(0:mpol,-ntor:ntor,ns), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
      nblock=SIZE(bsubsmnsh,1)*SIZE(bsubsmnsh,2)
      CALL to_full_mesh(bsubsmnsh, bsubsmnsf, nblock)                       
      CALL to_full_mesh(bsubumnch, bsubumncf, nblock)
      CALL to_full_mesh(bsubvmnch, bsubvmncf, nblock)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      bsubsmnsf(:,:,1) = 0
      bsubumncf(:,:,1) = 0
      bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,1)
      bsubumncf(m1,:,1) = bsubumnch(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS
!      bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,2)
!      bsubumncf(m1,:,1) = bsubumnch(m1,:,2)

      r0 = hs_i/2
!      bsubumncf(m1,-ntor:-1,1) = r0*bsubsmnsf(m1,-ntor:-1,1)
!      bsubumncf(m1,1:ntor,1)   = r0*bsubsmnsf(m1,1:ntor,1)

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
      bsubvmncf(:,:,1) = 0
      bsubvmncf(m0:m1,:,1) = bsubvmnch(m0:m1,:,2)

      bsubsmnsf(:,:,ns) = bsubsmnsh(:,:,ns)
      bsubumncf(:,:,ns) = bsubumnch(:,:,ns)
      bsubvmncf(:,:,ns) = bsubvmnch(:,:,ns)

!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,ns),                             &
               db_vds(0:mpol,-ntor:ntor,ns), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation4 failed in CV_CURRENTS'

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull(db_uds, bsubumnch, nblock)
      CALL GradientFull(db_vds, bsubvmnch, nblock)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      db_uds(m0:m1,:,1) = 2*ohs*bsubumnch(m0:m1,:,2)
      db_vds(m1,:,1) = 2*ohs*bsubvmnch(m1,:,2)
      bsubvmncf(m1,:,1) = bsubvmnch(m1,:,2)
      bsubumncf(m1,:,1) = r0*bsubsmnsf(m1,:,1)

!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
      DO m = 0, mpol
         DO n = -ntor, ntor
           ksupsmnsf(m,n,:) = -m*bsubvmncf(m,n,:)     + n*nfp*bsubumncf(m,n,:) ! dB_v/du - dB_u/dv  ! dCos -> -Sin
           ksupumncf(m,n,:) =  n*nfp*bsubsmnsf(m,n,:) - db_vds(m,n,:)          ! dB_s/dv - dB_v/ds  ! dSin ->  Cos
           ksupvmncf(m,n,:) = -m*bsubsmnsf(m,n,:)     + db_uds(m,n,:)          ! dB_u/ds - dB_s/du  ! dSin ->  Cos
         ENDDO
      ENDDO

!Diagnostic output: K^s-r0*K^u at origin BEFORE bc applied
      temp = (ksupsmnsf(m1,:,1) + r0*ksupumncf(m1,:,1))/2
      IF (l_update_state) THEN
        jsju_ratio = SUM((2*temp)**2)
        IF (jsju_ratio .GT. zero) THEN
          jsju_ratio = SUM((ksupsmnsf(m1,:,1) - r0*ksupumncf(m1,:,1))**2) / jsju_ratio
          jsju_ratio = SQRT(jsju_ratio)
        END IF
      END IF

!Apply exact bc K^s = r0*K^u for m=1 at origin
      IF (l_natural) ksupsmnsf(m1,:,1) = r0*ksupumncf(m1,:,1)

      ksupsmnsf(m0,:,1) = 0;  ksupsmnsf(m2:,:,1) = 0
      ksupumncf(m0,:,1) = 0;  ksupumncf(m2:,:,1) = 0
      ksupvmncf(m1:,:,1)= 0

      DEALLOCATE(bsubsmnsf, bsubumncf, bsubvmncf)

#if defined(MRC)
      IF (l_force_asym) THEN
!
! Calculate BsubXmnh Fourier coefficients
!
         ALLOCATE(bsubsmnch(0:mpol,-ntor:ntor,ns),                         & ! First, get BSUBMN on half mesh
                  bsubumnsh(0:mpol,-ntor:ntor,ns),                         &
                  bsubvmnsh(0:mpol,-ntor:ntor,ns), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

         CALL tomnsp(bsubsijh, bsubsmnch, fourier_parity_even, fourier_full_grid)                         ! Calculate BSUB?MN?H
         CALL tomnsp(bsubuijh, bsubumnsh, fourier_parity_odd, fourier_full_grid)
         CALL tomnsp(bsubvijh, bsubvmnsh, fourier_parity_odd, fourier_full_grid)
!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF
!
         ALLOCATE(bsubsmncf(0:mpol,-ntor:ntor,ns),                         &  ! Next, get BSUBMN on full mesh
                  bsubumnsf(0:mpol,-ntor:ntor,ns),                         &
                  bsubvmnsf(0:mpol,-ntor:ntor,ns), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'

         CALL to_full_mesh(bsubsmnch, bsubsmncf, nblock)
         CALL to_full_mesh(bsubumnsh, bsubumnsf, nblock)
         CALL to_full_mesh(bsubvmnsh, bsubvmnsf, nblock)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
         bsubsmncf(:,:,1) = 0
         bsubumnsf(:,:,1) = 0
         bsubsmncf(m1,:,1) = bsubsmnch(m1,:,1)
         bsubumnsf(m1,:,1) = bsubumnsh(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS
!      bsubsmncf(m1,:,1) = bsubsmnch(m1,:,2)
!      bsubumnsf(m1,:,1) = bsubumnsh(m1,:,2)

!      bsubumnsf(m1,-ntor:-1,1) = r0*bsubsmncf(m1,-ntor:-1,1)
!      bsubumnsf(m1,1:ntor,1)   = r0*bsubsmncf(m1,1:ntor,1)

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
         bsubvmnsf(:,:,1) = 0
         bsubvmnsf(m0:m1,:,1) = bsubvmnsh(m0:m1,:,2)

         bsubsmncf(:,:,ns) = bsubsmnch(:,:,ns)
         bsubumnsf(:,:,ns) = bsubumnsh(:,:,ns)
         bsubvmnsf(:,:,ns) = bsubvmnsh(:,:,ns)

!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!
!

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
         CALL GradientFull(db_uds, bsubumnsh, nblock)
         CALL GradientFull(db_vds, bsubvmnsh, nblock)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
         db_uds(m0:m1,:,1) = 2*ohs*bsubumnsh(m0:m1,:,2)
         db_vds(m1,:,1) = 2*ohs*bsubvmnsh(m1,:,2)
         bsubvmnsf(m1,:,1) = bsubvmnsh(m1,:,2)
         bsubumnsf(m1,:,1) = r0*bsubsmncf(m1,:,1)

!     NOTE: ksupsmncf(ns) is REALLY ksupsmnch(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
         DO m = 0, mpol
            DO n = -ntor, ntor
               ksupsmncf(m,n,:) =  m*bsubvmnsf(m,n,:)     - n*nfp*bsubumnsf(m,n,:) ! dB_v/du - dB_u/dv  ! dSin ->  Cos
               ksupumnsf(m,n,:) = -n*nfp*bsubsmncf(m,n,:) - db_vds(m,n,:)          ! dB_s/dv - dB_v/ds  ! dCos -> -Sin
               ksupvmnsf(m,n,:) =  m*bsubsmncf(m,n,:)     + db_uds(m,n,:)          ! dB_u/ds - dB_s/du  ! dCos -> -Sin
            ENDDO
         ENDDO

!Apply exact bc K^s = r0*K^u for m=1 at origin
         IF (l_natural) ksupsmncf(m1,:,1) = r0*ksupumnsf(m1,:,1)

         ksupsmncf(m0,:,1) = 0;  ksupsmncf(m2:,:,1) = 0
         ksupumnsf(m0,:,1) = 0;  ksupumnsf(m2:,:,1) = 0
         ksupvmnsf(m1:,:,1)= 0

         DEALLOCATE(bsubsmncf, bsubumnsf, bsubvmnsf)
      END IF
#endif

      DEALLOCATE (db_uds, db_vds)
!
!     Compute spectral truncation error measure
!     (Use ksupsijf for temporary FILTERED bsubXijh values)
!
      IF (l_update_state) THEN
         CALL toijsp(bsubsmnsh, ksupsijf, 0, 0, 1, 0)
         CALL toijsp(bsubumnch, ksupuijf, 0, 0, 0, 0)
         CALL toijsp(bsubvmnch, ksupvijf, 0, 0, 0, 0)

#if defined(MRC)
         IF (l_force_asym) THEN
            ALLOCATE(tempij(ntheta, nzeta, ns))

            CALL toijsp(bsubsmnch, tempij, fourier_none, fourier_none, &
                        fourier_parity_even, fourier_full_grid)
            ksupsijf = ksupsijf + tempij

            CALL toijsp(bsubumnsh, tempij, fourier_none, fourier_none, &
                        fourier_parity_odd, fourier_full_grid)
            ksupuijf = ksupuijf + tempij

            CALL toijsp(bsubvmnsh, tempij, fourier_none, fourier_none, &
                        fourier_parity_odd, fourier_full_grid)
            ksupvijf = ksupvijf + tempij
         END IF
#endif

         r0 = SUM((ksupsijf+bsubsijh)**2)
         IF (r0 .NE. zero) ste(2) = SQRT(SUM((ksupsijf-bsubsijh)**2)/r0)
         r0 = SUM((ksupuijf+bsubuijh)**2)
         IF (r0 .NE. zero) ste(3) = SQRT(SUM((ksupuijf-bsubuijh)**2)/r0)
         r0 = SUM((ksupvijf+bsubvijh)**2)         
         IF (r0 .NE. zero) ste(4) = SQRT(SUM((ksupvijf-bsubvijh)**2)/r0)
      END IF

      DEALLOCATE(bsubsmnsh, bsubumnch, bsubvmnch)
#if defined(MRC)
      IF (l_force_asym) THEN
         DEALLOCATE(bsubsmnch, bsubumnsh, bsubvmnsh)
      END IF
#endif
      DEALLOCATE(bsubsijh, bsubuijh, bsubvijh) 

      CALL toijsp(ksupsmnsf, ksupsijf, 0, 0, 1, 0) 
      CALL toijsp(ksupumncf, ksupuijf, 0, 0, 0, 0)
      CALL toijsp(ksupvmncf, ksupvijf, 0, 0, 0, 0)

#if defined(MRC)
      IF (l_force_asym) THEN
         IF (.not.ALLOCATED(tempij)) ALLOCATE(tempij(ntheta, nzeta, ns))

         CALL toijsp(ksupsmncf, tempij, fourier_none, fourier_none, &
                     fourier_parity_even, fourier_full_grid)
         ksupsijf = ksupsijf + tempij

         CALL toijsp(ksupumnsf, tempij, fourier_none, fourier_none, &
                     fourier_parity_odd, fourier_full_grid)
         ksupuijf = ksupuijf + tempij

         CALL toijsp(ksupvmnsf, tempij, fourier_none, fourier_none, &
                     fourier_parity_odd, fourier_full_grid)
         ksupvijf = ksupvijf + tempij

         DEALLOCATE(tempij)
      END IF
#endif

!     NEED THESE FOR RESISTIVE CONTRIBUTIONS (DIVIDE OUT jac FACTOR LATER)
!     ksubX = jac*JsubX;  ksupX = jac*JsupX
!
      IF (lcurr) THEN
         CALL tolowerf(ksupsijf, ksupuijf, ksupvijf,                  &
                       ksubsijf, ksubuijf, ksubvijf)

!        ksubsijf(:,:,ns) = 0; ksubuijf(:,:,ns) = 0
         ksubuijf(:,:,1)  = 0  
!        ksubvijf(:,:,ns) = 0
      ENDIF

      CALL second0(skstoff)
      time_current = time_current+(skstoff-skston)

      END SUBROUTINE cv_currents
