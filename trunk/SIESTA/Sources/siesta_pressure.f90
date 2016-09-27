!*******************************************************************************
!>  \brief Module contained subroutines for updating pressure as part of the
!>  SIESTA project.
!>  \author S. P. Hirshman and R. Sanchez
!>  \date Aug 18, 2006
!>
!>  Changes in pressure are governed by the continuity equation.
!>     
!>     dn/dt + Div(nv) = dn/dt + n*Div(v) + v.Grad(n) = 0                    (1)
!>
!>  where n is the density and v is the velocity. From the adiabatic principle,
!>
!>     pV^gamma = Constant                                                   (2)
!>
!>  where p is the pressure and V is the volume. From the ideal gas law V scales
!>  as V ~ 1/n. Using this equation (2) can be rewritten in terms of 
!>
!>     p/n^gamma = Constant                                                  (3)
!>
!>  Putting equation (3) into equation (1), the continuty equation can be 
!>  written in terms of the pressure.
!>
!>     dp/dt + gamma*p*Div(v) + v.Grad(p) = 0                                (4)
!>
!>  Using the idenity 
!> 
!>     Div(fA) = f*Div(A) + A.Grad(f)                                        (5)
!>
!>  equation (4) can be written as.
!>
!>     dp/dt = (gamma - 1)v.Grad(p) - gamma*Div(pv)                          (6)
!>
!>  This leads to equation (2.3) in Hirshman et. al. doi:10.1063/1.3597155.
!>  Expanding out equation (6) gives
!>
!>     dp/dt = (gamma - 1)*(v^s*dp/ds + v^u*dp/du + v^v*dp/dv)
!>           - gamma/J*(dJpv^s/ds + dJpv^u/du + dJpv^v/dv)
!*******************************************************************************
      MODULE siesta_pressure

       USE stel_kinds
       USE stel_constants
       USE fourier       
       USE quantities
       USE hessian, ONLY: l_Compute_Hessian
       USE shared_data, ONLY: l_linearize, l_push_s, l_push_u, xPosDef
       USE timer_mod, ONLY: time_update_pres

       CONTAINS
#if defined(SKS) 
!>     \brief Parallel subroutine to advance pressure from t to t+delta_t by updating jpmnch (jacobh*presh)
       SUBROUTINE update_pres_par (delta_t)
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: delta_t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: ifull=0, ihalf=1, pcos=0, psin=1, m0=0, m1=1
       LOGICAL, PARAMETER     :: l_upwind=.FALSE.
       INTEGER                :: m, n, istat, js, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
         work1, work2, work3, work4, work5, work6,                       &
         jvsupuijh, jvsupvijh
       REAL(rprec)            :: ton, toff, gamm1
       INTEGER                :: nsmin, nsmax, ns_span
!-----------------------------------------------

       CALL second0(ton)
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
       ns_span=nsmax-nsmin+1
!
!      First, calculate first term (~gamma) on r.h.s. of pressure evolution equation       
! 
       ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                        & ! Allocate work space
                work2(ntheta,nzeta,nsmin:nsmax),                                 &   
                work3(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

!      Div(p*Flux) terms
!      Note that near axis, jvsups ~ s cos(u+nv), jvsupv -> s sin(nv) due to jacobian element
!      In particular, pijf(1)*jvsupsijf makes no contribution to forces (integrates to zero,
!      since jvsupsijf(1) ~ cos(u) + O(s)!)
!      jvsupu may be finite, although dp/du -> 0 at axis
!      Note: WORK1 has COSINE parity, nsh = ns_i-1


!Upwind to suppress grid separation in radial direction
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*hs_i*                &
             pijf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijf(:,:,nsmin:nsmax))
          work2(:,:,nsmin:nsmax) = work1(:,:,nsmin:nsmax) +             &
             pijh(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
       ELSE
          work2(:,:,nsmin:nsmax) = pijf(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
       END IF

       !!On exit, work1 valid on nsmin+1:nsmax grid
       nblock = SIZE(work1,1)*SIZE(work1,2)
       CALL GradientHalf_par(work1, work2, nblock, nsmin, nsmax)

       ALLOCATE (jvsupuijh(ntheta,nzeta,nsmin:nsmax),                  &   
                 jvsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

       !!On exit, jsupXh valid on nsmin+1:nsmax grid
       CALL to_half_mesh_par(jvsupuijf(:,:,nsmin:nsmax), jvsupuijh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(jvsupvijf(:,:,nsmin:nsmax), jvsupvijh, nblock, nsmin, nsmax)

       ALLOCATE(work4(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work5(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation2 failed in UPDATE_PRES'

       nsmin=nsmin+1   !MAX(1,startglobrow)
       work2(:,:,nsmin:nsmax) = pijh(:,:,nsmin:nsmax)*jvsupuijh(:,:,nsmin:nsmax)
       work3(:,:,nsmin:nsmax) = pijh(:,:,nsmin:nsmax)*jvsupvijh(:,:,nsmin:nsmax)
       
!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
       CALL tomnsp_par(work2, work5, psin, ihalf, 2, ns_span)
       CALL tomnsp_par(work3, work6, psin, ihalf, 2, ns_span)
       CALL tomnsp_par(work1, work4, pcos, ihalf, 2, ns_span)
              
       DO m = 0, mpol
          DO n = -ntor, ntor
             work4(m,n,nsmin:nsmax) = work4(m,n,nsmin:nsmax) +          &
             m*work5(m,n,nsmin:nsmax) + n*nfp*work6(m,n,nsmin:nsmax)      !m: Pol. deriv of SIN(m*theta+n*NFP*phi)
          ENDDO                                                           !n: Tor. deriv of SIN(m*theta+n*NFP*phi)
       ENDDO

       djpmnch(:,:,nsmin:nsmax) = -gamma*work4(:,:,nsmin:nsmax)           !Gamma = adiabatic factor

#if defined(MRC)
       IF (l_force_asym) THEN
!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
          CALL tomnsp_par(work2, work5, fourier_parity_even, fourier_half_grid, 2, ns_span)
          CALL tomnsp_par(work3, work6, fourier_parity_even, fourier_half_grid, 2, ns_span)
          CALL tomnsp_par(work1, work4, fourier_parity_odd, fourier_half_grid, 2, ns_span)

          DO m = 0, mpol
             DO n = -ntor, ntor
                work4(m,n,nsmin:nsmax) = work4(m,n,nsmin:nsmax) -          &
                m*work5(m,n,nsmin:nsmax) - n*nfp*work6(m,n,nsmin:nsmax)      !m: Pol. deriv of COS(m*theta+n*NFP*phi)
             ENDDO                                                           !n: Tor. deriv of COS(m*theta+n*NFP*phi)
          ENDDO

          djpmnsh(:,:,nsmin:nsmax) = -gamma*work4(:,:,nsmin:nsmax)           !Gamma = adiabatic factor
       END IF
#endif

       DEALLOCATE(work5, work6)
!
!      Now, calculate next term ~(gamma-1) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(v) = 0
       nsmin=MAX(1,startglobrow-1)
       work1(:,:,nsmin:nsmax) =                                         &
         pijf_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax) +        &
         pijf_du(:,:,nsmin:nsmax)*jvsupuijf(:,:,nsmin:nsmax) +        &
         pijf_dv(:,:,nsmin:nsmax)*jvsupvijf(:,:,nsmin:nsmax)

       !!On exit, work2 valid on nsmin+1:nsmax grid
       CALL to_half_mesh_par(work1, work2, nblock, nsmin, nsmax)

!      MIXED HALF/FULL: DIAGONAL IN ANGULAR DERIVATIVES
!!       work1 = pijf_ds*jvsupsijf
!!       CALL to_half_mesh_par(work1,work2,nblock,nsmin,nsmax)
!!       work2 = work2+jvsupuijh*pijh_du+jvsupvijh*pijh_dv

       DEALLOCATE(jvsupuijh, jvsupvijh)

!Upwinding correction (makes 1st order streaming dominantly diagonal)
!NOTE: since it enters conservatively, it doesn't effect MHD force, only HESSIAN
!Can multiply this by a positive factor to increase numerical diffusion
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*                     &
             pijf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijf(:,:,nsmin:nsmax))

          IF (nsmin .EQ. 1) work1(:,:,1) = 0
          CALL GradientHalf_par(work3, work1, nblock, nsmin, nsmax)
          work2(:,:,nsmin:nsmax) = work2(:,:,nsmin:nsmax) +             &
                              hs_i*work3(:,:,nsmin:nsmax)
       END IF
       
       nsmin=nsmin+1
       CALL tomnsp_par(work2, work4, pcos, ihalf, 2, ns_span)            ! Calculate Fourier transform of WORK2

!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
       gamm1 = gamma-1
!!       IF (l_Compute_Hessian) gamm1 = gamm1 + xPosDef
       IF (l_linearize) gamm1 = gamm1 + xPosDef

       djpmnch(:,:,nsmin:nsmax) = djpmnch(:,:,nsmin:nsmax) +            &
              gamm1*work4(:,:,nsmin:nsmax)                                ! Harmonics of r.h.s. of Pres. Eq. (full mesh) 

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnch(:,:,nsmin:nsmax) = delta_t*djpmnch(:,:,nsmin:nsmax)    ! calculate increment of pressure harmonics.
!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
       IF (nsmin .EQ. 2) THEN
!          djpmnch(m1:,:,2) = 0             
          djpmnch(:,:,1)  = 0
          djpmnch(m0,:,1) = djpmnch(m0,:,2)                           ! store js=1, m=0 value
       END IF

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL tomnsp_par(work2, work4, pcos, ihalf, 2, ns_span)            ! Calculate Fourier transform of WORK2

          djpmnsh(:,:,nsmin:nsmax) = djpmnsh(:,:,nsmin:nsmax)              &
                                   + gamm1*work4(:,:,nsmin:nsmax) ! Harmonics of r.h.s. of Pres. Eq. (full mesh)

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
          djpmnsh(:,:,nsmin:nsmax) = delta_t*djpmnsh(:,:,nsmin:nsmax)    ! calculate increment of pressure harmonics.
!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
          IF (nsmin .EQ. 2) THEN
!          djpmnsh(m1:,:,2) = 0
             djpmnsh(:,:,1)  = 0
             djpmnsh(m0,:,1) = djpmnsh(m0,:,2)                           ! store js=1, m=0 value
          END IF
       END IF
#endif

       DEALLOCATE(work1, work2, work3, work4)

       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)

       END SUBROUTINE update_pres_par
#endif

!>    \brief Serial subroutine to advance pressure from t to t+delta_t by updating jpmnch (jacobh*presh)
       SUBROUTINE update_pres (delta_t)
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: delta_t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: ifull=0, ihalf=1, pcos=0, psin=1, m0=0, m1=1
       LOGICAL, PARAMETER     :: l_upwind=.FALSE.
       INTEGER :: m, n, istat, js, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
         work1, work2, work3, work4, work5, work6,                       &
         jvsupuijh, jvsupvijh
       REAL(rprec)            :: ton, toff, gamm1
!-----------------------------------------------
       CALL second0(ton)
!
!      First, calculate first term (~gamma) on r.h.s. of pressure evolution equation       
! 
       ALLOCATE(work1(ntheta,nzeta,ns),                                 & ! Allocate work space
                work2(ntheta,nzeta,ns),                                 &   
                work3(ntheta,nzeta,ns), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

!      Div(p*Flux) terms
!      Note that near axis, jvsups ~ s cos(u+nv), jvsupv -> s sin(nv) due to jacobian element
!      In particular, pijf(1)*jvsupsijf makes no contribution to forces (integrates to zero,
!      since jvsupsijf(1) ~ cos(u) + O(s)!)
!      jvsupu may be finite, although dp/du -> 0 at axis
!      Note: WORK1 has COSINE parity, nsh = ns_i-1


!Upwind to suppress grid separation in radial direction
       IF (l_upwind) THEN
          work1 = -0.5_dp*signjac*hs_i*pijf_ds*ABS(jvsupsijf)
          work2 = pijh*jvsupsijf + work1
       ELSE
          work2 = pijf*jvsupsijf
       END IF

       nblock = SIZE(work1,1)*SIZE(work1,2)
       CALL GradientHalf(work1, work2, nblock)
!!     work1(:,:,2:ns) = ohs*(work2(:,:,2:ns) - work2(:,:,1:nsh))          ! BC: JVSUPS(NS) = 0

       ALLOCATE (jvsupuijh(ntheta,nzeta,ns),                           &   
                 jvsupvijh(ntheta,nzeta,ns), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation1 failed in UPDATE_PRES'

       CALL to_half_mesh(jvsupuijf, jvsupuijh, nblock)
       CALL to_half_mesh(jvsupvijf, jvsupvijh, nblock)
       work2 = pijh*jvsupuijh
       work3 = pijh*jvsupvijh
       
       ALLOCATE(work4(0:mpol,-ntor:ntor,ns),                            &
                work5(0:mpol,-ntor:ntor,ns),                            &
                work6(0:mpol,-ntor:ntor,ns), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation2 failed in UPDATE_PRES'


!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
       CALL tomnsp(work2, work5, psin, ihalf)
       CALL tomnsp(work3, work6, psin, ihalf)
       CALL tomnsp(work1, work4, pcos, ihalf)
              
       DO m = 0, mpol
          DO n = -ntor, ntor
             work4(m,n,:) = work4(m,n,:)                               &
                          + m*work5(m,n,:) + n*nfp*work6(m,n,:)         !m: Pol. deriv of SIN(m*theta+n*NFP*phi)
          ENDDO                                                         !n: Tor. deriv of SIN(m*theta+n*NFP*phi)
       ENDDO

       djpmnch = -gamma*work4                                           !Gamma = adiabatic factor

#if defined(MRC)
       IF (l_force_asym) THEN
!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
          CALL tomnsp(work2, work5, fourier_parity_even, fourier_half_grid)
          CALL tomnsp(work3, work6, fourier_parity_even, fourier_half_grid)
          CALL tomnsp(work1, work4, fourier_parity_odd, fourier_half_grid)

          DO m = 0, mpol
             DO n = -ntor, ntor
                work4(m,n,:) = work4(m,n,:)                               &
                             - m*work5(m,n,:) - n*nfp*work6(m,n,:)         !m: Pol. deriv of COS(m*theta+n*NFP*phi)
             ENDDO                                                         !n: Tor. deriv of COS(m*theta+n*NFP*phi)
          ENDDO

          djpmnsh = -gamma*work4                                           !Gamma = adiabatic factor
       END IF
#endif

       DEALLOCATE(work5)
!
!      Now, calculate next term ~(gamma-1) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(v) = 0

       work1 = pijf_ds*jvsupsijf+pijf_du*jvsupuijf+pijf_dv*jvsupvijf
       CALL to_half_mesh(work1, work2, nblock)

!      MIXED HALF/FULL: DIAGONAL IN ANGULAR DERIVATIVES
!!       work1 = pijf_ds*jvsupsijf
!!       CALL to_half_mesh(work1,work2)
!!       work2 = work2+jvsupuijh*pijh_du+jvsupvijh*pijh_dv

       DEALLOCATE(jvsupuijh, jvsupvijh)
       
!Upwinding correction (makes 1st order streaming dominantly diagonal)
!NOTE: since it enters conservatively, it doesn't effect MHD force, only Hessian
!Can multiply this by a positive factor to increase numerical diffusion
       IF (l_upwind) THEN
          work1 = -0.5_dp*signjac*pijf_ds*ABS(jvsupsijf)
          work1(:,:,1) = 0;  work1(:,:,ns) = 0
          CALL GradientHalf(work3, work1, nblock)
          work2(:,:,2:ns) = work2(:,:,2:ns) + hs_i*work3(:,:,2:ns)
       END IF

       CALL tomnsp(work2, work4, pcos, ihalf)                           ! Calculate Fourier transform of WORK2


!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
       gamm1 = gamma-1
!!       IF (l_Compute_Hessian) gamm1 = gamm1 + xPosDef
       IF (l_linearize) gamm1 = gamm1 + xPosDef
  
       djpmnch = djpmnch + gamm1*work4                                  ! Harmonics of r.h.s. of Pres. Eq. (full mesh) 

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnch = delta_t*djpmnch                                        ! calculate increment of pressure harmonics.
       djpmnch(:,:,1)  = 0
       djpmnch(m0,:,1) = djpmnch(m0,:,2)                                ! store js=1, m=0 value

!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
!       djpmnch(m1:,:,2) = 0

#if defined(MRC)
       IF (l_force_asym) THEN
          CALL tomnsp(work2, work4, fourier_parity_odd, fourier_half_grid)                           ! Calculate Fourier transform of WORK2

          djpmnsh = djpmnsh + gamm1*work4                                  ! Harmonics of r.h.s. of Pres. Eq. (full mesh)

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
          djpmnsh = delta_t*djpmnsh                                        ! calculate increment of pressure harmonics.
          djpmnsh(:,:,1)  = 0
          djpmnsh(m0,:,1) = djpmnsh(m0,:,2)                                ! store js=1, m=0 value

!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
!       djpmnch(m1:,:,2) = 0
       END IF
#endif

       DEALLOCATE(work1, work2, work3, work4, work6)

       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)

       END SUBROUTINE update_pres

       END MODULE siesta_pressure
