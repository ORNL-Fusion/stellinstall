!*******************************************************************************
!>  \brief Module contained subroutines for updating pressure as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 18, 2006
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
       USE nscalingtools, ONLY: startglobrow, endglobrow, ErrorAbort
       IMPLICIT NONE

       INTEGER, PARAMETER  :: m0=0, m1=1, isym=0, iasym=1
       INTEGER, PRIVATE    :: nsmin, nsmax, ns_span
       REAL(dp),PRIVATE    :: gamm1
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE ::               &
         work1, work2, vgradph, work3, work4, work5, work6, &
         jvsupuijh, jvsupvijh

       CONTAINS
!#if defined(SKS) 
!>     \brief Parallel subroutine to advance pressure from t to t+delta_t by updating jpmnc,sh (jacobh*presh)
       SUBROUTINE update_pres_par
       USE quantities
       USE timer_mod, ONLY: time_update_pres
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL, PARAMETER  :: l_upwind=.FALSE.
       INTEGER             :: istat, nblock
       REAL(dp)            :: ton, toff
!-----------------------------------------------
       CALL second0(ton)
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
       ns_span=nsmax-nsmin+1
!
!      First, calculate first term (~gamma) on r.h.s. of pressure evolution equation       
! 
       ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                         & ! Allocate work space
                work2(ntheta,nzeta,nsmin:nsmax),                         &   
                work3(ntheta,nzeta,nsmin:nsmax),                         &   
                vgradph(ntheta,nzeta,nsmin:nsmax),                       &   
                stat=istat)   
       IF (istat .NE. 0) CALL ErrorAbort('Allocation1 failed in UPDATE_PRES')
       nblock = SIZE(work1,1)*SIZE(work1,2)

       ALLOCATE(work4(0:mpol,-ntor:ntor,nsmin:nsmax),                    &
                work5(0:mpol,-ntor:ntor,nsmin:nsmax),                    &
                work6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) CALL ErrorAbort('Allocation2 failed in UPDATE_PRES')

       ALLOCATE (jvsupuijh(ntheta,nzeta,nsmin:nsmax),                    &   
                 jvsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) CALL ErrorAbort('Allocation3 failed in UPDATE_PRES')
       !!On exit, jsupXh valid on nsmin+1:nsmax grid
       CALL to_half_mesh_par(jvsupuijf(:,:,nsmin:nsmax), jvsupuijh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(jvsupvijf(:,:,nsmin:nsmax), jvsupvijh, nblock, nsmin, nsmax)

!=====================================
!   
!      Compute first term ~(gamma-1)V * grad(p) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(v) = 0
       gamm1 = gamma-1
!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
!      IF (l_linearize) gamm1 = gamm1+1

       IF (.TRUE.) THEN
          work1(:,:,nsmin:nsmax) =                                       &
          pijf_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax) +          &
          pijf_du(:,:,nsmin:nsmax)*jvsupuijf(:,:,nsmin:nsmax) +          &
          pijf_dv(:,:,nsmin:nsmax)*jvsupvijf(:,:,nsmin:nsmax)

       !!On exit, vgradph is valid on nsmin+1:nsmax grid
          CALL to_half_mesh_par(work1, vgradph, nblock, nsmin, nsmax)

       ELSE

!      MIXED HALF/FULL: DIAGONAL IN ANGULAR DERIVATIVES AND MORE ACCURATE AT EDGE js=ns
          work1 = pijf_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
          CALL to_half_mesh_par(work1, vgradph, nblock, nsmin, nsmax)
          vgradph(:,:,nsmin:nsmax) = vgradph(:,:,nsmin:nsmax)            &
               + jvsupuijh(:,:,nsmin:nsmax)*pijh_du(:,:,nsmin:nsmax)     &
               + jvsupvijh(:,:,nsmin:nsmax)*pijh_dv(:,:,nsmin:nsmax)
       END IF

!Upwinding correction (makes 1st order streaming dominantly diagonal)
!NOTE: since it enters conservatively, it doesn't effect MHD force, only HESSIAN
!Can multiply this by a positive factor to increase numerical diffusion
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*                      &
             pijf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijf(:,:,nsmin:nsmax))
          IF (nsmin .EQ. 1) work1(:,:,1) = 0
          CALL GradientHalf_par(work3, work1, nblock, nsmin, nsmax)
          vgradph(:,:,nsmin:nsmax) = vgradph(:,:,nsmin:nsmax) +          &
                              hs_i*work3(:,:,nsmin:nsmax)
       END IF

!=====================================
!      Compute second term on r.h.s. ~gamma*Div(p*V)
!      Note that near axis, jvsups ~ s cos(u+nv), jvsupv -> s sin(nv) due to jacobian element
!      In particular, pijf(1)*jvsupsijf makes no contribution to forces (integrates to zero,
!      since jvsupsijf(1) ~ cos(u) + O(s)!)
!      jvsupu may be finite, although dp/du -> 0 at axis
!      Note: WORK1 has COSINE parity, nsh = ns_i-1

!Upwind to suppress grid separation in radial direction
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*hs_i*                 &
             pijf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijf(:,:,nsmin:nsmax))
          work2(:,:,nsmin:nsmax) = work1(:,:,nsmin:nsmax) +              &
             pijh(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax) 
       ELSE
          work2(:,:,nsmin:nsmax) = pijf(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
       END IF

       !!On exit, work1 is valid on nsmin+1:nsmax grid
       CALL GradientHalf_par(work1, work2, nblock, nsmin, nsmax)

       nsmin=MAX(1,startglobrow)
       work2(:,:,nsmin:nsmax) = pijh(:,:,nsmin:nsmax)*jvsupuijh(:,:,nsmin:nsmax)
       work3(:,:,nsmin:nsmax) = pijh(:,:,nsmin:nsmax)*jvsupvijh(:,:,nsmin:nsmax)
       DEALLOCATE(jvsupuijh, jvsupvijh)

       CALL Continuity(djpmnch, isym)
       IF (lasym) CALL Continuity(djpmnsh, iasym)

       DEALLOCATE(work1, work2, vgradph, work3, work4, work5, work6)

       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)

       END SUBROUTINE update_pres_par
!#endif

       SUBROUTINE update_pres

       CALL update_pres_par

       END SUBROUTINE update_pres


       SUBROUTINE Continuity(djpmnh, iparity)
       USE shared_data, ONLY: delta_t
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: iparity
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
            djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER   :: fours, fouru, fourv, m, n, mp, np, MPI_ERR
       REAL(dp)  :: sparity
!-----------------------------------------------

       IF (iparity .EQ. isym) THEN
          fours=f_cos; fouru=f_sin; fourv=f_sin
          sparity = 1
       ELSE
          fours=f_sin; fouru=f_cos; fourv=f_cos
          sparity = -1
       END IF

       m = LBOUND(djpmnh,3);  n = UBOUND(djpmnh,3)
       IF (m.NE.nsmin .OR. n.NE.nsmax) CALL ErrorAbort('djpmnh wrong bounds - Continuity')

!      FFT of gradp*Flux terms
       CALL tomnsp_par(vgradph, work4, fours, 1, ns_span)                ! Calculate Fourier transform of v_dot_gradph
       djpmnh = gamm1*work4(:,:,nsmin:nsmax)                             ! Harmonics of r.h.s. of Pres. Eq. (full mesh) 

!      FFT of Div(p*Flux) terms (multiply angular derivatives by m,n*nfp)
       CALL tomnsp_par(work1, work4, fours, 1, ns_span)
       CALL tomnsp_par(work2, work5, fouru, 1, ns_span)
       CALL tomnsp_par(work3, work6, fourv, 1, ns_span)

       DO m = 0, mpol
          mp = m*sparity
          DO n = -ntor, ntor
             np = n*nfp*sparity
             work5(m,n,nsmin:nsmax) = mp*work5(m,n,nsmin:nsmax)          &
                                    + np*work6(m,n,nsmin:nsmax)          ! mp: d/du; np: d/dv
          ENDDO
       ENDDO
     
       work4 = work4 + work5

       djpmnh = djpmnh - gamma*work4(:,:,nsmin:nsmax)                    ! Gamma = adiabatic factor

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnh = delta_t*djpmnh                                           ! calculate increment of pressure harmonics.

!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics

       IF (nsmin .EQ. 1) THEN
          djpmnh(:,:,1)  = 0
          djpmnh(m0,:,1) = djpmnh(m0,:,2)                                ! store js=1, m=0 value
       END IF

       END SUBROUTINE Continuity


       END MODULE siesta_pressure
