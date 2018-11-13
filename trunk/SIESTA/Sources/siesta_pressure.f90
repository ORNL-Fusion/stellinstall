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
       USE v3_utilities, ONLY: assert
       USE stel_kinds
       USE stel_constants
       USE quantities
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE utilities, ONLY: GradientHalf, to_half_mesh
       IMPLICIT NONE

       PRIVATE
       INTEGER, PARAMETER  :: m0=0, m1=1, n0=0, isym=0, iasym=1
       INTEGER             :: nsmin, nsmax
       REAL(dp)            :: gamm1
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::                       &
         work1, work2, vgradph, work3, work4, work5, work6,             &
         jvsupuijh, jvsupvijh, pijf0_du, pijf0_dv
         
       PUBLIC :: update_pres 

       CONTAINS
!>     \brief Advance pressure on half radial grid from t to t+delta_t 
!>     by updating jpmnc,sh (jacobh*presh)
       SUBROUTINE update_pres
       USE timer_mod, ONLY: time_update_pres
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER             :: istat, i, j, k
       REAL(dp)            :: ton, toff
       LOGICAL, PARAMETER  :: l_upwind=.FALSE.                            !does not work well, keep it off
!-----------------------------------------------
       CALL second0(ton)
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
!=====================================
!      Compute first term ~(gamma-1)V * grad(p) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(V) = 0
!=====================================
       ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                        &
                work2(ntheta,nzeta,nsmin:nsmax),                        &   
                work3(ntheta,nzeta,nsmin:nsmax),                        &   
                vgradph(ntheta,nzeta,nsmin:nsmax),  stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation1 failed in UPDATE_PRES')

       ALLOCATE(work4(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work5(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation2 failed in UPDATE_PRES')

       ALLOCATE (jvsupuijh(ntheta,nzeta,nsmin:nsmax),                   &   
                 jvsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
       CALL ASSERT(istat.EQ.0,'Allocation3 failed in UPDATE_PRES')

!On exit, jsupXh valid on [nsmin+1:nsmax]
       CALL to_half_mesh(jvsupuijf, jvsupuijh)
       CALL to_half_mesh(jvsupvijf, jvsupvijh)

       gamm1 = gamma-1
!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
!      IF (l_linearize) gamm1 = gamma

!On exit, vgradph is valid on [nsmin+1:nsmax]
!MIXED HALF/FULL FORM: DIAGONAL IN ANGULAR DERIVATIVES - EXACTLY CANCEL 
!div(pV) TERM WHEN div(V) = 0

       work1 = pijf0_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
       CALL to_half_mesh(work1, vgradph)
       vgradph = vgradph                                                &
               + jvsupuijh(:,:,nsmin:nsmax)*pijh0_du(:,:,nsmin:nsmax)   &
               + jvsupvijh(:,:,nsmin:nsmax)*pijh0_dv(:,:,nsmin:nsmax)

!Compute second term on r.h.s. ~gamma*Div(p*V)
!Note that near axis, MUST set jvsupsijf == 0 since it makes no
!contribution to the Div term
       work2 = pijf0(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)


!"Upwinding" to suppress grid separation of pressure in radial direction
!makes diagonal contribution to 1st order streaming. Conservative so 
!makes no contribution to MHD force (integrates to zero). NOTE: this is
!not REAL upwinding, since the ABS(jvsupsijf) does not lead to a good
!derivative around velocity = 0!
       IF (l_upwind) work2 = work2 - hs_i*work1
       
       IF (nsmin .EQ. 1) work2(:,:,1) = 0                                 ! enforce bc at origin in div(pv) term
!On exit, work1 is valid on [nsmin+1:nsmax]
       CALL GradientHalf(work1, work2)

       nsmin=MAX(2,startglobrow)
       work2(:,:,nsmin:nsmax) = pijh0(:,:,nsmin:nsmax)*jvsupuijh(:,:,nsmin:nsmax)
       work3(:,:,nsmin:nsmax) = pijh0(:,:,nsmin:nsmax)*jvsupvijh(:,:,nsmin:nsmax)
       DEALLOCATE(jvsupuijh, jvsupvijh)

       CALL Continuity(djpmnch, isym)
       IF (lasym) CALL Continuity(djpmnsh, iasym)

       DEALLOCATE(work1, work2, work3, work4, work5, work6, vgradph)
       
       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)

       END SUBROUTINE update_pres
       

!>     \brief Evolve p, using density (n, continuity) equation, with p = n^gamma, 
       SUBROUTINE Continuity(djpmnh, ipar)
       USE shared_data, ONLY: delta_t
       USE fourier       
       USE island_params, ONLY: ns=>ns_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: ipar
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::  djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER   :: fours, fouru, fourv, m, n, np, i, j, k
       REAL(dp)  :: sparity
!-----------------------------------------------
!      djpmnh needed up to nsmax=endglobrow+1 to compute perturbed pressure
!      in siesta_forces

       IF (ipar .EQ. isym) THEN
          fours=f_cos; fouru=f_sin; fourv=f_sin
          sparity = 1
       ELSE
          fours=f_sin; fouru=f_cos; fourv=f_cos
          sparity = -1
       END IF

       m = LBOUND(djpmnh,3);  n = UBOUND(djpmnh,3)
       CALL ASSERT(nsmin.GE.m, nsmax.LE.n,                              &
                  'djpmnh: wrong bounds - Continuity')

!      FFT of gradp*Flux terms
       CALL tomnsp_par(vgradph, work4, fours)                            ! Calculate Fourier transform of v*gradp
       djpmnh(:,:,nsmin:nsmax) = gamm1*work4(:,:,nsmin:nsmax)            ! Harmonics of r.h.s. of pres eq. (half mesh) 

!FFT of -Div(p*Flux) terms 
       CALL tomnsp_par(work1(:,:,nsmin:nsmax), work4(:,:,nsmin:nsmax), fours)
       CALL tomnsp_par(work2(:,:,nsmin:nsmax), work5(:,:,nsmin:nsmax), fouru)
       CALL tomnsp_par(work3(:,:,nsmin:nsmax), work6(:,:,nsmin:nsmax), fourv)

!multiply angular derivatives by m ~d/du, n*nfp ~d/dv (*sparity below)
       DO m = 0, mpol
          work5(m,:,nsmin:nsmax) = m *work5(m,:,nsmin:nsmax)         
       END DO
       DO n = -ntor, ntor
          np = n*nfp
          work6(:,n,nsmin:nsmax) = np*work6(:,n,nsmin:nsmax)
       ENDDO
     
!sum radial and angular div terms     
       work4(:,:,nsmin:nsmax) = work4(:,:,nsmin:nsmax) +                &
      (work5(:,:,nsmin:nsmax) + work6(:,:,nsmin:nsmax))*sparity
       
!Gamma = adiabatic factor
       djpmnh(:,:,nsmin:nsmax) = djpmnh(:,:,nsmin:nsmax)                &
                               - gamma*work4(:,:,nsmin:nsmax)

!SPH Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnh(:,:,nsmin:nsmax) = delta_t*djpmnh(:,:,nsmin:nsmax)         ! calculate increment of pressure harmonics.

       IF (startglobrow .EQ. 1) THEN
          djpmnh(:,:,1)  = 0
          djpmnh(m0,:,1) = djpmnh(m0,:,2)                                ! store js=1, m=0 value
          CALL ASSERT(ALL(djpmnh(m0,-ntor:-1,:).EQ.zero),'djpmn(m0,-n) NONZERO')
       END IF
       
      
       END SUBROUTINE Continuity

       END MODULE siesta_pressure
