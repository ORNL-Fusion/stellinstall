!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref siesta_pressure
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Updates pressure.
!*******************************************************************************

      MODULE siesta_pressure
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE quantities
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE utilities, ONLY: GradientHalf, to_half_mesh

      IMPLICIT NONE

!*******************************************************************************
!  utilities module variables
!*******************************************************************************
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workmn4
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workmn5
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workmn6

!*******************************************************************************
!  utilities module parameters
!*******************************************************************************
      PRIVATE !  Fix me refernce these from a central location.
      INTEGER, PARAMETER                       :: m0 = 0
      INTEGER, PARAMETER                       :: m1 = 1
      INTEGER, PARAMETER                       :: n0 = 0
      INTEGER, PARAMETER                       :: isym = 0
      INTEGER, PARAMETER                       :: iasym = 1

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PUBLIC :: update_pres

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief advance pressure on the half radial grid from t to t + delta_t.
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
!>     dp/dt = (gamma - 1)*(v^s*dp/ds + v^u*dp/du + v^v*dp/dv)               (7)
!>           - gamma/J*(dJpv^s/ds + dJpv^u/du + dJpv^v/dv)
!-------------------------------------------------------------------------------
      SUBROUTINE update_pres
      USE timer_mod, ONLY: time_update_pres

      IMPLICIT NONE

!  local variables
      INTEGER                                  :: istat
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jvsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jvsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij1
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij2
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij3
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: vgradph
      INTEGER                                  :: nsmin
      INTEGER                                  :: nsmax

!  local parameters
!  Does not work well, keep it off.
      LOGICAL, PARAMETER                       :: l_upwind = .FALSE.

!  Start of executable code
      CALL second0(ton)
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 1)

!  Compute first term ~(gamma-1)V * grad(p) on r.h.s. of pressure equation the
!  u,v derivatives should be consistent with div(pv) term and combine to cancel
!  if div(V) = 0
      ALLOCATE(workij1(ntheta,nzeta,nsmin:nsmax),                              &
               workij2(ntheta,nzeta,nsmin:nsmax),                              &
               workij3(ntheta,nzeta,nsmin:nsmax),                              &
               vgradph(ntheta,nzeta,nsmin:nsmax),  stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in update_pres')

      ALLOCATE(workmn4(0:mpol,-ntor:ntor,nsmin:nsmax),                         &
               workmn5(0:mpol,-ntor:ntor,nsmin:nsmax),                         &
               workmn6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in update_pres')

      ALLOCATE (jvsupuijh(ntheta,nzeta,nsmin:nsmax),                           &
                jvsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation3 failed in update_pres')

!  On exit, jsupXh valid on [nsmin+1:nsmax]
      CALL to_half_mesh(jvsupuijf, jvsupuijh)
      CALL to_half_mesh(jvsupvijf, jvsupvijh)

!  On exit, vgradph is valid on [nsmin+1:nsmax]. Mixed half/full form: Diagonal
!  in angular derivatives. Exactly cancel Div(pV) term when div(V) = 0.
      workij1 = pijf0_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
      CALL to_half_mesh(workij1, vgradph)
      vgradph = vgradph                                                        &
              + jvsupuijh(:,:,nsmin:nsmax)*pijh0_du(:,:,nsmin:nsmax)           &
              + jvsupvijh(:,:,nsmin:nsmax)*pijh0_dv(:,:,nsmin:nsmax)

!  Compute second term on r.h.s. ~gamma*Div(p*V). Note that near axis, must set
!  jvsupsijf == 0 since it makes no contribution to the Div term. Compute the
!  pv term.
      workij2 = pijf0(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)

!  "Upwinding" to suppress grid separation of pressure in radial direction makes
!  diagonal contribution to 1st order streaming. Conservative so makes no
!  contribution to MHD force (integrates to zero). Note this is not real
!  upwinding, since the ABS(jvsupsijf) does not lead to a good derivative around
!  velocity = 0!
      IF (l_upwind) THEN
         workij2 = workij2 - hs_i*workij1
      END IF

      IF (nsmin .EQ. 1) THEN
!  Enforce bc at origin in Div(pv) term.
         workij2(:,:,1) = 0
      END IF
!  Compute d/ds p*jv where j is the jacobian, v is the displacement vector and
!  p is the pressure. On exit, work1 is valid on [nsmin + 1:nsmax].
      CALL GradientHalf(workij1, workij2)

!  Work 2 and 3 now represent half mesh quantities.
      nsmin = MAX(2, startglobrow)
      workij2(:,:,nsmin:nsmax) = pijh0(:,:,nsmin:nsmax)*jvsupuijh(:,:,nsmin:nsmax)
      workij3(:,:,nsmin:nsmax) = pijh0(:,:,nsmin:nsmax)*jvsupvijh(:,:,nsmin:nsmax)
      DEALLOCATE(jvsupuijh, jvsupvijh)

      CALL Continuity(djpmnch, vgradph, workij1, workij2, workij3, isym,       &
                      nsmin, nsmax)
      IF (lasym) THEN
         CALL Continuity(djpmnsh, vgradph, workij1, workij2, workij3, iasym,   &
                         nsmin, nsmax)
      END IF

      DEALLOCATE(workij1, workij2, workij3, workmn4, workmn5, workmn6, vgradph)
       
      CALL second0(toff)
      time_update_pres = time_update_pres + (toff - ton)

      END SUBROUTINE update_pres
       
!-------------------------------------------------------------------------------
!>  @brief Evolve p, using density (n, continuity) equation, with p = n^gamma.
!>
!>  djpmnh needed up to nsmax = endglobrow + 1 to compute perturbed pressure in
!>  @ref siesta_forces.
!>
!>  @param[inout] djpmnh
!>  @param[in]    vgradph     Diplacement time gradient of pressure on half
!>                            grid.
!>  @param[in]    dpjsupsdsij Real space d/ds p*jvsups.
!>  @param[in]    pjsupuij    Real space p*jvsupu.
!>  @param[in]    pjsupvij    Real space p*jvsupu.
!>  @param[in]    iparity     Parity flag of fouier components.
!>  @param[in]    nsmin       Minimum s index.
!>  @param[in]    nsmax       Maximum s index.
!-------------------------------------------------------------------------------
      SUBROUTINE Continuity(djpmnh, vgradph, dpjsupsdsij, pjsupuij, pjsupvij,  &
                            iparity, nsmin, nsmax)
      USE shared_data, ONLY: delta_t
      USE fourier
      USE island_params, ONLY: ns=>ns_i

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: djpmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vgradph
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: dpjsupsdsij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: pjsupuij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: pjsupvij
      INTEGER, INTENT(in)                                     :: iparity
      INTEGER, INTENT(in)                                     :: nsmin
      INTEGER, INTENT(in)                                     :: nsmax

!  local variables
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      REAL (dp)                                               :: sparity

!  Start of executable code
      IF (iparity .EQ. isym) THEN
         fours = f_cos
         fouruv = f_sin
         sparity = 1
      ELSE
         fours = f_sin
         fouruv = f_cos
         sparity = -1
      END IF

      m = LBOUND(djpmnh,3)
      n = UBOUND(djpmnh,3)
      CALL assert(nsmin .ge. m, nsmax .le. n,                                  &
                  'djpmnh: wrong bounds - Continuity')

!  Calculate Fourier transform of v*gradp, r.h.s. of pres eq. (half mesh).
      CALL tomnsp(vgradph, workmn4, fours)
      djpmnh(:,:,nsmin:nsmax) = (gamma - 1)*workmn4(:,:,nsmin:nsmax)

!  FFT of -Div(p*Flux) terms
      CALL tomnsp(dpjsupsdsij(:,:,nsmin:nsmax), workmn4(:,:,nsmin:nsmax),      &
                  fours)
      CALL tomnsp(pjsupuij(:,:,nsmin:nsmax), workmn5(:,:,nsmin:nsmax), fouruv)
      CALL tomnsp(pjsupvij(:,:,nsmin:nsmax), workmn6(:,:,nsmin:nsmax), fouruv)

!  Multiply angular derivatives by m ~d/du, n*nfp ~d/dv (*sparity below).
      DO m = 0, mpol
         workmn5(m,:,nsmin:nsmax) = m*workmn5(m,:,nsmin:nsmax)
      END DO
      DO n = -ntor, ntor
         workmn6(:,n,nsmin:nsmax) = n*nfp*workmn6(:,n,nsmin:nsmax)
      ENDDO
     
!  Sum radial and angular div terms
      workmn4(:,:,nsmin:nsmax) = workmn4(:,:,nsmin:nsmax)                      &
                               + (workmn5(:,:,nsmin:nsmax) +                   &
                                  workmn6(:,:,nsmin:nsmax))*sparity
       
!  Gamma = adiabatic factor
      djpmnh(:,:,nsmin:nsmax) = djpmnh(:,:,nsmin:nsmax)                        &
                              - gamma*workmn4(:,:,nsmin:nsmax)

!  Evolve harmonics of (jacobh * pres)mn to conserve pressure. Calculate
!  increment of pressure harmonics.
      djpmnh(:,:,nsmin:nsmax) = delta_t*djpmnh(:,:,nsmin:nsmax)

!  Store js=1, m=0 value.
      IF (startglobrow .EQ. 1) THEN
         djpmnh(:,:,1)  = 0
         djpmnh(m0,:,1) = djpmnh(m0,:,2)! store js=1, m=0 value
         CALL assert(ALL(djpmnh(m0,-ntor:-1,:) .eq. zero),                     &
                         'djpmn(m0,-n) NONZERO')
      END IF

      END SUBROUTINE

      END MODULE
