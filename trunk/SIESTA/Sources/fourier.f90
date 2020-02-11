!*******************************************************************************
!>  @file fourier.f90
!>  @brief Contains module @ref fourier
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module contains subroutines for computing FFT on parallel radial intervals.
!>  Converts quantities between real and fourier space. Note @ref fixarray must
!>  be called once before call any of the Fourier subroutines to calculate all
!>  necessary sine and cosine factors on a fixed mesh of collocation angles.
!*******************************************************************************
      MODULE fourier
      USE v3_utilities, ONLY: assert, assert_eq
      USE stel_kinds
      USE stel_constants
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,                  &
                         mpol=>mpol_i, ntor=>ntor_i, nfp=>nfp_i
      USE timer_mod
      USE descriptor_mod, ONLY: INHESSIAN, DIAGONALDONE
      USE shared_data, ONLY: unit_out, lasym
      
      IMPLICIT NONE

!*******************************************************************************
!  fourier module variables
!*******************************************************************************
!>  Orthnormalization factor.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE :: orthonorm

!*******************************************************************************
!  fourier module parameters
!*******************************************************************************
!>  Cosine parity.
      INTEGER, PARAMETER   :: f_cos = 0
!>  Sine parity.
      INTEGER, PARAMETER   :: f_sin = 1

!  Bit locations for fourier control flags.
!>  Do not sum fouier real space transformed quantity. This is used when a
!>  stellarator symmetric parity is being converted to real space.
      INTEGER, PARAMETER   :: f_none = 0               !0000
!>  Radial derivative flag.
      INTEGER, PARAMETER   :: f_ds   = 1               !0001
!>  Poloidal derivative flag.
      INTEGER, PARAMETER   :: f_du   = 2               !0010
!>  Toroidal derivative flag.
      INTEGER, PARAMETER   :: f_dv   = 4               !0100
!>  Combined toroidal and poloidal derivatives.
      INTEGER, PARAMETER   :: f_dudv = IOR(f_du, f_dv) !0110
!>  Sum fouier real space transformed quantity. This is used when a
!>  non-stellarator symmetric parity is being converted to real space.
      INTEGER, PARAMETER   :: f_sum  = 8               !1000

!  Test bit positions for Derivative flags.
!>  Bit position of the @ref f_ds flag.
      INTEGER, PARAMETER   :: b_ds   = 0
!>  Bit position of the @ref f_du flag.
      INTEGER, PARAMETER   :: b_du   = 1
!>  Bit position of the @ref f_dv flag.
      INTEGER, PARAMETER   :: b_dv   = 2
!>  Bit position of the @ref f_sum flag.
      INTEGER, PARAMETER   :: b_sum  = 3

!>  Fouier test accuracy tolarance.
      REAL (dp), PARAMETER :: four_tol = 1.E-12_dp

!  Special Fourier modes.
!  Poloidal modes.
!>  m = 0 mode.
      INTEGER, PARAMETER   :: m0 = 0
!>  m = 1 mode.
      INTEGER, PARAMETER   :: m1 = 1
!>  m = 2 mode.
      INTEGER, PARAMETER   :: m2 = 2
!  Toroidal modes.
!>  n = 0 mode.
      INTEGER, PARAMETER   :: n0 = 0
!>  n = 1 mode.
      INTEGER, PARAMETER   :: n1 = 1

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Convert a 2D quantity from fourier to real space.
!>
!>  @param[in]    xmn    Fourier space quantity.
!>  @param[inout] xuv    Realspace quantity.
!>  @param[in]   dflag   Control flag for angle derivatives and summing with
!>                       previous quantity.
!>  @param[in]   iparity Fourier parity for the conversion. Either @ref isym or
!>                       @ref iasym.
!-------------------------------------------------------------------------------
      SUBROUTINE vtoijsp(xmn, xuv, dflag, iparity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), INTENT(in)    :: xmn
      REAL (dp), DIMENSION(:,:), INTENT(inout) :: xuv
      INTEGER, INTENT(in)                      :: dflag
      INTEGER, INTENT(in)                      :: iparity

!  local variables
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: xmn_3
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: xuv_3
      INTEGER                                  :: mpol_s
      INTEGER                                  :: ntor_s
      INTEGER                                  :: n1
      INTEGER                                  :: n2
      INTEGER                                  :: i1
      INTEGER                                  :: i2

!  Start of executable code
!  Save the number of poloial and toroidal modes.
      mpol_s = mpol
      ntor_s = ntor
      mpol = mpol_v
      ntor = ntor_v

!  Reshape the arrays to 3 dimensions so the regular toij routine can be used.
      n1 = SIZE(xmn,1)
      n2 = SIZE(xmn,2)
	  i1 = SIZE(xuv,1)
      i2 = SIZE(xuv,2)

      ALLOCATE(xmn_3(n1,n2,1), xuv_3(i1,i2,1))

	  xmn_3 = RESHAPE(xmn, SHAPE(xmn_3))
	  xuv_3 = RESHAPE(xuv, SHAPE(xuv_3))

	  CALL toijsp(xmn_3, xuv_3, dflag, iparity)

	  xuv = RESHAPE(xuv_3, SHAPE(xuv))

      DEALLOCATE(xmn_3, xuv_3)

!  Restore the poloidal and toroidal mode size.
      mpol = mpol_s
      ntor = ntor_s

	  END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert a quantity from fourier to real space in parallel.
!>
!>  This subroutine moves a quantity to real space by summing over its Fourier
!>  harmonics. If the @ref f_du flag is set take the poloidal derivative. If the
!>  @ref f_dv flag is set take the toroidal derivative. If the f_sum flag is set
!>  add the real quanity to the previous value.
!>
!>  @param[in]    xmn     Fourier space quantity.
!>  @param[inout] xuv     Real space quantity.
!>  @param[in]    dflag   Derivative and sum control flags.
!>  @param[in]    iparity Fourier parity flag.
!-------------------------------------------------------------------------------
      SUBROUTINE toijsp(xmn, xuv, dflag, iparity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: xmn
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: xuv
      INTEGER, INTENT(in)                        :: dflag
      INTEGER, INTENT(in)                        :: iparity

!  local variables
      INTEGER                                    :: m
      INTEGER                                    :: n
      INTEGER                                    :: jk
      INTEGER                                    :: lk
      INTEGER                                    :: js
      REAL (dp), DIMENSION(:,:), ALLOCATABLE     :: work0
      REAL (dp), DIMENSION(:,:), ALLOCATABLE     :: work1
      REAL (dp), DIMENSION(:,:), ALLOCATABLE     :: work2
      REAL (dp), DIMENSION(:,:), ALLOCATABLE     :: xuv1
      REAL (dp)                                  :: skston
      REAL (dp)                                  :: skstoff

!  Start of executable code
      CALL second0(skston)

      ALLOCATE(work0(0:mpol,-ntor:ntor), xuv1(ntheta,nzeta),                   &
               work1(ntheta,-ntor:ntor), work2(ntheta,-ntor:ntor),             &
               stat=lk)
      CALL assert_eq(0, lk, 'Allocation error in TOIJSP')

      NS_LOOP: DO js = 1, MIN(SIZE(xmn,3), SIZE(xuv,3))
         work0(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
         work1 = zero
         work2 = zero

         DO jk = 1, ntheta
!  First sum over poloidal modes
            DO m = 0, mpol
!  Poloidal derivative requested use COSMUM, SINMUM
               IF (BTEST(dflag, b_du)) THEN
                  work1(jk,:) = work1(jk,:) - work0(m,:)*sinmum(jk, m)
                  work2(jk,:) = work2(jk,:) + work0(m,:)*cosmum(jk, m)
               ELSE
                  work1(jk,:) = work1(jk,:) + work0(m,:)*cosmu(jk, m)
                  work2(jk,:) = work2(jk,:) + work0(m,:)*sinmu(jk, m)
               END IF
            END DO
         END DO
!  Make sure that spatial variable is zeroed.
         xuv1 = zero

!  First do the n = 0 mode.
         DO lk = 1, nzeta
            IF (iparity .eq. f_cos) THEN
               IF (BTEST(dflag, b_dv)) THEN
                  xuv1(:,lk) = zero
               ELSE
                  xuv1(:,lk) = work1(:,0)
               END IF
            ELSE
               IF (BTEST(dflag, b_dv)) THEN
                  xuv1(:,lk) = zero
               ELSE
                  xuv1(:,lk) = work2(:,0)
               END IF
            END IF

!  Sum over remaining toroidal modes.
            DO n = 1, ntor
               IF (iparity .eq. f_cos) THEN
                  IF (BTEST(dflag, b_dv)) THEN
                     xuv1(:,lk) =  xuv1(:,lk)                                  &
                                - (work1(:,n) + work1(:,-n))*sinnvn(lk, n)     &
                                - (work2(:,n) - work2(:,-n))*cosnvn(lk, n)
                  ELSE
                     xuv1(:,lk) =  xuv1(:,lk)                                  &
                                + (work1(:,n) + work1(:,-n))*cosnv(lk, n)      &
                                - (work2(:,n) - work2(:,-n))*sinnv(lk, n)
                  END IF
               ELSE !  Sine series.
                  IF (BTEST(dflag, b_dv)) THEN
                     xuv1(:,lk) =  xuv1(:,lk)                                  &
                                + (work1(:,n) - work1(:,-n))*cosnvn(lk, n)     &
                                - (work2(:,n) + work2(:,-n))*sinnvn(lk, n)
                  ELSE
                     xuv1(:,lk) =  xuv1(:,lk)                                  &
                                + (work1(:,n) - work1(:,-n))*sinnv(lk, n)      &
                                + (work2(:,n) + work2(:,-n))*cosnv(lk, n)
                  END IF
               END IF
            END DO
         END DO

         IF (BTEST(dflag, b_sum)) THEN
            xuv(:,:,js) = xuv(:,:,js) + xuv1(:,:)
         ELSE
            xuv(:,:,js) = xuv1(:,:)
         END IF
      END DO NS_LOOP

      DEALLOCATE(work0, work1, work2, xuv1)

      CALL second0(skstoff)
      IF (DIAGONALDONE .AND. INHESSIAN) THEN
         toijsp_time = toijsp_time + (skstoff - skston)
      END IF

      time_toijsp = time_toijsp + (skstoff - skston)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert a quantity from real to fourier space.
!>
!>  This subroutine moves a quantity to Fourier space by summing over its real
!>  space positions.
!>
!>  @param[in]  xuv     Real space quantity.
!>  @param[out] xmn     Fourier space quantity.
!>  @param[in]  iparity Fourier parity flag.
!-------------------------------------------------------------------------------
      SUBROUTINE tomnsp(xuv, xmn, iparity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:,:), INTENT(out) :: xmn
      INTEGER, INTENT(in)                      :: iparity

!  local variables
      INTEGER                                  :: m
      INTEGER                                  :: n
      INTEGER                                  :: jk
      INTEGER                                  :: lk
      INTEGER                                  :: js
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: work1
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: work2
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: xmn1
      REAL (dp), DIMENSION(:), ALLOCATABLE     :: work3
      REAL (dp), DIMENSION(:), ALLOCATABLE     :: work4
      REAL (dp)                                :: skston
      REAL (dp)                                :: skstoff
      INTEGER                                  :: istat

!  Start of executable code
      CALL second0(skston)
      ALLOCATE(work1(0:mpol,nzeta), work2(0:mpol,nzeta),                       &
               work3(0:mpol), work4(0:mpol),                                   &
               xmn1(0:mpol,-ntor:ntor), stat=istat)
      CALL assert_eq(0, istat, 'Allocation error in tomnsp')

      NS_LOOP: DO js = 1, MIN(UBOUND(xmn,3), UBOUND(xuv,3))
         xmn1      = 0

!  First, add over poloidal collocation angles. Note use of cosmuI/sinmuI
!  include normalization factors.
         DO lk = 1, nzeta
            DO m = 0, mpol
               work1(m,:) = 0
               work2(m,:) = 0
               DO jk = 1, ntheta
                  work1(m,lk) = work1(m,lk) + xuv(jk,lk,js)*cosmui(jk, m)
                  work2(m,lk) = work2(m,lk) + xuv(jk,lk,js)*sinmui(jk, m)
               END DO
            END DO
         END DO

!  Then add over toroidal collocation angles
         DO n = 0, ntor
            DO lk = 1, nzeta
               IF (iparity .eq. f_cos) THEN
                  work3(:) = work1(:,lk)*cosnv(lk, n)
                  work4(:) = work2(:,lk)*sinnv(lk, n)
                  xmn1(:,n) =  xmn1(:,n) + work3(:) - work4(:)
                  IF (n .ne. 0) THEN
                     xmn1(1:,-n) = xmn1(1:,-n) + work3(1:) + work4(1:)
                  END IF
               ELSE ! Sine series
                  work3(:) = work1(:,lk)*sinnv(lk, n)
                  work4(:) = work2(:,lk)*cosnv(lk, n)
                  xmn1(:,n) =  xmn1(:,n) + work3(:) + work4(:)
                  IF (n .ne. 0) THEN
                     xmn1(1:,-n) = xmn1(1:,-n) - work3(1:) + work4(1:)
                  END IF
               END IF
            END DO
         END DO

! Redundant: To avoid counting (0,-n) for non-zero n.
         IF (iparity .eq. f_cos) THEN
            xmn1(m0,-ntor:-n1) = zero
         ELSE
            xmn1(m0,-ntor:n0) = zero
         END IF
         xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn1(:,:)
      END DO NS_LOOP

      DEALLOCATE(work1, work2, work3, work4, xmn1)

      CALL second0(skstoff)
      IF (DIAGONALDONE.AND.INHESSIAN) THEN
         tomnsp_time = tomnsp_time + (skstoff - skston)
      END IF
      time_tomnsp = time_tomnsp + (skstoff - skston)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Computes the origin value of a half-mesh quantity.
!>
!>  The origin quantity is found from the m mode of the js = 2 value.
!>
!>  @param[inout] xuv  Real space quantity.
!>  @param[in]    mode Poloidal mode number to be retained.
!-------------------------------------------------------------------------------
      SUBROUTINE getorigin(xuv, mode)

      IMPLICIT NONE

!  Declare Arguments
       REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: xuv
       INTEGER, INTENT(in)                                    :: mode

!  local variables
       REAL(dp), DIMENSION(0:mpol,-ntor:ntor, 2)              :: xmn

!  local parameters
       INTEGER, PARAMETER                                     :: js1 = 1
       INTEGER, PARAMETER                                     :: js2 = 2

!  Start of executable code
       CALL tomnsp2(xuv, xmn(:,:,1), f_cos, js2)
       IF (lasym) THEN
          CALL tomnsp2(xuv, xmn(:,:,2), f_sin, js2)
       END IF

       IF (mode .GT. 0) THEN
          xmn(0:mode-1,:,1:2) = 0
       END IF

       IF (mode .LT. mpol) THEN
          xmn(mode+1:,:,1:2) = 0
       END IF

       CALL toijsp2(xmn(:,:,1), xuv, f_cos, f_none, js1)
       IF (lasym) THEN
          CALL toijsp2(xmn(:,:,2), xuv, f_sin, f_sum, js1)
       END IF

       END SUBROUTINE

!SPH 10/15/14
!-------------------------------------------------------------------------------
!>  @brief  Convert to Fourier pest coordinates.
!>
!>  This subroutine moves a quantity X to Fourier space producing its harmonics
!>  X_MN by summing over all the presribed (toroidal and poloidal) collocation
!>  points:
!>
!>    theta_j = j*pi/M,        j = 0,..., M
!>    zeta_k  = k*2*pi/(2N+1), k = 0,..., 2N
!>
!>  where M = mpol - 1 and N = ntor. XUV is in VMEC coordinates. the XMN are
!>  harmonics in PEST coordinates, where
!>
!>    UP(s,u,v) = U + LAMBDA(s,u,v)
!>
!>  @param[in]  xuv  Real space quantity.
!>  @param[out] xmn  Fourier space quantity.
!>  @param[in]  lmns Lambda for even parity.
!>  @param[in]  lmnc Lambda for odd parity.
!>  @param[in]  ipar Partity flag.
!>                   * 0 Even parity
!>                   * 1 Odd parity.
!-------------------------------------------------------------------------------
      SUBROUTINE tomnsp_pest(xuv, xmn, lmns, lmnc, IPAR)
      USE stel_kinds

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(ntheta,nzeta,ns), INTENT(IN)       :: xuv
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT) :: xmn
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: lmns
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: lmnc
      INTEGER, INTENT(IN)                                     :: ipar

!  Local variables
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: jk
      INTEGER                                                 :: istat
      INTEGER                                                 :: ns_save
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: work1
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: work2
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: work3
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: work4
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: lambda
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: dupdu
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: cos1l
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: sin1l
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: cosml
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: sinml
      REAL (dp), ALLOCATABLE, DIMENSION(:,:)                  :: cosmuip
      REAL (dp), ALLOCATABLE, DIMENSION(:,:)                  :: sinmuip
      REAL (dp)                                               :: ton
      REAL (dp)                                               :: toff

!  Start of executable code
      CALL second0(ton)
      ns_save = ns
      ns = MIN(SIZE(lmns, 3), SIZE(xmn,3)) ! only up to s=1, NOT vacuum

!  Compute lambda
      ALLOCATE (lambda(ntheta,nzeta,ns), dupdu(ntheta,nzeta,ns),               &
                cos1l(ntheta,nzeta,ns), cosml(ntheta,nzeta,ns),                &
                sin1l(ntheta,nzeta,ns), sinml(ntheta,nzeta,ns),                &
                stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error #1 in TOMNSP_PEST')
        
      CALL TOIJSP(lmns, lambda, f_none, f_sin)
      CALL TOIJSP(lmns, dupdu,  f_du, f_sin)
      IF (lasym) THEN
         CALL TOIJSP(lmnc, lambda, f_sum, f_cos)
         CALL TOIJSP(lmnc, dupdu,  IOR(f_sum,f_du), f_cos)
      END IF

!  Transform from up to u
      dupdu = 1 + dupdu                                              

!       INTERPOLATE TO FULL GRID IF ihalf=0
!        IF (ihalf .EQ. 0) THEN
           lambda(:,:,1) = lambda(:,:,2)                                  !could do better, why bother?
           dupdu (:,:,1) = dupdu (:,:,2)
           DO jk=2,ns-1
              lambda(:,:,jk) = .5_dp*(lambda(:,:,jk)+lambda(:,:,jk+1))
              dupdu (:,:,jk) = .5_dp*(dupdu (:,:,jk)+dupdu (:,:,jk+1))
           END DO
           lambda(:,:,ns) = 2*lambda(:,:,ns)-lambda(:,:,ns-1)           
           dupdu (:,:,ns) = 2*dupdu (:,:,ns)-dupdu (:,:,ns-1)           
!        END IF

        cos1l = COS(lambda); sin1l = SIN(lambda)
!
!       2. DO Fourier Transform in PEST UP=u+lambda, VP=v Coordinates
!          Use recurrence relations below for EFFICIENCY
!
!          cosmuiP==cos(mu_p) = cosmui(u,m)*cos(m*lambda) - sinmui(u,m)*sin(m*lambda)
!          sinmuiP==sin(mu_p) = sinmui(u,m)*cos(m*lambda) + cosmui(u,m)*sin(m*lambda)
!          Let cos(m*lambda) == cosml(s,u,v,m)
!              sin(m*lambda) == sinml(s,u,v,m)
!          use recurrence formulae to compute these efficiently
!              cosml(m=0) = 1;  sinml(m=0) = 0
!              compute cosml(m=1), sinml(m=1)
!          then
!              cosml(m+1) = cosml(m)*cosml(1) - sinml(m)*sinml(1)
!              sinml(m+1) = sinml(m)*cosml(1) + cosml(m)*sinml(1)
!
        ALLOCATE(work1(0:mpol,nzeta,ns), work2(0:mpol,nzeta,ns),        &
                 work3(ntheta,nzeta,ns), cosmuiP(nzeta,ns),             &
                 sinmuiP(nzeta,ns), stat=istat)
        CALL ASSERT(istat.eq.0,'Allocation error #2 in TOMNSP_PEST')
        work1 = zero; work2 = zero

        DO m = 0, mpol
           IF (m .EQ. 0) THEN
              cosml = dupdu;  sinml=0
           ELSE IF (m .EQ. 1) THEN
              cosml = cos1l*dupdu; sinml = sin1l*dupdu
           ELSE
              work3 = cosml
              cosml = cosml*cos1l - sinml*sin1l
              sinml = sinml*cos1l + work3*sin1l
           END IF
           DO jk = 1, ntheta                                             ! First, add over poloidal collocation angles
              cosmuiP(:,1:ns) = cosmui(jk,m)*cosml(jk,:,1:ns)           &
                      - sinmui(jk,m)*sinml(jk,:,1:ns)
              sinmuiP(:,1:ns) = sinmui(jk,m)*cosml(jk,:,1:ns)           &
                      + cosmui(jk,m)*sinml(jk,:,1:ns)
              work1(m,:,1:ns) = work1(m,:,1:ns)                         &
                               + xuv(jk,:,1:ns)*cosmuiP(:,1:ns)          ! Note use of cosmuI/sinmuI; include norm. factors
              work2(m,:,1:ns) = work2(m,:,1:ns)                         &
                               + xuv(jk,:,1:ns)*sinmuiP(:,1:ns)
           ENDDO
        END DO

        DEALLOCATE (work3, cosmuip, sinmuip)
        xmn = zero                                                       ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta,ns), work4(0:mpol,nzeta, ns),       &
                 stat=istat)
        CALL ASSERT(istat.eq.0,'Allocation error #3 in TOMNSP_PEST')
        work3 = zero; work4 = zero

        DO n = 0, ntor
          DO jk = 1, nzeta                                               ! Then, add over toroidal collocation angles
            IF (ipar == f_cos) THEN                                      ! COSINE series
              work3(:,jk,1:ns) = work1(:,jk,1:ns)*cosnv(jk, n)
              work4(:,jk,1:ns) = work2(:,jk,1:ns)*sinnv(jk, n)
              xmn(:,n,1:ns) =  xmn(:,n,1:ns)                            &
                 + work3(:,jk,1:ns) - work4(:,jk,1:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,1:ns) = xmn(1:,-n,1:ns)                       &
                 + work3(1:,jk,1:ns) + work4(1:,jk,1:ns)
              ENDIF
            ELSE                                                         ! SINE series
              work3(:,jk,1:ns) = work1(:,jk,1:ns)*sinnv(jk, n)
              work4(:,jk,1:ns) = work2(:,jk,1:ns)*cosnv(jk, n)
              xmn(:,n,1:ns) =  xmn(:,n,1:ns)                            &
                + work3(:,jk,1:ns) + work4(:,jk,1:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,1:ns) = xmn(1:,-n,1:ns)                       &
                - work3(1:,jk,1:ns) + work4(1:,jk,1:ns)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(0,-ntor:-1,:) = zero                                         ! Redundant: To avoid counting (0,-n) for non-zero n.

        DO jk = 1,ns
           xmn(:,:,jk) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,jk)
        END DO
         
        DEALLOCATE(work1, work2, work3, work4, lambda, cos1l, sin1l,    &
                   sinml, cosml, dupdu)

        CALL second0(toff)
        IF (DIAGONALDONE.AND.INHESSIAN) THEN
           tomnsp_time = tomnsp_time + (toff - ton)
        END IF

        time_tomnsp = time_tomnsp + (toff - ton)
        
        ns = ns_save

        END SUBROUTINE

        SUBROUTINE TOMNSP2(XUV, XMN, IPAR, JS)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics XMN by
!         summing over all the prescribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPAR = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!       JS:  Surface where xmn is computed
!
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), DIMENSION(ntheta,nzeta,ns),   INTENT(IN)   :: xuv
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor), INTENT(OUT)  :: xmn
        INTEGER, INTENT(IN):: ipar, js
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, in
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work1, work2,          &
          work3, work4
        REAL(dp) :: ton, toff
!-----------------------------------------------
        CALL second0(ton)

        ALLOCATE(work1(0:mpol,nzeta), work2(0:mpol,nzeta))
        work1 = zero; work2 = zero

        DO m = 0, mpol
          DO jk = 1, ntheta                                               ! First, add over poloidal collocation angles
            work1(m,:) = work1(m,:) + xuv(jk,:,js)*cosmui(jk, m)          ! Note use of cosmuI/sinmuI; include norm. factors
            work2(m,:) = work2(m,:) + xuv(jk,:,js)*sinmui(jk, m)
          ENDDO
        ENDDO

        xmn = zero                                                        ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta), work4(0:mpol,nzeta))
        work3 =zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                                ! Then, add over toroidal collocation angles
            IF (ipar .EQ. f_cos) THEN                                       ! COSINE series
              work3(:,lk) = work1(:,lk)*cosnv(lk, n)
              work4(:,lk) = work2(:,lk)*sinnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) - work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                         &
                 + work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ELSE                                                          ! SINE series
              work3(:,lk) = work1(:,lk)*sinnv(lk, n)
              work4(:,lk) = work2(:,lk)*cosnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) + work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                         &
                 - work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(0,-ntor:-1) = zero                                            ! Redundant: To avoid counting (0,-n) for non-zero n.
        xmn(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:)
         
        DEALLOCATE(work1, work2, work3, work4)

        CALL second0(toff)
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP2

        
        SUBROUTINE TOIJSP2(XMN, XUV, IPAR, FCOMB, JS)
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), DIMENSION(:,:), INTENT(IN)     :: xmn
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  :: xuv
        INTEGER, INTENT(IN):: ipar, js, fcomb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(dp), ALLOCATABLE  :: xmns(:,:,:)
!-----------------------------------------------
        ALLOCATE(xmns(0:mpol,-ntor:ntor,1))
        xmns(:,:,1) = xmn
        
        CALL toijsp(xmns, xuv(:,:,js:js), fcomb, ipar)
        
        DEALLOCATE(xmns)
        
        END SUBROUTINE TOIJSP2

!-------------------------------------------------------------------------------
!>  @brief Initalize fourier routines.
!>
!>  This subroutine computes the cosine-sine factors that will be needed when
!>  moving between Fourier and real space. All normalizations are contained in
!>  the poloidal quantities used in the Fourier to real space transformation:
!>  * sinmui
!>  * cosmui
!>
!>  Fourier representations are assumed to have stellarator symmetry:
!>  *# cosine (iparity=0):      C(u, v) = sum_u sum_v (C_mn COS(mu + n*nfp*v))
!>  *# sine   (iparity=1):      S(u, v) = sum_u sum_v (S_mn SIN(mu + n*nfp*v))
!>
!>  The number of collocation points have been set initalially equal to the
!>  number of modes:
!>    theta_j = j*pi/M,         j = 0,...,M   Where M = mpol - 1
!>    zeta_k  = k*2*pi/(2N + 1) k = 0,...,2N  Where N = ntor.
!-------------------------------------------------------------------------------
      SUBROUTINE init_fourier
      USE shared_data, ONLY: lasym

      IMPLICIT NONE

!  Local variables
      INTEGER             :: jk
      INTEGER             :: m
      INTEGER             :: n
      INTEGER             :: lk
      INTEGER             :: istat
      INTEGER             :: ntheta1
      INTEGER             :: mpol1
      INTEGER             :: ntor1
      REAL (dp)           :: arg
      REAL (dp)           :: cosi
      REAL (dp)           :: sini
      REAL (dp)           :: pinorm

!  Local parameters
      REAL(dp), PARAMETER :: nfactor = 1.41421356237310_dp ! RS: Normalization of m=0/n=0 Fourier modes.

!  Start of executable code

!  When using stellarator symmetry, we only need half of the poloidal grid.
      IF (lasym) THEN
         ntheta1 = ntheta
         pinorm = twopi
      ELSE
         ntheta1 = ntheta - 1
         pinorm = pi
      END IF

      mpol1 = MAX(mpol_v, mpol)
      ntor1 = MAX(ntor_v, ntor)
       
!       EXTEND mpol, ntor FOR POST-PROCESSING OF PRESSURE B dot GRAD P

      ALLOCATE(cosmu(ntheta,0:mpol1), cosmum(ntheta,0:mpol1),                  &
               cosmui(ntheta,0:mpol1), sinmu(ntheta,0:mpol1),                  &
               sinmum(ntheta,0:mpol1), sinmui(ntheta,0:mpol1),                 &
               cosnv(nzeta,0:ntor1), cosnvn(nzeta,0:ntor1),                    &
               sinnv(nzeta,0:ntor1), sinnvn(nzeta,0:ntor1), stat=istat)
      CALL assert_eq(0, istat, 'Allocation error in init_fourier')
        
      dnorm_i = one/(ntheta1*nzeta) ! ST:Norm for Fourier-Space -> UV-Space

      DO jk = 1, ntheta
         arg = pinorm*(jk-1)/ntheta1  ! Poloidal collocation points*m
         DO m = 0, mpol1
            cosi = COS(m*arg)
            sini = SIN(m*arg)

            CALL Round_Trig(cosi, sini)

            cosmu (jk, m) = cosi         ! Used in UV-SPace to Fourier-space
            cosmum(jk, m) = m*cosi       ! Used in UV-SPace to Fourier-space
            cosmui(jk, m) = dnorm_i*cosi ! Used in Fourier-Space to UV-space
            sinmu (jk, m) = sini         ! Used in UV-SPace to Fourier-space
            sinmum(jk, m) = m*sini       ! Used in UV-SPace to Fourier-space
            sinmui(jk, m) = dnorm_i*sini ! Used in Fourier-Space to UV-space
            IF (.NOT.lasym .AND. (jk.EQ.1 .OR. jk.EQ.ntheta)) THEN
               cosmui(jk, m) = 0.5_dp*cosmui(jk, m) ! 0, pi poloidal points are divided by 2 when integrating
               sinmui(jk, m) = 0.5_dp*sinmui(jk, m)
            END IF
         END DO
      END DO

!  For now we are in one field period.
      DO lk = 1, nzeta
         arg = twopi*(lk-1)/nzeta                                ! Toroidal collocation points*ABS(n*nfp_i)
         DO n = 0, ntor1
            cosi = COS(n*arg)                                     ! sign of "n" taken care later through ISIGN
            sini = SIN(n*arg)
            
            CALL Round_Trig(cosi, sini)

            cosnv (lk, n) = cosi                                 ! All normalization factors goes in poloidal quantities
            cosnvn(lk, n) = n*nfp*cosi                           ! NFP is number of field periods
            sinnv (lk, n) = sini
            sinnvn(lk, n) = n*nfp*sini
         END DO
      END DO

!  Compute orthonorm factor for cos(mu + nv), sin(mu + nv) basis
!  (Noe cos(mu)cos(nv) basis!) so that
!
!     <cos(mu + nv)*cos(m'u + n'v)>*orthonorm(m,n) = 1 (for m=m',n=n')
!
!  Orthonormal factor is Sqrt(2) for all modes except (0,0) mode. The (0,0) mode
!  is 1.
      ALLOCATE (orthonorm(0:mpol1,-ntor1:ntor1), stat=istat)
      CALL assert_eq(0, istat,'Allocate ORTHONORM failed in FIXARRAY!')
      orthonorm = nfactor ! SQRT(2)
      orthonorm(0,0) = 1
      IF (ntheta .EQ. mpol + 1) THEN
         orthonorm(mpol,:) = 1
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Round trig values to whole number for near numbers.
!>
!>  @param[inout] cosi Cosine value to round.
!>  @param[inout] sini Sine value to round.
!-------------------------------------------------------------------------------
      SUBROUTINE Round_Trig(cosi, sini)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(inout) :: cosi
      REAL(dp), INTENT(inout) :: sini

!  Local variables
      REAL(dp)                :: eps

!  Start of executable code
      eps = 10*EPSILON(eps)

      IF (ABS(cosi - 1)      .LE. eps) THEN
         cosi = 1
         sini = 0
      ELSE IF (ABS(cosi + 1) .LE. eps) THEN
         cosi = -1
         sini = 0
      ELSE IF (ABS(sini - 1) .LE. eps) THEN
         cosi = 0
         sini = 1
      ELSE IF (ABS(sini + 1) .LE. eps) THEN
         cosi = 0
         sini = -1
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocate Fourier routines.
!>
!>  Fouier arrays were allocated in @ref init_fourier.
!-------------------------------------------------------------------------------
      SUBROUTINE dealloc_fourier

      IMPLICIT NONE

!  Local variables
      INTEGER :: istat

      DEALLOCATE(cosmu, cosmum, cosmui, sinmu, sinmum, sinmui,                 &
                 cosnv, cosnvn, sinnv, sinnvn, orthonorm, stat=istat) 
      CALL assert_eq(0, istat,'Deallocation error in dealloc_fourier')
         
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Test fourier routines.
!>
!>  This tests the routines by fourier transforming form flux to real space then
!>  back again.
!-------------------------------------------------------------------------------
      SUBROUTINE test_fourier
      USE descriptor_mod, ONLY: iam, nprocs
      USE shared_data, ONLY: lasym, lverbose
      USE timer_mod, ONLY: test_fourier_time

      IMPLICIT NONE

!  Local variables
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: lk
      INTEGER                                 :: istat
      INTEGER                                 :: n1
      INTEGER                                 :: iblock
      INTEGER                                 :: ns_save
      REAL (dp)                               :: arg
      REAL (dp)                               :: cosi
      REAL (dp)                               :: sini
      REAL (dp)                               :: testcc
      REAL (dp)                               :: testsc
      REAL (dp)                               :: testcs
      REAL (dp)                               :: testss
      REAL (dp)                               :: ton
      REAL (dp)                               :: toff
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: testr
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: testi
!-----------------------------------------------
      CALL second0(ton)

      ns_save = ns
      ns = 1
!CHECK ORTHONORMALITY OF BASIS FUNCTIONS
      ALLOCATE (testr(ntheta,nzeta,ns), testi(0:mpol,-ntor:ntor,ns),           &
                stat=istat)
      testr = 0
      testi = 0
      iblock= 0

      NTOR_LOOP: DO n = -ntor, ntor
         MPOL_LOOP: DO m = 0, mpol
            iblock = iblock + 1
            IF (MOD(iblock - 1, nprocs) .NE. iam) CYCLE
            IF (m.EQ.0 .AND. n.LT.0) CYCLE

!TEST cos(mu+nv)
            testi = 0;  testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_none, f_cos)
            CALL tomnsp(testr, testi, f_cos)
            testi(m,n,ns) = testi(m,n,ns) - 1
            IF (ANY(testi(:,:,ns) .GT. four_tol)) GOTO 100
            IF (lasym) THEN
               CALL tomnsp(testr, testi, f_sin)
               IF (ANY(testi(:,:,ns) .GT. four_tol)) GOTO 100
            END IF

!TEST sin(mu+nv)
            IF (m.EQ.0 .AND. n.EQ.0) CYCLE
            testi = 0;  testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_none, f_sin)
            CALL tomnsp(testr, testi, f_sin)
            testi(m,n,ns) = testi(m,n,ns) - 1
            IF (ANY(testi(:,:,ns) .GT. four_tol)) GOTO 100
            IF (lasym) THEN
               CALL tomnsp(testr, testi, f_cos)
               IF (ANY(testi(:,:,ns) .GT. four_tol)) GOTO 100
            END IF

!TEST U and V DERIVATIVES OF cos(mu+nv)
            testi = 0
            testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_du, f_cos)
            CALL tomnsp(testr, testi, f_sin)
            testi(m,n,ns) = testi(m,n,ns) + m
            IF (ANY(ABS(testi(:,:,ns)) .GT. four_tol)) GOTO 100
            testi = 0
            testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_dv, f_cos)
            CALL tomnsp(testr, testi, f_sin)
            testi(m,n,ns) = testi(m,n,ns) + n*nfp
            IF (ANY(ABS(testi(:,:,ns)) .GT. four_tol)) GOTO 100

!TEST U and V DERIVATIVES (sin(mu+nv))
            testi = 0
            testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_du, f_sin)
            CALL tomnsp(testr, testi, f_cos)
            testi(m,n,ns) = testi(m,n,ns) - m
            IF (ANY(ABS(testi(:,:,ns)) .GT. four_tol)) GOTO 100
            testi = 0
            testi(m,n,ns) = 1
            CALL toijsp(testi, testr, f_dv, f_sin)
            CALL tomnsp(testr, testi, f_cos)
            testi(m,n,ns) = testi(m,n,ns) - n*nfp
            IF (ANY(ABS(testi(:,:,ns)) .GT. four_tol)) GOTO 100

         END DO MPOL_LOOP
      END DO NTOR_LOOP

      GOTO 200

 100  CONTINUE

      WRITE (*,1000) m, n
      CALL ASSERT(.false., 'FOURIER ORTHOGONALITY TEST FAILED!')

 200  CONTINUE

      ns = ns_save
      DEALLOCATE (testr, testi)

      CALL second0(toff)
      test_fourier_time=toff-ton

      IF (iam .EQ. 0) THEN
         IF (.NOT.lverbose) THEN
            WRITE(*,1001) toff - ton
         END IF
         WRITE(unit_out,1001) toff - ton
      END IF

1000  FORMAT('ORTHOGONALITY TEST FAILED FOR MODES ',i3,1x,i3)
1001  FORMAT(1x,'ORTHOGONALITY TEST SUCCEEDED - TIME: ',1p,e12.2,' s')

      END SUBROUTINE TEST_FOURIER


      END MODULE fourier
