      MODULE fourier
!     
!     WRITTEN 08-09-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: CONVERTS QUANTITIES FROM FOURIER SPACE TO REAL SPACE AND VICEVERSA
!
!       NOTE: CALL FIXARRAY must be used once before calling any of the Fourier subroutines
!       to calculate all necessary cosine/sine factors on fixed mesh of collocation angles
!       NS, NTOR, MPOL and NFP used as in the wout.file
!
      USE stel_kinds
      USE stel_constants
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,            &
                   mpol=>mpol_i, ntor=>ntor_i, nfp=>nfp_i
      USE timer_mod
      USE descriptor_mod, ONLY: INHESSIAN, DIAGONALDONE
      IMPLICIT NONE
!
!     PARAMETER DECLARATIONS
!     Partity flags.
      INTEGER, PARAMETER :: f_cos = 0
      INTEGER, PARAMETER :: f_sin = 1

!     Derivative and sum flags.
      INTEGER, PARAMETER :: f_none = 0   ! 000
      INTEGER, PARAMETER :: f_ds   = 1   ! 001
      INTEGER, PARAMETER :: f_du   = 2   ! 010
      INTEGER, PARAMETER :: f_dv   = 4   ! 100
      INTEGER, PARAMETER :: f_dudv = 6   ! 110
      INTEGER, PARAMETER :: f_sum  = 8   !1000

!     Bit positions for Derivative flags.
      INTEGER, PARAMETER :: b_ds   = 0
      INTEGER, PARAMETER :: b_du   = 1
      INTEGER, PARAMETER :: b_dv   = 2
      INTEGER, PARAMETER :: b_sum  = 3

!
!     VARIABLE DECLARATIONS  
!
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: orthonorm
      CONTAINS
 
!#if defined(SKS)
      SUBROUTINE TOIJSP_PAR(XMN, XUV, DFLAG, IPARITY, nsmin, nsmax)
!
!       DESCRIPTION:  This subroutine moves a quantity X to real space by summing over its Fourier harmonics X_MN.
!                     IPARITY = 0: DFLAG=0, COS(mu+nv) ; DFLAG&1b, -mSIN(mu+nv); DFLAG&10b, -nSIN(mu+nv)
!                             = 1: DFLAG=0, SIN(mu+nv) ; DFLAG&1b,  mCOS(mu+nv); DFLAG&10b,  nCOS(mu+nv)
!       IF DFLAG = 0b01, the first poloidal derivative of that quantity is obtained
!       IF DFLAG = 0b10, the first toroidal derivative of that quantity is obtained
!       IF DFLAG = 0b11, it gives the joint poloidal, toroidal derivative
!
      USE stel_kinds
      USE mpi_inc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      INTEGER, INTENT(IN):: ihalf                                          ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
      INTEGER, INTENT(IN):: nsmin, nsmax 
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN) :: xmn
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(INOUT)   :: xuv
      INTEGER, INTENT(IN):: dflag                                          ! DFLAG = order of poloidal/toroidal derivative
      INTEGER, INTENT(IN):: iparity                                        ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, jk, lk, js
      INTEGER :: in 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work0, work1, work2, xuv1
      REAL(dp) :: skston, skstoff
      INTEGER :: counter, nx
      INTEGER, SAVE :: PTOIJSPPASS=1
!-----------------------------------------------
      CALL second0(skston)

      ALLOCATE(work0(0:mpol,-ntor:ntor), xuv1(ntheta,nzeta),             &         
               work1(ntheta,-ntor:ntor), work2(ntheta,-ntor:ntor),       &
               stat=lk)
      IF (lk .NE. 0) STOP 'Allocation error in TOIJSP_PAR'
 
      NS_LOOP: DO js = nsmin, nsmax

      work0(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
      work1 = zero; work2 = zero 

      DO jk = 1, ntheta
        DO m = 0, mpol                                                     ! First DO sum over poloidal modes
          IF (BTEST(dflag, b_du)) THEN                                     ! Poloidal derivative requested: USE COSMUM, SINMUM
            work1(jk,:) = work1(jk,:) - work0(m,:)*sinmum(jk, m)
            work2(jk,:) = work2(jk,:) + work0(m,:)*cosmum(jk, m)
          ELSE
            work1(jk,:) = work1(jk,:) + work0(m,:)*cosmu(jk, m)
            work2(jk,:) = work2(jk,:) + work0(m,:)*sinmu(jk, m)
          END IF
        ENDDO
      ENDDO
 
      xuv1 = zero                                                          ! Make sure that SPATIAL variable is zeroed

      DO lk = 1, nzeta

        IF (iparity == f_cos) THEN                                         ! Do first N=0 mode

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
          ENDIF

        ENDIF          

        DO n = 1, ntor                                                     ! Then sum over N>0 and N<0 toroidal modes                                                     
          IF (iparity == f_cos) THEN                                       ! COSINE series

            IF (BTEST(dflag, b_dv)) THEN                                   ! First toroidal derivative requested
              xuv1(:,lk) =  xuv1(:,lk)                                   &
              - (work1(:,n) + work1(:,-n))*sinnvn(lk, n)                 &
              - (work2(:,n) - work2(:,-n))*cosnvn(lk, n)
            ELSE
              xuv1(:,lk) =  xuv1(:,lk)                                   &
              + (work1(:,n) + work1(:,-n))*cosnv(lk, n)                  &
              - (work2(:,n) - work2(:,-n))*sinnv(lk, n)
            ENDIF

          ELSE                                                             ! SINE series

            IF (BTEST(dflag, b_dv)) THEN                                   ! First toroidal derivative requested
              xuv1(:,lk) =  xuv1(:,lk)                                   &
              + (work1(:,n) - work1(:,-n))*cosnvn(lk, n)                 &
              - (work2(:,n) + work2(:,-n))*sinnvn(lk, n)
            ELSE
              xuv1(:,lk) =  xuv1(:,lk)                                   &
              + (work1(:,n) - work1(:,-n))*sinnv(lk, n)                  &
              + (work2(:,n) + work2(:,-n))*cosnv(lk, n)
            ENDIF

          ENDIF
        ENDDO
      ENDDO

      IF (BTEST(dflag, b_sum)) THEN
         xuv(:,:,js) = xuv(:,:,js) + xuv1(:,:)
      ELSE
         xuv(:,:,js) = xuv1(:,:)
      END IF

      END DO NS_LOOP

      DEALLOCATE(work0, work1, work2, xuv1)

      CALL second0(skstoff)
      IF(DIAGONALDONE.AND.INHESSIAN) toijsp_time=toijsp_time+(skstoff-skston)
      time_toijsp = time_toijsp + (skstoff-skston)

      PTOIJSPPASS = PTOIJSPPASS+1

      END SUBROUTINE TOIJSP_PAR

      SUBROUTINE TOMNSP_PAR(XUV, XMN, IPARITY, nsmin, nsmax)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!
      USE stel_kinds
      USE mpi_inc
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      INTEGER, INTENT(IN) :: ihalf                                        ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(IN):: xuv
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(OUT):: xmn
      INTEGER, INTENT(IN):: iparity                                       ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, jk, lk, js
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work1, work2, xuv1, xmn1
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: t1, t2, work3, work4
      REAL(dp) :: skston, skstoff
      INTEGER :: ii, jj, kk 
      INTEGER :: counter, nx
      INTEGER, SAVE :: PTOMNSPPASS=1
!-----------------------------------------------

      CALL second0(skston)

      ALLOCATE(work1(0:mpol,nzeta), work2(0:mpol,nzeta),                 &
               work3(0:mpol), work4(0:mpol),                             &
               xuv1(ntheta,nzeta),  xmn1(0:mpol,-ntor:ntor),             &
               t1(nzeta), t2(nzeta), stat=ii)
      IF (ii .NE. 0) STOP 'Allocation error in tomnsp_par'

      NS_LOOP: DO js = nsmin, nsmax

      xuv1(:,:) = xuv(:,:,js)
      xmn1      = 0

      DO m = 0, mpol
        t1 = 0; t2 = 0
        DO jk = 1, ntheta                                                 ! First, add over poloidal collocation angles
          t1 = t1 + xuv1(jk,:)*cosmui(jk, m)                              ! Note use of cosmuI/sinmuI; include norm. factors
          t2 = t2 + xuv1(jk,:)*sinmui(jk, m)
        ENDDO
        work1(m,:) = t1;  work2(m,:) = t2
      ENDDO


      DO n = 0, ntor  
        DO lk = 1, nzeta                                                  ! Then, add over toroidal collocation angles
          IF (iparity == f_cos) THEN                                      ! COSINE series
            work3(:) = work1(:,lk)*cosnv(lk, n)
            work4(:) = work2(:,lk)*sinnv(lk, n)
            xmn1(:,n) =  xmn1(:,n) + work3(:) - work4(:)
            IF (n .NE. 0) THEN
              xmn1(1:,-n) = xmn1(1:,-n) + work3(1:) + work4(1:)
            ENDIF
          ELSE                                                            ! SINE series
            work3(:) = work1(:,lk)*sinnv(lk, n)
            work4(:) = work2(:,lk)*cosnv(lk, n)
            xmn1(:,n) =  xmn1(:,n) + work3(:) + work4(:)
            IF (n .NE. 0) THEN
              xmn1(1:,-n) = xmn1(1:,-n) - work3(1:) + work4(1:)
            ENDIF
          ENDIF
        ENDDO
      ENDDO 

!NOT NECESSARY      xmn1(0,-ntor:-1) = zero                               ! Redundant: To avoid counting (0,-n) for non-zero n.

      xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn1(:,:)

      END DO NS_LOOP

      DEALLOCATE(work1, work2, work3, work4, xuv1, xmn1, t1, t2)

      CALL second0(skstoff)
      IF(DIAGONALDONE.AND.INHESSIAN) tomnsp_time=tomnsp_time+(skstoff-skston)

      time_tomnsp = time_tomnsp + (skstoff-skston)

      PTOMNSPPASS = PTOMNSPPASS+1
      END SUBROUTINE TOMNSP_PAR
!#endif

        SUBROUTINE GETORIGIN (XUV, MODE, IPARITY)
!
!       DESCRIPTION: Computes the origin value of a half-mesh quantity (stores it at index 1)
!                    from the m=mode (cos,sin(mu)-average) of the quantity at js=2
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       MODE         the poloidal mode number to be retained
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(in) :: iparity, mode
        REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(INOUT) :: xuv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: js1=1, js2=2
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor) :: xmn
!-----------------------------------------------
        CALL tomnsp2(xuv, xmn, iparity, js2)
        
        IF (mode .GT. 0) xmn(0:mode-1,:) = 0
        IF (mode .LT. mpol) xmn(mode+1:,:) = 0

        CALL toijsp2(xmn, xuv, iparity, js1)

        END SUBROUTINE GETORIGIN

        SUBROUTINE TOIJSP(XMN, XUV, DFLAG, IPARITY)
!
!       DESCRIPTION:  This subroutine moves a quantity X to real space by summing over its Fourier harmonics X_MN.
!       IF DFLAG = 0b01, the first poloidal derivative of that quantity is obtained
!       IF DFLAG = 0b10, the first toroidal derivative of that quantity is obtained
!       IF DFLAG = 0b11, it gives the joint poloidal, toroidal derivative
!
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!        INTEGER, INTENT(IN):: ihalf                                        ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns),                       &
          INTENT(IN) :: xmn
        REAL(dp), DIMENSION(ntheta,nzeta,ns),                            &
          INTENT(INOUT) :: xuv
        INTEGER, INTENT(IN) :: dflag                                       ! DFLAG = order of poloidal/toroidal derivative
        INTEGER, INTENT(IN) :: iparity                                     ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)

        CALL TOIJSP_PAR(XMN, XUV, DFLAG, IPARITY, 1, ns)

        END SUBROUTINE TOIJSP


        SUBROUTINE TOMNSP(XUV, XMN, IPARITY)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even in mu+nv; = 1, X is odd
!
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!        INTEGER, INTENT(IN):: ihalf                                        ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
        REAL(dp), DIMENSION(ntheta,nzeta,ns),                            &
          INTENT(IN):: xuv
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns),                       &
          INTENT(OUT):: xmn
        INTEGER, INTENT(IN):: iparity                                      ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------

        CALL TOMNSP_PAR(XUV, XMN, IPARITY, 1, ns)

        END SUBROUTINE TOMNSP
      
!SPH 10/15/14
        SUBROUTINE TOMNSP_PEST(XUV, XMN, LMNS, IPARITY, LMNC)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics X_MN by
!         summing over all the presribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       XUV is in VMEC coordinates; the XMN are harmonics in PEST coordinates, where
!              UP(s,u,v) = U + LAMBDA(s,u,v)
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even in mu+nv; = 1, X is odd
!
        USE stel_kinds
        USE shared_data, ONLY: lasym

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!        INTEGER, INTENT(IN):: ihalf                                        ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH
        REAL(dp), DIMENSION(ntheta,nzeta,ns),                            &
          INTENT(IN):: xuv
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns),                       &
          INTENT(OUT):: xmn
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns),                       &
          INTENT(IN) :: lmns
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns),                       &
          INTENT(IN), OPTIONAL :: lmnc
        INTEGER, INTENT(IN):: iparity                                      ! IPARITY = 0, cosine (EVEN); = 1; sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, js, istat
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: work1, work2,         &
          work3, work4, lambda, dupdu, cos1l, sin1l, cosml, sinml
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: cosmuip, sinmuip
        REAL(dp) :: ton, toff, skston, skstoff
!-----------------------------------------------
        CALL second0(ton)
        skston=ton

!
!       1. COMPUTE LAMBDA
!
        ALLOCATE (lambda(ntheta,nzeta,ns), dupdu(ntheta,nzeta,ns),       &  
                  cos1l(ntheta,nzeta,ns), cosml(ntheta,nzeta,ns),        &
                  sin1l(ntheta,nzeta,ns), sinml(ntheta,nzeta,ns),        &
                  stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #1 in TOMNSP_PEST'
        
        CALL TOIJSP(lmns, lambda, f_none, f_sin)
        CALL TOIJSP(lmns, dupdu,  f_du, f_sin)
        IF (PRESENT(lmnc)) THEN
           CALL TOIJSP(lmns, lambda, f_sum, f_cos)
           CALL TOIJSP(lmns, dupdu,  IOR(f_sum,f_du), f_cos)
        END IF

!       transform from up to u
        dupdu = 1 + dupdu                                              

!       INTERPOLATE TO FULL GRID IF ihalf=0
!        IF (ihalf .EQ. 0) THEN
           lambda(:,:,1) = lambda(:,:,2)                                  !could do better, why bother?
           dupdu (:,:,1) = dupdu (:,:,2)
           DO js=2,ns-1
              lambda(:,:,js) = .5_dp*(lambda(:,:,js)+lambda(:,:,js+1))
              dupdu (:,:,js) = .5_dp*(dupdu (:,:,js)+dupdu (:,:,js+1))
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
        ALLOCATE(work1(0:mpol,nzeta,ns), work2(0:mpol,nzeta,ns),         &
                 work3(ntheta,nzeta,ns), cosmuiP(nzeta,ns),              &
                 sinmuiP(nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #2 in TOMNSP_PEST'
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
              cosmuiP(:,1:ns) = cosmui(jk,m)*cosml(jk,:,1:ns)            &
                      - sinmui(jk,m)*sinml(jk,:,1:ns)
              sinmuiP(:,1:ns) = sinmui(jk,m)*cosml(jk,:,1:ns)            &
                      + cosmui(jk,m)*sinml(jk,:,1:ns)
              work1(m,:,1:ns) = work1(m,:,1:ns)                          &
                               + xuv(jk,:,1:ns)*cosmuiP(:,1:ns)          ! Note use of cosmuI/sinmuI; include norm. factors
              work2(m,:,1:ns) = work2(m,:,1:ns)                          &
                               + xuv(jk,:,1:ns)*sinmuiP(:,1:ns)
           ENDDO
        END DO

        DEALLOCATE (work3, cosmuip, sinmuip)
        xmn = zero                                                       ! Make sure that FOURIER variable is zeroed
        ALLOCATE(work3(0:mpol,nzeta,ns), work4(0:mpol,nzeta, ns),        &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation error #3 in TOMNSP_PEST'
        work3 = zero; work4 = zero

        DO n = 0, ntor  
          DO lk = 1, nzeta                                               ! Then, add over toroidal collocation angles
            IF (iparity == f_cos) THEN                                   ! COSINE series
              work3(:,lk,1:ns) = work1(:,lk,1:ns)*cosnv(lk, n)
              work4(:,lk,1:ns) = work2(:,lk,1:ns)*sinnv(lk, n)
              xmn(:,n,1:ns) =  xmn(:,n,1:ns)                             &
                 + work3(:,lk,1:ns) - work4(:,lk,1:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,1:ns) = xmn(1:,-n,1:ns)                        &
                 + work3(1:,lk,1:ns) + work4(1:,lk,1:ns)
              ENDIF
            ELSE                                                         ! SINE series
              work3(:,lk,1:ns) = work1(:,lk,1:ns)*sinnv(lk, n)
              work4(:,lk,1:ns) = work2(:,lk,1:ns)*cosnv(lk, n)
              xmn(:,n,1:ns) =  xmn(:,n,1:ns)                             &
                + work3(:,lk,1:ns) + work4(:,lk,1:ns)
              IF (n .ne. 0) THEN
                xmn(1:,-n,1:ns) = xmn(1:,-n,1:ns)                        &
                - work3(1:,lk,1:ns) + work4(1:,lk,1:ns)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 
        xmn(0,-ntor:-1,:) = zero                                         ! Redundant: To avoid counting (0,-n) for non-zero n.

        DO js = 1,ns
           xmn(:,:,js) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:,js)
        END DO
         
        DEALLOCATE(work1, work2, work3, work4, lambda, cos1l, sin1l,     &
                   sinml, cosml, dupdu)

        CALL second0(toff)
#if defined(SKS)
        skstoff=toff
        IF(DIAGONALDONE.AND.INHESSIAN) tomnsp_time=tomnsp_time+(skstoff-skston)
#endif
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP_PEST


        SUBROUTINE TOMNSP2(XUV, XMN, IPARITY, JS)
!
!       Description: This subroutine moves a quantity X to Fourier space producing its harmonics XMN by
!         summing over all the prescribed (toroidal and poloidal) collocation points:
!           theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!         where M = mpol - 1 and N = ntor.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives) is even; = 1, X is odd
!       JS:  Surface where xmn is computed
!
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), DIMENSION(ntheta,nzeta,ns),                            &
          INTENT(IN):: xuv
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor),                          &
          INTENT(OUT):: xmn
        INTEGER, INTENT(IN):: iparity, js
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk, in
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work1, work2,           &
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
          DO lk = 1, nzeta                                      ! Then, add over toroidal collocation angles
            IF (iparity == f_cos) THEN                          ! COSINE series
              work3(:,lk) = work1(:,lk)*cosnv(lk, n)
              work4(:,lk) = work2(:,lk)*sinnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) - work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                          &
                 + work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ELSE                                                ! SINE series
              work3(:,lk) = work1(:,lk)*sinnv(lk, n)
              work4(:,lk) = work2(:,lk)*cosnv(lk, n)
              xmn(:,n) =  xmn(:,n) + work3(:,lk) + work4(:,lk)
              IF (n .ne. 0) THEN
                xmn(1:mpol,-n) = xmn(1:mpol,-n)                          &
                 - work3(1:mpol,lk) + work4(1:mpol,lk)
              ENDIF
            ENDIF
          ENDDO
        ENDDO 
 

        xmn(0,-ntor:-1) = zero              ! Redundant: To avoid counting (0,-n) for non-zero n.
        xmn(:,:) = orthonorm(0:mpol,-ntor:ntor)*xmn(:,:)
         
        DEALLOCATE(work1, work2, work3, work4)

        CALL second0(toff)
        time_tomnsp = time_tomnsp + (toff-ton)

        END SUBROUTINE TOMNSP2

        
        SUBROUTINE TOIJSP2(XMN, XUV, IPARITY, JS)
        USE stel_kinds
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), DIMENSION(0:mpol,-ntor:ntor),                          &
          INTENT(IN):: xmn
        REAL(dp), DIMENSION(ntheta,nzeta,ns),                            &
          INTENT(INOUT):: xuv
        INTEGER, INTENT(IN):: iparity, js                                 ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n, jk, lk
        INTEGER :: in 
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work0,                  &
          work1, work2
        REAL(dp) :: ton, toff
!-----------------------------------------------
        CALL second0(ton)

        ALLOCATE(work0(0:mpol,-ntor:ntor))                                 ! RS (10/03/06): Added WORK0 for n=0/m=0 normalization

        work0 = xmn
        work0(:,:) = orthonorm(0:mpol,-ntor:ntor)*work0(:,:)
             
        ALLOCATE(work1(ntheta,-ntor:ntor), work2(ntheta,-ntor:ntor))
        work1 = zero; work2 = zero 

        DO jk = 1, ntheta
          DO m = 0, mpol                                                   ! First DO sum over poloidal modes
              work1(jk,:) = work1(jk,:) + work0(m,:)*cosmu(jk, m)
              work2(jk,:) = work2(jk,:) + work0(m,:)*sinmu(jk, m)
          ENDDO
        ENDDO
        DEALLOCATE(work0)

        xuv(:,:,js) = 0                                                    ! Make sure that SPATIAL variable is zeroed

        DO lk = 1, nzeta
          
          IF (iparity == 0) THEN                                           ! Do first N=0 mode
            xuv(:,lk,js) = work1(:,0) 
          ELSE
            xuv(:,lk,js) = work2(:,0) 
          ENDIF          
          
          DO n = 1, ntor                                                   ! Then sum over N>0 and N<0 toroidal modes                                                     
            IF (iparity == f_cos) THEN                                     ! COSINE series
                xuv(:,lk,js) =  xuv(:,lk,js)                             &
                  + (work1(:,n) + work1(:,-n))*cosnv(lk, n)              &
                  - (work2(:,n) - work2(:,-n))*sinnv(lk, n)
            ELSE                                                           ! SINE series
                xuv(:,lk,js) =  xuv(:,lk,js)                             &
                  + (work1(:,n) - work1(:,-n))*sinnv(lk, n)              &
                  + (work2(:,n) + work2(:,-n))*cosnv(lk, n)
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE(work1, work2)

        CALL second0(toff)
        time_toijsp = time_toijsp + (toff-ton)

        END SUBROUTINE TOIJSP2

      SUBROUTINE INIT_FOURIER
!
!     This subroutine computes the cosine-sine factors that will be needed when moving between
!     Fourier and real space. All normalizations are contained in the poloidal quantities 
!     used in the Fourier to real space transformation: SINMUI, COSMUI
!
!     Fourier representations are assumed to have STELLARATOR symmetry:
!         1. COSINE (iparity=0):      C(u, v) = sum_u sum_v (C_mn COS(mu + n*nfp*v))
!         2. SINE   (iparity=1):      S(u, v) = sum_u sum_v (S_mn SIN(mu + n*nfp*v))
!
!     The number of collocation points have been set initially equal to the number of modes:
!       theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!       where M = mpol - 1 and N = ntor.
!
        USE shared_data, ONLY: lasym
        USE descriptor_mod, ONLY: iam, nprocs
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(dp), PARAMETER :: nfactor = 1.41421356237310_dp          ! RS: Normalization of m=0/n=0 Fourier modes.
        INTEGER        :: jk, m, n, lk, istat, ntheta1
        REAL(dp)    :: arg, cosi, sini
        REAL(dp)    :: pinorm
!-----------------------------------------------
!       TO CHANGE TO FULL u-GRID
!       1) nu_i = 2*nu_i
!       2) ntheta1 = ntheta
!       3) arg = twopi*(jk-1)/ntheta1

        IF (lasym) THEN
           ntheta1 = ntheta
           pinorm = twopi
        ELSE
           ntheta1 = ntheta-1
           pinorm = pi
        END IF

        ALLOCATE(cosmu(ntheta,0:mpol),  cosmum(ntheta,0:mpol),           &
                 cosmui(ntheta,0:mpol), sinmu(ntheta,0:mpol),            &
                 sinmum(ntheta,0:mpol), sinmui(ntheta,0:mpol),           &
                 cosnv(nzeta,0:ntor), cosnvn(nzeta,0:ntor),              &
                 sinnv(nzeta,0:ntor), sinnvn(nzeta,0:ntor), stat=istat)

        IF (istat .NE. 0) STOP 'Allocation error in fixarray'
         
        dnorm_i = one/(ntheta1*nzeta)                                      ! ST:Norm for Fourier-Space -> UV-Space

        DO jk = 1, ntheta
          arg = (pinorm*(jk-1))/ntheta1                                    ! Poloidal collocation points
          DO m = 0, mpol
            cosi = COS(m*arg)
            sini = SIN(m*arg)

            CALL Round_Trig(cosi, sini)

            cosmu (jk, m)  = cosi                                          ! Used in UV-SPace to Fourier-space
            cosmum(jk, m)  = m*cosi                                        ! Used in UV-SPace to Fourier-space
            cosmui(jk, m) = dnorm_i*cosi                                   ! Used in Fourier-Space to UV-space
            sinmu (jk, m)  = sini                                          ! Used in UV-SPace to Fourier-space
            sinmum(jk, m)  = m*sini                                        ! Used in UV-SPace to Fourier-space
            sinmui(jk, m) = dnorm_i*sini                                   ! Used in Fourier-Space to UV-space
            IF (.NOT.lasym .AND. (jk.EQ.1 .OR. jk.EQ.ntheta)) THEN
               cosmui(jk, m) = cosmui(jk, m)/2                             ! end poloidal points are divided by 2 when integrating
               sinmui(jk, m) = sinmui(jk, m)/2
            ENDIF
          ENDDO
        ENDDO

!
!FOR NOW, WE ARE IN 1 FIELD PERIOD
!
        DO lk = 1, nzeta
          arg = twopi*(lk-1)/nzeta                                         ! Toroidal collocation points*ABS(n*nfp_i)
          DO n = 0, ntor
            cosi = COS(n*arg)                                              ! sign of "n" taken care later through ISIGN
            sini = SIN(n*arg)
            
            CALL Round_Trig(cosi, sini)

            cosnv (lk, n)  = cosi                                          ! All normalization factors goes in poloidal quantities
            cosnvn(lk, n)  = n*nfp*cosi                                    ! NFP is number of field periods
            sinnv (lk, n)  = sini
            sinnvn(lk, n)  = n*nfp*sini
          ENDDO
        ENDDO

!COMPUTE orthonorm FACTOR FOR cos(mu+nv), sin(mu+nv) BASIS (NOT cos(mu)cos(nv) basis!)
!SO THAT <cos(mu+nv)*cos(m'u+n'v)> * orthonorm(m,n) = 1 (for m=m',n=n')
!BE SURE TO FIX THIS WHEN DOING ASYMMETRIC CASE!

        ALLOCATE (orthonorm(0:mpol,-ntor:ntor), stat=istat)
        IF (istat .NE. 0) STOP 'Allocate ORTHONORM failed in FIXARRAY!'
        orthonorm = nfactor                                       ! SQRT(2)
        orthonorm(0,0) = 1
        IF (ntheta .EQ. (mpol+1)) orthonorm(mpol,:) = 1

        CALL Test_Fourier

      END SUBROUTINE INIT_FOURIER

      SUBROUTINE ROUND_TRIG(cosi, sini)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), INTENT(inout) :: cosi, sini
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        REAL(dp) :: eps
!-----------------------------------------------
        eps = 10*EPSILON(eps)

        IF (ABS(cosi-1) .le. eps) THEN
           cosi = 1
           sini = 0
        ELSE IF (ABS(cosi+1) .le. eps) THEN
           cosi = -1
           sini = 0
        ELSE IF (ABS(sini-1) .le. eps) THEN
           cosi = 0
           sini = 1
        ELSE IF (ABS(sini+1) .le. eps) THEN
           cosi = 0
           sini = -1
        END IF

      END SUBROUTINE Round_Trig


      SUBROUTINE DEALLOC_FOURIER
      INTEGER :: istat

      DEALLOCATE(cosmu, cosmum, cosmui, sinmu, sinmum, sinmui,           &
                 cosnv, cosnvn, sinnv, sinnvn, orthonorm, stat=istat) 
      IF (istat .ne. 0) STOP 'Allocation error in dealloc_fixarray'
         
      END SUBROUTINE DEALLOC_FOURIER

      SUBROUTINE TEST_FOURIER
      USE descriptor_mod, ONLY: iam, nprocs
      USE shared_data, ONLY: lasym, lverbose
#if defined(SKS)
      USE timer_mod, ONLY: test_fourier_time
#endif      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER             :: m, n, lk, istat, n1, iblock, ns_save
	  REAL(dp), PARAMETER :: FOUR_TOL = 1.E-12_dp
	  REAL(dp)            :: ton, toff
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: testr, testi
!-----------------------------------------------
      CALL second0(ton)

      ns_save = ns
      ns = 1
!CHECK ORTHONORMALITY OF BASIS FUNCTIONS
      ALLOCATE (testr(ntheta,nzeta,ns), testi(0:mpol,-ntor:ntor,ns),     &
                stat=istat)
      testr = 0
      testi = 0
      iblock= 0

      NTOR_LOOP: DO n = -ntor, ntor
      MPOL_LOOP: DO m = 0, mpol
          iblock = iblock+1
          IF (MOD(iblock-1, nprocs) .NE. iam) CYCLE
          IF (m.EQ.0 .AND. n.LT.0) CYCLE

!TEST cos(mu+nv)
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_none, f_cos)
          CALL tomnsp(testr, testi, f_cos)
          testi(m,n,ns) = testi(m,n,ns) - 1
          IF (ANY(testi(:,:,ns) .GT. FOUR_TOL)) GOTO 100
          IF (lasym) THEN
          CALL tomnsp(testr, testi, f_sin)
          IF (ANY(testi(:,:,ns) .GT. FOUR_TOL)) GOTO 100
          END IF

          IF (m.EQ.0 .AND. n.EQ.0) CYCLE

!TEST sin(mu+nv)
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_none, f_sin)
          CALL tomnsp(testr, testi, f_sin)
          testi(m,n,ns) = testi(m,n,ns) - 1
          IF (ANY(testi(:,:,ns) .GT. FOUR_TOL)) GOTO 100
          IF (lasym) THEN
          CALL tomnsp(testr, testi, f_cos)
          IF (ANY(testi(:,:,ns) .GT. FOUR_TOL)) GOTO 100
          END IF

!TEST U and V DERIVATIVES OF cos(mu+nv)
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_du, f_cos)
          CALL tomnsp(testr, testi, f_sin)
          testi(m,n,ns) = testi(m,n,ns) + m
          IF (ANY(ABS(testi(:,:,ns)) .GT. FOUR_TOL)) GOTO 100
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_dv, f_cos)
          CALL tomnsp(testr, testi, f_sin)
          testi(m,n,ns) = testi(m,n,ns) + n*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. FOUR_TOL)) GOTO 100

!TEST U and V DERIVATIVES (sin(mu+nv))
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_du, f_sin)
          CALL tomnsp(testr, testi, f_cos)
          testi(m,n,ns) = testi(m,n,ns) - m
          IF (ANY(ABS(testi(:,:,ns)) .GT. FOUR_TOL)) GOTO 100
          testi = 0;  testi(m,n,ns) = 1
          CALL toijsp(testi, testr, f_dv, f_sin)
          CALL tomnsp(testr, testi, f_cos)
          testi(m,n,ns) = testi(m,n,ns) - n*nfp
          IF (ANY(ABS(testi(:,:,ns)) .GT. FOUR_TOL)) GOTO 100

         END DO MPOL_LOOP
       END DO NTOR_LOOP

       GOTO 200

 100   CONTINUE
       PRINT *,'ORTHOGONALITY TEST FAILED FOR MODES ',m,n
       STOP

 200   CONTINUE

       ns = ns_save
       DEALLOCATE (testr, testi)

       CALL second0(toff)
#if defined(SKS)
       test_fourier_time=toff-ton
#endif       
       IF (iam .EQ. 0) THEN
          IF (lverbose)PRINT 300,'ORTHOGONALITY TEST SUCCEEDED - TIME: ', toff-ton,' s'
          WRITE (33, 300) 'ORTHOGONALITY TEST SUCCEEDED - TIME: ', toff-ton,' s'
       END IF

 300   FORMAT(1x,a,1p,e12.2,a)

       END SUBROUTINE TEST_FOURIER


      END MODULE fourier
