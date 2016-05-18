      SUBROUTINE scalpot(bvec, amatrix, wint, ns, ivacskip)
      USE vacmod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ns, ivacskip
      REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      REAL(rprec), INTENT(in) :: wint(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ip, istore, istart, ndim, i
      REAL(rprec), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
      REAL(rprec), ALLOCATABLE :: greenp(:,:)
!-----------------------------------------------
      IF (.NOT.ALLOCATED(amatsav))
     1   STOP 'AMATSAV: Allocation error in scalpot'

      ALLOCATE (grpmn(nuv2*mnpd2), stat=ip)
      IF (ip .NE. 0) STOP 'GRPMN: Allocation error in scalpot'

!
!     COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!     ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!     AND GRPMN CONTAINS TRANSFORM OF NORMAL DERIVATIVE
!     OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
!
!     FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED
!
      ndim = mnpd2/mnpd
      CALL analyt (grpmn, bvec, ivacskip, ndim)

      IF (ivacskip .NE. 0) THEN
         bvec = bvec + bvecsav
      ELSE

         ALLOCATE (gstore(nuv), stat=ip)
         IF (ip .NE. 0) STOP 'Allocation error in SCALPOT'

!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP& PRIVATE(green,greenp)

!$OMP WORKSHARE
         bvecsav = bvec
         gstore  = 0
!$OMP END WORKSHARE

         ALLOCATE(green(nuv), greenp(nuv,istore_max))
!
!        COMPUTE SURFACE INTEGRALS OF SOURCE AND GREENS FUNCTION NEEDED
!        FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
!        NOTE: SOURCE IS THE RHS OF EQ.(3.2), KERNEL IS THE LHS OF EQ (3.2).
!        IP IS THE INDEX OF THE PRIMED VARIABLE MESH.
!

!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(i,ip,istart,istore)
         DO istart = 0, nuv2 - 1, istore_max
            DO istore = 1, istore_max

!        COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT
!        [FIRST TERMS IN EQ.(2.14, 2.16)].
!
!        BECAUSE OF THE LARGE SIZES INVOLVED (NUV2*NUV2), THIS IS DONE HERE BY STORING
!        THESE QUANTITIES - FOR ALL VALUES OF THE UNPRIMED U,V COORDINATE MESH - ON A
!        LIMITED SUBSET (ISTORE_max) OF PRIMED MESH VALUES (INDEXED BY ISTORE~IP), 
!        THE FOURIER TRANSFORM OVER THE PRIMED MESH IS "BUILT-UP" BY MULTIPLE CALLS TO FOURP
!        WITHIN THIS LOOP.

               ip = istore + istart
               CALL greenf(green, greenp(:,istore), ip)

               DO i = 1, nuv
!
!        PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM
!        [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
!        AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE

!  Note there can appear to be randomness when running multiple times. This
!  occurs from the different treats updating the gstore variable in a random
!  order causing differences in the floating point round off. The differences
!  in the results are on the order of machine precision.
!$OMP ATOMIC
                  gstore(i) = gstore(i) + bexni(ip)*green(i)
               END DO

               IF (ip .eq. nuv2 .or. istore .eq. istore_max) THEN
                  EXIT
               END IF
            END DO
!
!        PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
!        (FOR )
!        AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
!        - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)

            CALL fourp(grpmn, greenp(:,1:istore), istore, istart, ip,          &
     &                 ndim)
         END DO
!$OMP END DO
         DEALLOCATE(green, greenp)
!$OMP END PARALLEL
!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!
         CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, 
     1               ndim, ns)
         DEALLOCATE (gstore)

      ENDIF

      DEALLOCATE (grpmn)

!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP WORKSHARE
      amatrix(:) = amatsav(:)
!$OMP END WORKSHARE

!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV)
      IF (ivacskip .EQ. 0) THEN
!$OMP WORKSHARE
         bvecsav = bvec - bvecsav
!$OMP END WORKSHARE
      END IF
!$OMP END PARALLEL


      END SUBROUTINE scalpot
