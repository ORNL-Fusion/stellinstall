      SUBROUTINE othersums(irho)
c--
c     calculates other1
c     of integral of W(lambda) which is not dependant on lambda
c--
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m
      REAL(rprec) :: denom, qn
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: sumit
C-----------------------------------------------
c
c-  m - poloidal, n - toroidal
c
      IF (isymm0 .ne. 0) THEN
         other1(irho) = zero
         RETURN
      ENDIF
c
c- Form the "other" SUMs.  The first USEs alpha1(m,n), calculated in WOFLAM,
c  and d(m,n) from DENM.
c
c- Load  rfmn with
c
c  (m*R+n*periods*S)/(m-n*periods*q) *
c         EXP(m*thetamax-n*zetahmax) * (1.5*alpha1(m,n)+d(m,n))
c
c- for only those harmonics that are going to be used in the sum.
c
      qn = periods*qsafety(irho)*zetasign
      DO m = -mbuse, mbuse
         DO n = 0, nbuse
            denom = m + n*qn
            IF (n.ne.0 .or. m.ne.0) THEN
               denom = denom/(denom**2 + (damp_bs*m)**2)
               rfmn(m,n) = (m*capr(irho)+n*periods*caps(irho))*denom*(
     1            COS(m*thetamax(irho)-n*zetahmax(irho))*REAL(1.5_dp*
     2            alpha1mn(m,n)+dmn(m,n))-SIN(m*thetamax(irho)-n*
     3            zetahmax(irho))*AIMAG(1.5_dp*alpha1mn(m,n)+dmn(m,n)))
            ELSE
               rfmn(m,n) = zero
            ENDIF
         END DO
      END DO

c  First SUM,

      other1(irho) = sumit(rfmn,mbuse,nbuse)
c
c- Then multiply by ALL the stuff out front.
c
      other1(irho) = -other1(irho)*qsafety(irho)/ftrapped(irho)*(one +
     1   aiogar(irho)/qsafety(irho))


      END SUBROUTINE othersums
