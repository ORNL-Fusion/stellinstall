      SUBROUTINE reorganz(coefs, mbuse, nbuse, factor, a1,
     1   ntheta, nzeta)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER mbuse, nbuse, ntheta, nzeta
      REAL(rprec) factor
      REAL(rprec), DIMENSION(ntheta + 2,nzeta) :: a1
      complex(rprec), DIMENSION(-mbuse:mbuse,0:nbuse) :: coefs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, m, i, n, ia
C-----------------------------------------------
c--
c  Version reorganz creates a COMPLEX coefs array.
c  Subroutine to reorganize and scale the output from fft991 and cft99 into complex
c  coefficients with INDEX (m,n).  The coefficients are scaled by the
c  factor FACTOR.  This factor should be 1/(number of points in the
c  forward complex transform direction).  The original coefficients appear
c  in array a1.  The scaled coefficients will be
c  copied to array COEFS in the calling list.  Here, m is the poloidal
c  mode number and n is the toroidal mode number/periods. Because we USE two
c  FORWARD transforms, we must flip the sign of m to get the desired nu
c  argument (u=theta, v=zeta)
c--
c
c- Because of (anti)symmetry, ONLY nonnegative values of n are needed.  We also
c  only fill to m = mbuse and n = nbuse since this is all that will be
c  used in the sums over m and n.
c  therefore, only i = 1 to mbuse+1 (for m = 0 to mbuse) and i = nth+1-mbuse
c  to nth (for m = -mbuse to -1) and ONLY j = 1 to nbuse+1 (for n = 0 to nbuse)
c  are needed.
c
      DO j = 1, nbuse + 1
         n = j - 1
         DO i = 1, mbuse + 1
            m = i - 1
            ia = 2*i - 1
            coefs(-m,n) = factor*CMPLX(a1(ia,j),a1(ia+1,j))
         END DO
      END DO

      DO j = nzeta + 1 - nbuse, nzeta
         n = -(nzeta + 1) + j
         DO i = 1, mbuse + 1
            m = i - 1
            ia = 2*i - 1
            coefs(m,-n) = factor*CMPLX(a1(ia,j),(-a1(ia+1,j)))
         END DO
      END DO

      coefs(1:mbuse,0) = coefs(-1:-mbuse:-1,0)

      END SUBROUTINE reorganz
