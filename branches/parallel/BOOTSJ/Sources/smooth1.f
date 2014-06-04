      SUBROUTINE smooth1(ya, n1, n2, wk, frac)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n1, n2
      REAL(rprec) frac
      REAL(rprec), DIMENSION(*) :: ya, wk
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n11, n21, n, i
      REAL(rprec) :: as, a, a1, a2, a3
C-----------------------------------------------
c--
c       smoothing REAL array ya(*) for i: n1.le.i.le.n2
c     frac - defines how strong is smoothing;
c     IF frac=0 THEN ONLY smoothes the peaks
c--
      n11 = n1 + 1
      n21 = n2 - 1
      n = n21 - n11
      IF (n .le. 0) RETURN                    !too little number of points
c
c- check the edges
c
      as = SUM(ABS(ya(n11:n21-1)-ya(n11+1:n21)))
      as = as/n
c
      IF (as .le. zero) THEN
         ya(n1) = ya(n11)
         ya(n2) = ya(n21)
         RETURN
      ELSE
         a = 3*(as + ABS(ya(n11)-ya(n11+1)))
         a1 = ABS(ya(n1)-ya(n11))
         IF (a1 > a) ya(n1) = 2*ya(n11) - ya(n11+1)
         a = 3*(as + ABS(ya(n21)-ya(n21-1)))
         a1 = ABS(ya(n2)-ya(n21))
         IF (a1 > a) ya(n2) = 2*ya(n21) - ya(n21-1)
      ENDIF
c
c- work array
c
      wk(n1:n2) = ya(n1:n2)
c
c- check for strong peaks and remove
c
      DO i = n11, n21
         a1 = .5_dp*(ya(i+1)+ya(i-1))
         a2 = 3*ABS(ya(i+1)-ya(i-1))
         a3 = ABS(ya(i)-a1)
         IF (a3 > a2) ya(i) = a1
      END DO
c
c- smoothing with a factor frac
c
      IF (frac .le. zero) RETURN
      ya(n11:n21) =(ya(n11:n21)+.5_dp*(wk(n11+1:n21+1)+wk(n11-1:n21-1))*
     1   frac)/(one + frac)

      END SUBROUTINE smooth1
