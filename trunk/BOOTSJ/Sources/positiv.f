      SUBROUTINE positiv(ya, n, ivar)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, ivar
      REAL(rprec), DIMENSION(*) :: ya
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       REAL(rprec), PARAMETER :: zero = 0, D36 = 1.E-36_dp
C-----------------------------------------------
c--
c     make corrections to ensure a positive array
c--
c
      IF (ivar .eq. 1) THEN
         ya(:n) = ABS(ya(:n)) + D36
      ELSE
         WHERE (ya(:n) .le. zero) ya(:n) = D36
      ENDIF

      END SUBROUTINE positiv
