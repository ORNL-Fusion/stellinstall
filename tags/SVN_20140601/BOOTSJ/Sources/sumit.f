      FUNCTION sumit (f, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER mbuse, nbuse
      REAL(rprec), DIMENSION(-mbuse:mbuse,0:nbuse) :: f
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m
      REAL(rprec) :: sumit
C-----------------------------------------------
c-
c specific SUM of terms in f(i,j) table
c-
c
      sumit = 0
      DO m = -mbuse, mbuse
         sumit = sumit + SUM(f(m,1:nbuse))
      END DO
      sumit = 2*sumit
      sumit = sumit + SUM(f(-mbuse:mbuse,0))
      sumit = sumit - f(0,0)

      END FUNCTION sumit
