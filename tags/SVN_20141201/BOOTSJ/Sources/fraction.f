      SUBROUTINE fraction(irho)
c--
c  This calculates the fraction passing and the fraction trapped.
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
      INTEGER irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      INTEGER, PARAMETER :: n_lambda = 41
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l
      REAL(rprec), DIMENSION(n_lambda) ::
     1   alambda, avgvpov, antgrand, answer

C-----------------------------------------------

      dlambda = one/(n_lambda - 1)

c  Fill lambda mesh, evaluate V||/V

      DO l = 1, n_lambda
         alambda(l) = (l - 1)*dlambda
         avgvpov(l) = SUM(SQRT(ABS(one - alambda(l)*bfield))*gsqrt_b)
      END DO


c  Form integrand

      WHERE (avgvpov .gt. zero) antgrand = alambda/avgvpov

c  Do integral for the passing fraction

      CALL simpun (alambda, antgrand, n_lambda, answer)

      fpassing(irho) = .75_dp*avgbobm2*answer(n_lambda)*sum_gsqrt_b
      ftrapped(irho) = one - fpassing(irho)

      END SUBROUTINE fraction
