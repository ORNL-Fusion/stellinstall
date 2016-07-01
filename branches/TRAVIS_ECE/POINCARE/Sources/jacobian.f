      SUBROUTINE jacobian(neq, phi, yin, ml, mu, pd, nrowd)                   ! dummy jacobian subroutine (not currently used; only
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE
      INTEGER:: neq, ml, mu, nrowd
      REAL(KIND=rprec):: phi
      REAL(KIND=rprec), DIMENSION(neq):: yin
      REAL(KIND=rprec), DIMENSION(nrowd, neq) :: pd
      END SUBROUTINE jacobian
