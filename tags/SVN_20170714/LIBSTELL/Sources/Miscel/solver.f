      SUBROUTINE solver(amat, b, m, nrhs, info)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)  :: m, nrhs
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(inout) :: amat(m,m), b(m,nrhs)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, ALLOCATABLE :: ipiv(:)
!-----------------------------------------------
      info = 0
      ALLOCATE (ipiv(m))

!     Compute the solution to a REAL system of linear equations
!       AMAT * X = B,
!     where AMAT is an M-by-M matrix and X and B are N-by-NRHS matrices.
!       X is returned in B
!
!     FACTOR AMATRIX INTO LU FORM
!     AND SOLVE BY GAUSSIAN ELIMINATION
!
      CALL dgesv (m, nrhs, amat, m, ipiv, b, m, info)
!     IF (info .ne. 0) PRINT *, ' Condition No. = 0 in Solver'

      DEALLOCATE (ipiv)

      END SUBROUTINE solver
