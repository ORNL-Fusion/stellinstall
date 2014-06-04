      SUBROUTINE do_fft(a11, answer_mn, trigsu, trigsv, ifaxu, ifaxv,
     1   ntheta, nzeta, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ntheta, nzeta, mbuse, nbuse
      INTEGER, DIMENSION(*) :: ifaxu, ifaxv
      REAL(rprec), DIMENSION(ntheta + 2,nzeta) :: a11
      REAL(rprec), DIMENSION(3*ntheta/2 + 1) :: trigsu
      REAL(rprec), DIMENSION(2*nzeta) :: trigsv
      complex(rprec), DIMENSION(-mbuse:mbuse,0:nbuse) :: answer_mn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: inc, jump, isign, jumpv, incv
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work11
      REAL(rprec) :: factor
C-----------------------------------------------

      ALLOCATE (work11(nzeta*(ntheta+2)))
c
c- Forward REAL to complex transform in the theta direction with INDEX n.
c  i.e., integral of EXP(-i*n*theta) * f(theta,zetah).
c
      inc = 1
      jump = ntheta + 2
      isign = -1

      CALL fft991(a11,work11,trigsu,ifaxu,inc,jump,ntheta,nzeta,isign)

c
c- now forward transform in the zetah direction with INDEX m.
c  i.e., integral of EXP(-i*m*zetah) * [theta transform of f(theta,zetah)]
c
      jumpv = 1
      incv = jump/2
      CALL cfft99(a11,work11,trigsv,ifaxv,incv,jumpv,nzeta,incv,isign)

      DEALLOCATE (work11)

c
c- Now reorganize and scale these to get the complex d(n,m)
c  factor = 1 / (number of points used in the forward transform direction).
c  because of (anti)symmetry, only nonnegative m values are needed.  we also
c  only fill to n = mbuse and m = nbuse since this is all that will be
c  used in the sums over n and m.
c
      factor = one/nzeta
c  store a11 in answer_mn array
      CALL reorganz (answer_mn, mbuse, nbuse, factor, a11, ntheta,
     1   nzeta)

      END SUBROUTINE do_fft
