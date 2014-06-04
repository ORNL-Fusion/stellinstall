      SUBROUTINE caprsh2(irho)
c--
c  Evaluate CAPR, CAPS, and H2.  Here, m is the poloidal mode number
c  and n is the toroidal mode number/periods.
c--LAB--changed calculation of H2 to USE proper 1/b**2 Jacobian for
c  flux surface average
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
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nh, m, i, j
      REAL(rprec) :: h2top_th, h2top_phi, den, qn
      REAL(rprec) ::  b2, h2top
      REAL(rprec) :: den_tmp, sin_eps
C-----------------------------------------------------------------------

c  the following code was added to evaluate H2 LAB
      h2top = zero
      den   = zero
      DO i=1, nthetah
         DO j=1, nzetah
            h2top_th = zero
            h2top_phi = zero
            den_tmp = zero
            DO nh = 0, nbuse            !  evaluate i,j terms in num and denom of H2
               qn = qsafety(irho)*periods*nh*zetasign
               DO m = -mbuse, mbuse
                  sin_eps = amnfit(irho,m,nh)*
     1            (sinmi(m,i)*cosnj(nh,j)+cosmi(m,i)*sinnj(nh,j))
                  h2top_th = h2top_th - m*sin_eps
                  h2top_phi = h2top_phi - qn*sin_eps
                  den_tmp = den_tmp - (m + qn)*sin_eps
               ENDDO
            ENDDO
            b2 = bfield(i,j)**2
            h2top = h2top + (h2top_th**2 - h2top_phi**2)/b2
            den  = den + den_tmp**2/b2
         ENDDO
      ENDDO

      IF (den .eq. zero) STOP 'den = 0 in caprsh2'
      h2(irho) = h2top/den
      capr(irho) = (one - h2(irho))/(2*qsafety(irho))
      caps(irho) = (one + h2(irho))/2

      END SUBROUTINE caprsh2
