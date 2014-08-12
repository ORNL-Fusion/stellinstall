      SUBROUTINE tok_fraction(ft_tok)
c--
c  This calculates the fraction passing and the fraction trapped
c  for a tokamak using the Lin-Liu Miller approximation with Houlberg
c  weighting factors.  (Phys. Plas. 2, 1966 (1995).
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
      REAL(rprec) :: ft_tok
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i
      REAL(rprec) fub, flb, bavg, b2av, sum_gsqrt_tok
      REAL(rprec) bmax
C-----------------------------------------------

      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    bfield_tok,  gsqrt_tok, b2obm_tok, one_b2
      ALLOCATE(bfield_tok(nthetah),  gsqrt_tok(nthetah),
     1   b2obm_tok(nthetah), one_b2(nthetah), stat = i)
      IF(i .ne. 0) STOP 'allocation error of tokamak fields'


c  make symetric copies of field quantites--average over zeta to extract
c  symmetric component.

      DO i =1,nthetah
         bfield_tok(i) = SUM(bfield(i,:nzetah))/nzetah
      ENDDO

c  renormalize to MAX tok field equivalent

      bmax = MAXVAL(bfield_tok)
      bfield_tok = bfield_tok/bmax
      WHERE(bfield_tok .gt. one) bfield_tok = one

c  calculate 1/b**2, b**2

      one_b2 = one/bfield_tok**2
      b2obm_tok = bfield_tok**2

c  jacobian only includes 1/b**2 component and is normalized to bmax tok
c  integrate over jacobian to obtain normalization for averages
      sum_gsqrt_tok= SUM(one_b2)


c  find average b, b**2

      bavg = SUM(bfield_tok*one_b2)/sum_gsqrt_tok
      b2av = SUM(b2obm_tok*one_b2)/sum_gsqrt_tok

c  find upper bound

      fub = one-(one-SQRT(one-bavg)*(one+0.5_dp*bavg))*b2av/bavg**2

c  find lower bound


c  minus <1/b**2>

      flb = - SUM(one/b2obm_tok**2)/sum_gsqrt_tok

c  plus <SQRT(1-b)/b**2>

      flb = flb + SUM(SQRT(one-bfield_tok)
     1    /b2obm_tok**2)/sum_gsqrt_tok

c  plus <SQRT(1-b)/b)>/2

      flb = flb + 0.5_dp*SUM(SQRT(one-bfield_tok)
     1    /bfield_tok/b2obm_tok)/sum_gsqrt_tok

c  1+(previous stuff)*<b**2>


      flb = one + flb*b2av

c finally combine upper and lower limits with Houlberg factors
      ft_tok = 0.25_dp*flb + 0.75_dp*fub

      DEALLOCATE(bfield_tok, gsqrt_tok, b2obm_tok, one_b2)

      RETURN
      END SUBROUTINE tok_fraction
