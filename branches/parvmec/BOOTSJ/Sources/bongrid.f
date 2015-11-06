      SUBROUTINE bongrid(irho, ians)

c  First time through, form the THETA-ZETAH grid
c  Every time through, evaluate B and QSAFETY on the grid for this RHO.
c  Here, m is the poloidal mode number and n (and nh, for n-hat) is the
c  toroidal mode number.

C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE parambs
      USE vmec0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER irho, ians
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, nh, m, mbuse1, imax, jmax
      INTEGER, SAVE :: ihere = 0
      INTEGER :: ij_max(2)
      REAL(rprec), SAVE :: twopi, dth, dzetah, dthm, dzetahm
      REAL(rprec) :: b
C-----------------------------------------------


      IF (ihere .eq. 0) THEN

c  ALLOCATE ALL of the variables needed on the theta-zeta grid

         CALL allocate_angles

c-----------------------------------------------------------------------
c  Form the THETA-ZETAH grid.

         twopi = 8*ATAN(one)
         dth = twopi/nthetah
         dzetah = twopi/nzetah

         DO i = 1, nthetah
            theta(i) = (i - 1)*dth
         END DO
         DO j = 1, nzetah
            zetah(j) = (j - 1)*dzetah
         END DO

!   Form a finer mesh to evaluate fine grained bmax.

         dthm = dth/(nthetahm-1)
         dzetahm = dzetah/(nzetahm-1)

!  load sin and cosine arrays

         DO j = 1, nzetah
            DO nh = 0, nbuse
               sinnj(nh,j) = SIN(zetasign*nh*zetah(j))
               cosnj(nh,j) = COS(zetasign*nh*zetah(j))
            ENDDO
         ENDDO
         DO i = 1, nthetah
            DO m = -mbuse, mbuse
               sinmi(m,i) = SIN(m*theta(i))
               cosmi(m,i) = COS(m*theta(i))
            END DO
         END DO
         ihere = 1
      ENDIF          !end of ihere initial evaluation

c  Now evaluate B on the theta-zetah grid for this RHO using the EPSILONs just
c  found.  Loop over theta and phihat, summing up ALL the (m,nh) terms in B.

      DO j = 1, nzetah
         DO i = 1, nthetah
            b = zero
            DO m = -mbuse, mbuse
               DO nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1               (cosmi(m,i)*cosnj(nh,j)-sinmi(m,i)*sinnj(nh,j))
               END DO
            END DO
            bfield(i,j) = ABS(b)
         END DO
      END DO

c   find MAX of b on global mesh

      ij_max = MAXLOC(bfield)
      imax = ij_max(1)
      jmax = ij_max(2)

c  USE the theta and zeta from this search as the center of a finer search

c  first form the grid

      thetam(1) = theta(imax) - dth/2
      zetahm(1) = zetah(jmax) - dzetah/2

      DO i = 2, nthetahm
         thetam(i) = thetam(i-1) + dthm
      ENDDO
      DO j = 2, nzetahm
         zetahm(j) = zetahm(j-1) + dzetahm
      ENDDO

c  load the sines and cosines on the finer mesh

      DO j = 1, nzetahm
         DO nh = 0, nbuse
            sinnjm(nh,j) = SIN(zetasign*nh*zetahm(j))
            cosnjm(nh,j) = COS(zetasign*nh*zetahm(j))
         ENDDO
      ENDDO
      DO i = 1, nthetahm
         DO m = -mbuse, mbuse
            sinmim(m,i) = SIN(m*thetam(i))
            cosmim(m,i) = COS(m*thetam(i))
         END DO
      END DO

c  evaluate b on the finer mesh

      DO j = 1, nzetahm
         DO i = 1, nthetahm
            b = zero
            DO m = -mbuse, mbuse
               DO nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1              (cosmim(m,i)*cosnjm(nh,j)-sinmim(m,i)*sinnjm(nh,j))
               END DO
            END DO
            bfieldm(i,j) = ABS(b)
         END DO
      END DO


c- evaluate bmax1(irho), thetamax(irho), b2avg(irho), and zetahmax(irho)
c  based on finer mesh evaluation

      ij_max = MAXLOC(bfieldm)
      imax = ij_max(1)
      jmax = ij_max(2)
      thetamax(irho) = thetam(imax)
      zetahmax(irho) = zetahm(jmax)
      bmax1(irho) = bfieldm(imax,jmax)

c  evaluate jacobian.  Leave off flux surface quantites.  Evaluate
c  the SUM for later USE.  Note that the value is not scaled to
c  Bmax.

      gsqrt_b(:nthetah,:nzetah) = one/bfield(:nthetah,:nzetah)**2
      sum_gsqrt_b = SUM(gsqrt_b)

c  find b2avg LAB--boozer paper uses both bzero2 (1/leading coefficient of 1/b**2 expansion
c  and <b**2>.  They are the same (in boozer coordinates).  Both result from
c  flux surface averages. I will use <b2> everywhere

      b2avg(irho) = SUM(bfield**2 * gsqrt_b)/sum_gsqrt_b

c  Scale the array BFIELD so that it CONTAINS B/BMAX instead of B.

      bfield(:nthetah,:nzetah) = bfield(:nthetah,:nzetah)/bmax1(irho)
      WHERE(bfield .gt. one) bfield = one
      b2obm = bfield**2

c   pressure related derivatives are needed
c   first calculate difference denominators for pressure related derivatives
c   difference array has differences on the full mesh

c   USE a parabolic fit near rho = 0 to get derivatives at 1st half mesh

      IF(irho .eq. 1) THEN
         drho = (rhoar(2)**2 - rhoar(1)**2)/(2*rhoar(1))

c   USE slope of last two points at outer rho point

      ELSEIF(irho .eq. irup) THEN
         drho = 0.5_dp*(d_rho(irho)+d_rho(irho-1))

c  ALL other points

      ELSE
         drho = d_rho(irho) + 0.5_dp*(d_rho(irho+1)+d_rho(irho-1))
      ENDIF

c  evaluate Electron temperature gradients in Kev

      IF (irho .ne. 1 .and. irho .ne. irup) temperho1 =
     1   (tempe1(irho+1)-tempe1(irho-1))/drho
      IF (irho .eq. 1) temperho1 = (tempe1(irho+1)-tempe1(irho))/drho
      IF (irho .eq. irup) temperho1=(tempe1(irho)-tempe1(irho-1))/drho

c evaluate Ion temperature gradients in Kev

      IF (irho .ne. 1 .and. irho .ne. irup) tempirho1 =
     1  (tempi1(irho+1)-tempi1(irho-1))/drho
      IF (irho .eq. 1)tempirho1 = (tempi1(irho+1)-tempi1(irho))/drho
      IF (irho .eq. irup)temperho1=(tempi1(irho)-tempi1(irho-1))/drho

c  evaluate electron density gradients in 10**20 m-3

      IF (irho .ne. 1 .and. irho .ne. irup) densrho1 =
     1  (dense(irho+1)-dense(irho-1))/drho
      IF (irho .eq. 1) densrho1 = (dense(irho+1)-dense(irho))/drho
      IF (irho .eq. irup) densrho1 = (dense(irho)-dense(irho-1))/drho

c  WRITE out the lower order Boozer coefficients

      mbuse1 = mbuse
      IF (mbuse > 5) mbuse1 = 5            !mbuse1 > 5 will not fit on page

      WRITE (ians, 400, advance='no')
  400 FORMAT(/' nh ')
      DO m = -mbuse1, mbuse1
         WRITE (ians, 402, advance='no') m
      END DO
  402 FORMAT('   m=',i2,'   ')
      WRITE (ians, '(a1)') ' '
c
      DO nh = 0, nbuse
         WRITE (ians, 406) nh, (amnfit(irho,m,nh),m=(-mbuse1),mbuse1)
      END DO
  406 FORMAT(1x,i2,1p,13e10.3)

c  Calculate the fraction trapped and fraction passing for the "equivalent"

      CALL tok_fraction(fttok(irho))

      fptok(irho) = one - fttok(irho)

      END SUBROUTINE bongrid
