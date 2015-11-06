      SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv2), INTENT(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: bvec_temp1, bvec_temp2
      REAL(rprec) :: sin_temp, cos_temp
C-----------------------------------------------

      IF (n .lt. 0) STOP 'error calling analysum!'

      bvec_temp1 = zero
      bvec_temp2 = zero

      DO i = 1, nuv2
         sin_temp = (sinu1(i,m)*cosv1(i,n) - sinv1(i,n)*cosu1(i,m))            &
     &            * cmns(l,m,n) !SIN(mu - |n|v)*cmns
         bvec_temp1 = bvec_temp1 + (tl(i)*bexni(i)*sin_temp)

         IF (ivacskip .eq. 0) THEN
            grpmn(m,n,i,1) = grpmn(m,n,i,1) + sl(i)*sin_temp
         END IF

         IF (lasym) THEN
            cos_temp = (cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m))         &
     &               * cmns(l,m,n) !COS(mu - |n|v)*cmns
            bvec_temp2 = bvec_temp2 + tl(i)*bexni(i)*cos_temp

            IF (ivacskip .eq. 0) THEN
               grpmn(m,n,i,2) = grpmn(m,n,i,2) + sl(i)*cos_temp
            END IF
         END IF
      END DO

      bvec(m,n,1) = bvec(m,n,1) + bvec_temp1
      IF (lasym) THEN
         bvec(m,n,2) = bvec(m,n,2) + bvec_temp2
      END IF

      END SUBROUTINE analysum
