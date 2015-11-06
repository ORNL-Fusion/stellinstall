      SUBROUTINE analysum2(grpmn, bvec, slp, tlp, slm, tlm,
     1    m, n, l, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv2), INTENT(in) ::
     1   slp, tlp, slm, tlm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: bvec_temp_p1, bvec_temp_p2
      REAL(rprec) :: bvec_temp_n1, bvec_temp_n2
      REAL(rprec) :: sin_temp, cos_temp, temp
      REAL(rprec) :: sin_temp_m, cos_temp_m
C-----------------------------------------------

      IF (n .lt. 0) STOP 'error calling analysum2!'

      bvec_temp_p1 = zero
      bvec_temp_p2 = zero
      bvec_temp_n1 = zero
      bvec_temp_n2 = zero

      DO i = 1, nuv2

         sin_temp = sinu1(i,m)*cosv1(i,n)*cmns(l,m,n)
         temp = -cosu1(i,m)*sinv1(i,n)*cmns(l,m,n)
         sin_temp_m = sin_temp - temp !SIN(mu + |n|v) * cmns (l,m,|n|)
         sin_temp = sin_temp + temp   !SIN(mu - |n|v) * cmns (l,m,|n|)
         bvec_temp_p1 = bvec_temp_p1 + tlp(i)*bexni(i)*sin_temp
         bvec_temp_n1 = bvec_temp_n1 + tlm(i)*bexni(i)*sin_temp_m

         IF (ivacskip .eq. 0) THEN
            grpmn(m,n,i,1)  = grpmn(m,n,i,1)  + slp(i)*sin_temp
            grpmn(m,-n,i,1) = grpmn(m,-n,i,1) + slm(i)*sin_temp_m
         END IF


         IF (lasym) THEN
            cos_temp = cosu1(i,m)*cosv1(i,n)*cmns(l,m,n)
            temp = sinu1(i,m)*sinv1(i,n)*cmns(l,m,n)
            cos_temp_m = cos_temp - temp !COS(mu + |n|v) * cmns (l,m,|n|)
            cos_temp = cos_temp + temp   !COS(mu - |n|v) * cmns (l,m,|n|)
            bvec_temp_p2 = bvec_temp_p2 + tlp(i)*bexni(i)*cos_temp
            bvec_temp_n2 = bvec_temp_n2 + tlm(i)*bexni(i)*cos_temp_m

            IF (ivacskip .eq. 0) THEN
               grpmn(m,n,i,2)  = grpmn(m,n,i,2)  + slp(i)*cos_temp
               grpmn(m,-n,i,2) = grpmn(m,-n,i,2) + slm(i)*cos_temp_m
            END IF

         END IF
      END DO

      bvec(m,n,1)  = bvec(m,n,1)  + bvec_temp_p1
      bvec(m,-n,1) = bvec(m,-n,1) + bvec_temp_n1
      IF (lasym) THEN
         bvec(m,n,2)  = bvec(m,n,2)  + bvec_temp_p2
         bvec(m,-n,2) = bvec(m,-n,2) + bvec_temp_n2
      END IF

      END SUBROUTINE analysum2
