      SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv2), INTENT(in) :: sl, tl
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: bvec_
      REAL(rprec) :: sinp, cosp
!-----------------------------------------------
      IF (n .LT. 0) STOP 'error calling analysum!'
      
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP& PRIVATE(bvec_)
      ALLOCATE(bvec_(0:mf,-nf:nf,ndim), stat=i)
      IF (i. NE. 0) STOP 'Allocation error in analysum'
      bvec_ = 0
!$OMP DO PRIVATE(i,sinp,cosp) SCHEDULE(STATIC) 
      DO i = 1,nuv2
         sinp = (sinu1(i,m)*cosv1(i,n) - sinv1(i,n)*cosu1(i,m))
     1            *  cmns(l,m,n)                            !SIN(mu - |n|v)*cmns
         IF (ivacskip .EQ. 0) grpmn(m,n,i,1) = grpmn(m,n,i,1) 
     1                                        + sl(i)*sinp
         bvec_(m,n,1) = bvec_(m,n,1) + tl(i)*bexni(i)*sinp

         IF (lasym) THEN
            cosp = (cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m))
     1           *  cmns(l,m,n)                                !COS(mu - |n|v)*cmns
            IF (ivacskip .EQ. 0) grpmn(m,n,i,2) = grpmn(m,n,i,2) 
     1                                           + sl(i)*cosp
            bvec_(m,n,2) = bvec_(m,n,2) + tl(i)*bexni(i)*cosp
         END IF
      END DO
!$OMP END DO NOWAIT
!$OMP CRITICAL
      bvec  = bvec  + bvec_
!$OMP END CRITICAL
      DEALLOCATE(bvec_)
!$OMP END PARALLEL


      END SUBROUTINE analysum
