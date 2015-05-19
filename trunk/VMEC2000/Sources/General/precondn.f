      SUBROUTINE precondn(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1                    xodd, axm, axd, bxm, bxd, 
#ifdef _HBANGLE
     2                    axm2, axd2, bxm2, bxd2,
#endif
     3                    cx, eqfactor, trigmult)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
#ifdef _HBANGLE
      USE angle_constraints, ONLY: cos2u, sin2u
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm, axd, bxm, bxd
#ifdef _HBANGLE
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm2, axd2, bxm2, bxd2
#endif
      REAL(rprec), DIMENSION(ns+1), INTENT(out) :: cx
      REAL(rprec), DIMENSION(ns), INTENT(out) :: eqfactor
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: trigmult
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, lb, lb0
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
#ifdef _HBANGLE
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax2, bx2
#endif
      REAL(rprec) :: pfactor, temp(ns+1), temp2(ns), cx2(ns)
      REAL(rprec):: t1, t2, t3, ptau, ptau2, osh
C-----------------------------------------------
!
!     COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z
!     FORCE. NOTE THAT THE NORMALIZATION IS:
!
!     AX(off-diag) ~ <(cosmui cosmu cosnv cosnv) 2(R**2*Xu**2*bsq/gsqrt)>
!     Factor of 2 arising from 1/gsqrt**2 in bsq
!
!     Now, cosmui cosmu ~ mscale(0)**2, cosnv**2 ~ nscale(0)**2
!     Therefore, AX ~ (mscale(0)*nscale(0))**2 2<R**2*Xu**2*bsq/gsqrt>
!                   ~ 2*r0scale**2 <...>
!
      ALLOCATE (ax(ns+1,4), bx(ns+1,4))
      ax = 0;  bx = 0;  cx = 0;  temp = 0
!      pfactor = -2*r0scale**2       !v8.50
      pfactor = -4*r0scale**2        !restored in v8.51
#ifdef _HBANGLE
      ALLOCATE (ax2(ns+1,4), bx2(ns+1,4))
      ax2 = 0; bx2 = 0
#endif
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP DO PRIVATE(js,l,lb,lb0,t1,t2,t3,osh,ptau,ptau2)
      DO js = 2,ns
         osh = one/shalf(js)
      DO l=1,nznt
         lb0=ns*(l-1) 
         lb=lb0+js 
         t1 = pfactor*r12(lb)*bsq(lb)
         ptau2 = r12(lb)*t1/gsqrt(lb)
         t1 = t1*wint(lb)
         temp(js) = temp(js) + t1*trigmult(l)*xu12(lb)
         ptau = r12(lb)*t1/gsqrt(lb)
         t1 = xu12(lb)*ohs
         t2 = cp25*(xue(lb)*osh + xuo(lb))*osh
         t3 = cp25*(xue(lb-1)*osh + xuo(lb-1))*osh
         ax(js,1) = ax(js,1) + ptau*t1*t1
         ax(js,2) = ax(js,2) + ptau*(-t1+t3)*(t1+t2)
         ax(js,3) = ax(js,3) + ptau*(t1+t2)*(t1+t2)
         ax(js,4) = ax(js,4) + ptau*(-t1+t3)*(-t1+t3)
#ifdef _HBANGLE
         ptau2 = ptau2*cos2u(l)
         ax2(js,1) = ax2(js,1) + ptau2*t1*t1
         ax2(js,2) = ax2(js,2) + ptau2*(-t1+t3)*(t1+t2)
         ax2(js,3) = ax2(js,3) + ptau2*(t1+t2)*(t1+t2)
         ax2(js,4) = ax2(js,4) + ptau2*(-t1+t3)*(-t1+t3)
#endif
!
!       COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR M**2, N**2 TERMS
!
         t1 = cp5*(xs(lb) + cp5*xodd(lb)*osh)
         t2 = cp5*(xs(lb) + cp5*xodd(lb-1)*osh)
         bx(js,1) = bx(js,1) + ptau*t1*t2
         bx(js,2) = bx(js,2) + ptau*t1*t1
         bx(js,3) = bx(js,3) + ptau*t2*t2
         cx(js) = cx(js) + cp25*pfactor*lu1(lb)**2*gsqrt(lb)
     1                   * wint(lb)
#ifdef _HBANGLE
         bx2(js,1) = bx2(js,1) + ptau2*t1*t2
         bx2(js,2) = bx2(js,2) + ptau2*t1*t1
         bx2(js,3) = bx2(js,3) + ptau2*t2*t2
#endif
        END DO
      END DO
!$OMP END DO

!$OMP WORKSHARE
      temp = temp/vp
!$OMP END WORKSHARE

!$OMP DO PRIVATE(js)
      DO js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
#ifdef _HBANGLE
        axm2(js,1) =-ax2(js,1)
        axd2(js,1) = ax2(js,1) + ax2(js+1,1)
        axm2(js,2) = ax2(js,2) * sm(js) * sp(js-1)
        axd2(js,2) = ax2(js,3)*sm(js)**2 + ax2(js+1,4)*sp(js)**2
        bxm2(js,1) = bx2(js,1)
        bxm2(js,2) = bx2(js,1) * sm(js) * sp(js-1)
        bxd2(js,1) = bx2(js,2) + bx2(js+1,3)
        bxd2(js,2) = bx2(js,2)*sm(js)**2 + bx2(js+1,3)*sp(js)**2
#endif
        cx2(js)    = cx(js) + cx(js+1)
        temp2(js)  = signgs*(temp(js) + temp(js+1))
        eqfactor(js) = axd(js,2)*hs*hs/temp2(js)
      END DO
!$OMP END DO
!$OMP WORKSHARE
      cx = cx2
!$OMP END WORKSHARE
!$OMP END PARALLEL

      eqfactor(1) = 0;  eqfactor(ns) = 0
      axm(ns+1,:) = 0;  axd(ns+1,:) = 0
      bxm(ns+1,:) = 0;  bxd(ns+1,:) = 0
      DEALLOCATE (ax, bx)
#ifdef _HBANGLE
      axm2(ns+1,:) = 0;  axd2(ns+1,:) = 0
      bxm2(ns+1,:) = 0;  bxd2(ns+1,:) = 0
      DEALLOCATE (ax2, bx2)
#endif
      END SUBROUTINE precondn
