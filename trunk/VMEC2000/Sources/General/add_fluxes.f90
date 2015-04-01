      SUBROUTINE add_fluxes(overg, bsupu, bsupv, lcurrent)
      USE vmec_main
      USE realspace, ONLY: wint, guu, guv, chip, phip, sigma_an 
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in)    :: overg
      REAL(rprec), DIMENSION(nrzt), INTENT(inout) :: bsupu, bsupv
      LOGICAL, INTENT(in) :: lcurrent
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
      REAL(rprec), PARAMETER :: iotaped = 0.10_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, ub, lb
      REAL(rprec), DIMENSION(2:ns) :: top, bot
      REAL(rprec) :: tops, bots
!-----------------------------------------------
!
!     ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
!     COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
!     IF ncurr == 1
!
      IF (.not.lcurrent .or. ncurr.eq.0) GOTO 100

      top = icurv(2:ns)
      bot = 0
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED) 
!$OMP DO PRIVATE(js,lb) 
      DO js = 2, ns
         DO l=1,nznt
         lb=js+ns*(l-1)
!         ub=ns*l
#ifdef _ANIMEC
         top(js) = top(js)-wint(lb)*(guu(lb)*bsupu(lb)                &
                               +guv(lb)*bsupv(lb))*sigma_an(lb)
         bot(js) = bot(js) + wint(lb)*sigma_an(lb)*overg(lb)*guu(lb)
#else
         top(js) = top(js) - wint(lb)*(guu(lb)*bsupu(lb)              &
             +                  guv(lb)*bsupv(lb))
         bot(js) = bot(js) + wint(lb)*overg(lb)*guu(lb)
#endif
         END DO
      END DO
!$OMP END DO
!$OMP WORKSHARE
      WHERE (bot .ne. zero) chips(2:ns) = top/bot
      WHERE (phips(2:ns) .ne. zero) iotas(2:ns) = chips(2:ns)/phips(2:ns)
!$OMP END WORKSHARE
!$OMP END PARALLEL

 100  CONTINUE

!     CHANGE THIS FOR lRFP = T  (solve for phips?)
      IF (ncurr .eq. 0) THEN
         chips = iotas*phips
      ELSE IF (.not.lcurrent) THEN
         WHERE (phips .ne. zero) iotas = chips/phips
      END IF
!!      IF (.not.lcurrent) chips = iotas*phips

      DO js = 2, ns
         chip(js:nrzt:ns) = chips(js)
      END DO

      chipf(2:ns1) = (chips(2:ns1) + chips(3:ns1+1))/2
      chipf(ns)    = 2*chips(ns)-chips(ns1)


!     Do not compute iota too near origin
      IF (lrfp) THEN
         iotaf(1)  = one/(c1p5/iotas(2) - p5/iotas(3))
         iotaf(ns) = one/(c1p5/iotas(ns) - p5/iotas(ns-1))
         DO js = 2, ns-1
            iotaf(js) = 2.0_dp/(one/iotas(js) + one/iotas(js+1))
         END DO

      ELSE
         iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)               !zero gradient near axis
         iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
         DO js = 2, ns-1
            iotaf(js) = p5*(iotas(js) + iotas(js+1))
         END DO
      END IF

!$OMP PARALLEL DO NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED) &
!$OMP PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l
         chip(lb+1:ub) = chips(2:ns)
         bsupu(lb:ub) = bsupu(lb:ub)+chip(lb:ub)*overg(lb:ub)
      END DO
!$OMP END PARALLEL DO

      END SUBROUTINE add_fluxes
