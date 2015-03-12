!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

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
!-----------------------------------------------
!
!     ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
!     COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
!     IF ncurr == 1
!
      IF (.not.lcurrent .or. ncurr.eq.0) GOTO 100

      top = icurv(2:ns)
      bot = 0
!$OMP PARALLEL DO &
!$OMP PRIVATE(l,lb,ub)&
!$OMP REDUCTION(+:top,bot)
      DO l=1,nznt
         lb=ns*(l-1)+2
         ub=ns*l
#ifdef _ANIMEC
         top = top-wint(lb:ub)*(guu(lb:ub)*bsupu(lb:ub)                &
                               +guv(lb:ub)*bsupv(lb:ub))*sigma_an(lb:ub)
         bot = bot + wint(lb:ub)*sigma_an(lb:ub)*overg(lb:ub)*guu(lb:ub)
#else
         top = top - wint(lb:ub)*(guu(lb:ub)*bsupu(lb:ub)              &
             +                  guv(lb:ub)*bsupv(lb:ub))
         bot = bot + wint(lb:ub)*overg(lb:ub)*guu(lb:ub)
#endif
      END DO
!$OMP END PARALLEL DO
      WHERE (bot .ne. zero) chips(2:ns) = top/bot
      WHERE (phips(2:ns) .ne. zero) iotas(2:ns) = chips(2:ns)/phips(2:ns)

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

!$OMP PARALLEL DO &
!$OMP PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l
         chip(lb+1:ub) = chips(2:ns)
         bsupu(lb:ub) = bsupu(lb:ub)+chip(lb:ub)*overg(lb:ub)
      END DO
!$OMP END PARALLEL DO

      END SUBROUTINE add_fluxes
