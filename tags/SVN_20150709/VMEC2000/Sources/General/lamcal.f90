      SUBROUTINE lamcal(overg, guu, guv, gvv)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, jlam, lamscale
      USE realspace, ONLY: sqrts, phip
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::                    &
         overg, guu, guv, gvv
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: damping_fac=2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m,n,js
      REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor
!-----------------------------------------------

      pfactor0 = damping_fac/(2*r0scale*lamscale)**2
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP WORKSHARE
      blam=0
      clam=0
      dlam=0
      faclam = 0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(js)
      DO js = 1, ns
         blam(js) = SUM(guu(js,:)*overg(js,:))
         clam(js) = SUM(gvv(js,:)*overg(js,:))
         dlam(js) = SUM(guv(js,:)*overg(js,:))
!         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
 
      blam(1) = blam(2)
      clam(1) = clam(2)
      dlam(1) = dlam(2)

      DO js = 2, ns
        blam(js) = cp5*(blam(js) + blam(js+1))
        clam(js) = cp5*(clam(js) + clam(js+1))
        dlam(js) = cp5*(dlam(js) + dlam(js+1))
      END DO


!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP DO PRIVATE(m,n,tmm,power,pfactor,tnn,tnm,js) COLLAPSE(2)
      DO m = 0, mpol1
         DO n = 0, ntor
         tmm = m*m
         power = MIN(tmm/256, 8._dp)
         pfactor = pfactor0
            IF (m.EQ.0 .AND. n.EQ.0) CYCLE
!            IF (n .gt. 1) pfactor = pfactor0/4                          !sometimes helps convergence
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            DO js = jlam(m), ns
               faclam(js,n,m,1) = (blam(js)*tnn                         &
               + SIGN(dlam(js),blam(js))*tnm                            &
               +       clam(js)*tmm)
              IF (faclam(js,n,m,1) .EQ. zero)                           &
                  faclam(js,n,m,1) = -1.E-10_dp
            faclam(js,n,m,1) = (pfactor/faclam(js,n,m,1))               &
               * sqrts(js)**power                                        !Damps m > 16 modes
            END DO
         END DO
      END DO
!$OMP END DO

      DO n = 2, ntmax
!$OMP WORKSHARE
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
!$OMP END WORKSHARE
      END DO

!
!     ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
!
!$OMP DO PRIVATE(js)
      DO js = 1, ns
         faclam(js,0,0,1) = (pfactor0*lamscale**2)/blam(js)
      END DO
!$OMP END DO
!$OMP END PARALLEL

      END SUBROUTINE lamcal
