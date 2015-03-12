!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

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
      INTEGER :: m,n,js,l
      REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor
!-----------------------------------------------
      blam=0
      clam=0
      dlam=0
!$OMP PARALLEL DO &
!$OMP PRIVATE(l) &
!$OMP REDUCTION(+:blam,clam,dlam)
      DO l=1,nznt
         blam(1:ns) = blam(1:ns) + guu(:,l)*overg(:,l)
         clam(1:ns) = clam(1:ns) + gvv(:,l)*overg(:,l)
         dlam(1:ns) = dlam(1:ns) + guv(:,l)*overg(:,l)
      END DO
!$OMP END PARALLEL DO
      blam(1) = blam(2)
      clam(1) = clam(2)
      dlam(1) = dlam(2)
      DO js = 2, ns
        blam(js) = cp5*(blam(js) + blam(js+1))
        clam(js) = cp5*(clam(js) + clam(js+1))
        dlam(js) = cp5*(dlam(js) + dlam(js+1))
      END DO

      faclam = 0
      pfactor0 = damping_fac/(2*r0scale*lamscale)**2

!$OMP PARALLEL DO &
!$OMP PRIVATE(m,tmm,power,pfactor,n,tnn,tnm,js) 
      DO m = 0, mpol1
         tmm = m*m
         power = MIN(tmm/256, 8._dp)
         pfactor = pfactor0
         DO n = 0, ntor
            IF (m.eq.0 .and. n.eq.0) CYCLE
!            IF (n .gt. 1) pfactor = pfactor0/4                          !sometimes helps convergence
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            DO js = jlam(m), ns
               faclam(js,n,m,1) = (blam(js)*tnn                         &
               + SIGN(dlam(js),blam(js))*tnm                            &
               +       clam(js)*tmm)
              IF (faclam(js,n,m,1) .eq. zero)                           &
                  faclam(js,n,m,1) = -1.E-10_dp
            faclam(js,n,m,1) = (pfactor/faclam(js,n,m,1))               &
               * sqrts(js)**power                                        !Damps m > 16 modes
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO 

      DO n = 2, ntmax
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
      END DO

!
!     ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
!
      DO js = 1, ns
         faclam(js,0,0,1) = (pfactor*lamscale**2)/blam(js)
      END DO

      END SUBROUTINE lamcal
