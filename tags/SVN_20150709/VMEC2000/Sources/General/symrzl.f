      SUBROUTINE symrzl(r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons,
     1   zcons, r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(inout) ::
     1   r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, zcons
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) ::
     1   r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mpar, ir, i, jk, jka, n2
C-----------------------------------------------
!
!     FIRST SUM SYMMETRIC, ANTISYMMETRIC PIECES ON EXTENDED INTERVAL, THETA = [PI,2*PI]
!
      DO mpar = 0, 1
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP DO PRIVATE(i,ir,jk,jka)
         DO i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i                 !-theta
            DO jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               r1s(jk,i,mpar) = r1s(jka,ir,mpar) - r1a(jka,ir,mpar)
               rus(jk,i,mpar) =-rus(jka,ir,mpar) + rua(jka,ir,mpar)
               z1s(jk,i,mpar) =-z1s(jka,ir,mpar) + z1a(jka,ir,mpar)
               zus(jk,i,mpar) = zus(jka,ir,mpar) - zua(jka,ir,mpar)
               lus(jk,i,mpar) = lus(jka,ir,mpar) - lua(jka,ir,mpar)
               rcons(jk,i,mpar)= rcons(jka,ir,mpar)-rcona(jka,ir,mpar)
               zcons(jk,i,mpar)=-zcons(jka,ir,mpar)+zcona(jka,ir,mpar)
            END DO
            IF (lthreed) THEN
               DO jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  rvs(jk,i,mpar)=(-rvs(jka,ir,mpar))+rva(jka,ir,mpar)
                  zvs(jk,i,mpar) = zvs(jka,ir,mpar) - zva(jka,ir,mpar)
                  lvs(jk,i,mpar) = lvs(jka,ir,mpar) - lva(jka,ir,mpar)
               END DO
            ENDIF
         END DO
!$OMP END DO
!
!        NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES FOR THETA = [0,PI]
!
!$OMP DO PRIVATE(i,jk)
         DO jk = 1, ns*nzeta
         DO i = 1, ntheta2
         r1s(jk,i,mpar) = r1s(jk,i,mpar) + r1a(jk,i,mpar)
         rus(jk,i,mpar) = rus(jk,i,mpar) + rua(jk,i,mpar)
         z1s(jk,i,mpar) = z1s(jk,i,mpar) + z1a(jk,i,mpar)
         zus(jk,i,mpar) = zus(jk,i,mpar) + zua(jk,i,mpar)
         lus(jk,i,mpar) = lus(jk,i,mpar) + lua(jk,i,mpar)
         rcons(jk,i,mpar) = rcons(jk,i,mpar) + rcona(jk,i,mpar)
         zcons(jk,i,mpar) = zcons(jk,i,mpar) + zcona(jk,i,mpar)
         IF (lthreed) THEN
            rvs(jk,i,mpar) = rvs(jk,i,mpar) + rva(jk,i,mpar)
            zvs(jk,i,mpar) = zvs(jk,i,mpar) + zva(jk,i,mpar)
            lvs(jk,i,mpar) = lvs(jk,i,mpar) + lva(jk,i,mpar)
         ENDIF
         END DO
         END DO
!$OMP END DO
!$OMP END PARALLEL
      END DO
      END SUBROUTINE symrzl
