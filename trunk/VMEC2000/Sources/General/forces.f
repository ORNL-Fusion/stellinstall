!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

      SUBROUTINE forces
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces, ru12 => azmn_e, zu12 => armn_e, 
     1             azmn_e => azmn_e, armn_e => armn_e,
     2             lv_e => crmn_e, lu_e => czmn_e, lu_o => czmn_o,
     3             crmn_e => crmn_e, czmn_e => czmn_e, czmn_o => czmn_o
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p25 = p5*p5, dshalfds=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim, lb, ub
      REAL(rprec), DIMENSION(:), POINTER :: 
     1    bsqr, gvvs, guvs, guus
C-----------------------------------------------
      ndim = 1+nrzt
!
!     POINTER ALIASES
!
      bsqr => extra1(:,1);  gvvs => extra2(:,1)
      guvs => extra3(:,1);  guus => extra4(:,1)

!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT RADIAL-DERIVATIVE TERMS IN JACOBIAN
!
!     LV :  VARIATION OF R-TERM IN JACOBIAN
!
!     GVV:  VARIATION OF R**2-TERM AND Rv**2,Zv**2 IN gvv
!
!     GUU, GUV: VARIATION OF Ru, Rv, Zu, Zv IN guu, guv
!
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC)
!$OMP& PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l

         lu_e(lb) = 0; lv_e(lb) = 0
         guu(lb)  = 0; guv(lb)  = 0; gvv(lb) = 0
         guus(lb:ub) = guu(lb:ub)*shalf(1:ns)
         gvvs(lb:ub) = gvv(lb:ub)*shalf(1:ns)

         armn_e(lb:ub)  = ohs*zu12(lb:ub) * lu_e(lb:ub)
         azmn_e(lb:ub)  =-ohs*ru12(lb:ub) * lu_e(lb:ub)
         brmn_e(lb:ub)  = brmn_e(lb:ub) * lu_e(lb:ub)
         bzmn_e(lb:ub)  =-bzmn_e(lb:ub) * lu_e(lb:ub)
         bsqr(lb:ub)    = dshalfds*lu_e(lb:ub)/shalf(1:ns)
         armn_o(lb:ub)  = armn_e(lb:ub) *shalf(1:ns)
         azmn_o(lb:ub)  = azmn_e(lb:ub) *shalf(1:ns)
         brmn_o(lb:ub)  = brmn_e(lb:ub) *shalf(1:ns)
         bzmn_o(lb:ub)  = bzmn_e(lb:ub) *shalf(1:ns)
!
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
         guu(lb:ub-1) = p5*(guu(lb:ub-1) + guu(lb+1:ub))
         guu(ub) = p5*guu(ub)
         gvv(lb:ub-1) = p5*(gvv(lb:ub-1) + gvv(lb+1:ub))
         gvv(ub) = p5*gvv(ub)
         bsqr(lb:ub-1) = bsqr(lb:ub-1) + bsqr(lb+1:ub)
         guus(lb:ub-1) = p5*(guus(lb:ub-1) + guus(lb+1:ub))
         guus(ub) = p5*guus(ub)
         gvvs(lb:ub-1) = p5*(gvvs(lb:ub-1) + gvvs(lb+1:ub))
         gvvs(ub) = p5*gvvs(ub)

         armn_e(lb:ub-1) = armn_e(lb+1:ub) - armn_e(lb:ub-1) 
     1                   + p5*(lv_e(lb:ub-1) + lv_e(lb+1:ub))
         armn_e(ub) = -armn_e(ub)+p5*lv_e(ub)
         azmn_e(lb:ub-1) = azmn_e(lb+1:ub) - azmn_e(lb:ub-1)
         azmn_e(ub) = -azmn_e(ub)
         brmn_e(lb:ub-1) = p5*(brmn_e(lb:ub-1) + brmn_e(lb+1:ub))
         brmn_e(ub) = p5*brmn_e(ub)
         bzmn_e(lb:ub-1) = p5*(bzmn_e(lb:ub-1) + bzmn_e(lb+1:ub))
         bzmn_e(ub) = p5*bzmn_e(ub)

         armn_e(lb:ub) = armn_e(lb:ub) - (gvvs(lb:ub)*r1(lb:ub,1)
     1              + gvv(lb:ub)*r1(lb:ub,0))
         brmn_e(lb:ub) = brmn_e(lb:ub) + bsqr(lb:ub)*z1(lb:ub,1)
     1              -(guus(lb:ub)*ru(lb:ub,1) + guu(lb:ub)*ru(lb:ub,0))
         bzmn_e(lb:ub) = bzmn_e(lb:ub) - (bsqr(lb:ub)*r1(lb:ub,1)
     1              + guus(lb:ub)*zu(lb:ub,1) + guu(lb:ub)*zu(lb:ub,0))
         lv_e(lb:ub) = lv_e(lb:ub)*shalf(1:ns)
         lu_o(lb:ub) = dshalfds*lu_e(lb:ub)

         armn_o(lb:ub-1) = armn_o(lb+1:ub) - armn_o(lb:ub-1)
     1             - zu(lb:ub-1,0)*bsqr(lb:ub-1)
     1             + p5*(lv_e(lb:ub-1) + lv_e(lb+1:ub))
         armn_o(ub) = -armn_o(ub)-zu(ub,0)*bsqr(ub)+p5*lv_e(ub)
         azmn_o(lb:ub-1) = azmn_o(lb+1:ub) - azmn_o(lb:ub-1) 
     1             + ru(lb:ub-1,0)*bsqr(lb:ub-1)
         azmn_o(ub) = -azmn_o(ub)+ru(ub,0)*bsqr(ub)
         brmn_o(lb:ub-1) = p5*(brmn_o(lb:ub-1) + brmn_o(lb+1:ub))
         brmn_o(ub) = p5*brmn_o(ub)
         bzmn_o(lb:ub-1) = p5*(bzmn_o(lb:ub-1) + bzmn_o(lb+1:ub))
         bzmn_o(ub) = p5*bzmn_o(ub)
         lu_o(lb:ub-1)   = lu_o(lb:ub-1) + lu_o(lb+1:ub)

         guu(lb:ub)  = guu(lb:ub) * sqrts(1:ns)**2
         bsqr(lb:ub) = gvv(lb:ub) * sqrts(1:ns)**2

         armn_o(lb:ub) = armn_o(lb:ub) - (zu(lb:ub,1)*lu_o(lb:ub)
     1              + bsqr(lb:ub)*r1(lb:ub,1) + gvvs(lb:ub)*r1(lb:ub,0))
         azmn_o(lb:ub) = azmn_o(lb:ub) + ru(lb:ub,1)*lu_o(lb:ub)
         brmn_o(lb:ub) = brmn_o(lb:ub) + z1(lb:ub,1)*lu_o(lb:ub)
     1              -(guu(lb:ub)*ru(lb:ub,1) + guus(lb:ub)*ru(lb:ub,0))
         bzmn_o(lb:ub) = bzmn_o(lb:ub) - (r1(lb:ub,1)*lu_o(lb:ub)
     1              + guu(lb:ub)*zu(lb:ub,1) + guus(lb:ub)*zu(lb:ub,0))

         IF (lthreed) THEN

         guvs(lb:ub) = guv(lb:ub)*shalf(1:ns)
         guv(lb:ub-1)  = p5*(guv(lb:ub-1) + guv(lb+1:ub))
         guv(ub) = p5*guv(ub)
         guvs(lb:ub-1) = p5*(guvs(lb:ub-1) + guvs(lb+1:ub))
         guvs(ub) = p5*guvs(ub)

         brmn_e(lb:ub) = brmn_e(lb:ub) 
     1          - (guv(lb:ub)*rv(lb:ub,0) + guvs(lb:ub)*rv(lb:ub,1))
         bzmn_e(lb:ub) = bzmn_e(lb:ub) 
     1          - (guv(lb:ub)*zv(lb:ub,0) + guvs(lb:ub)*zv(lb:ub,1))
         crmn_e(lb:ub) = guv(lb:ub) *ru(lb:ub,0) 
     1                 + gvv(lb:ub) *rv(lb:ub,0)
     2          + gvvs(lb:ub)*rv(lb:ub,1) + guvs(lb:ub)*ru(lb:ub,1)
         czmn_e(lb:ub) = guv(lb:ub) *zu(lb:ub,0) 
     1                 + gvv(lb:ub) *zv(lb:ub,0)
     2          + gvvs(lb:ub)*zv(lb:ub,1) + guvs(lb:ub)*zu(lb:ub,1)
         guv(lb:ub) = guv(lb:ub) *sqrts(lb:ub)*sqrts(lb:ub)
         brmn_o(lb:ub) = brmn_o(lb:ub) 
     1          - (guvs(lb:ub)*rv(lb:ub,0) + guv(lb:ub)*rv(lb:ub,1))
         bzmn_o(lb:ub) = bzmn_o(lb:ub) 
     1          - (guvs(lb:ub)*zv(lb:ub,0) + guv(lb:ub)*zv(lb:ub,1))
         crmn_o(lb:ub) = guvs(lb:ub)*ru(lb:ub,0) 
     1                 + gvvs(lb:ub)*rv(lb:ub,0)
     2          + bsqr(lb:ub)*rv(lb:ub,1) + guv(lb:ub) *ru(lb:ub,1)
         czmn_o(lb:ub) = guvs(lb:ub)*zu(lb:ub,0) 
     1                 + gvvs(lb:ub)*zv(lb:ub,0)
     2          + bsqr(lb:ub)*zv(lb:ub,1) + guv(lb:ub) *zu(lb:ub,1)
         END IF

!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
         IF (ivac .GE. 1) THEN
         armn_e(ub) = armn_e(ub) + zu0(ub)*rbsq(l)
         armn_o(ub) = armn_o(ub) + zu0(ub)*rbsq(l)
         azmn_e(ub) = azmn_e(ub) - ru0(ub)*rbsq(l)
         azmn_o(ub) = azmn_o(ub) - ru0(ub)*rbsq(l)
         ENDIF

!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
#ifndef _HBANGLE
         rcon(lb:ub,0) = (rcon(lb:ub,0) - rcon0(lb:ub)) * gcon(lb:ub)
         zcon(lb:ub,0) = (zcon(lb:ub,0) - zcon0(lb:ub)) * gcon(lb:ub)
         brmn_e(lb:ub) = brmn_e(lb:ub) + rcon(lb:ub,0)
         bzmn_e(lb:ub) = bzmn_e(lb:ub) + zcon(lb:ub,0)
         brmn_o(lb:ub) = brmn_o(lb:ub)+ rcon(lb:ub,0)*sqrts(lb:ub)
         bzmn_o(lb:ub) = bzmn_o(lb:ub)+ zcon(lb:ub,0)*sqrts(lb:ub)
         rcon(lb:ub,0) =  ru0(lb:ub) * gcon(lb:ub)
         zcon(lb:ub,0) =  zu0(lb:ub) * gcon(lb:ub)
         rcon(lb:ub,1) = rcon(lb:ub,0) * sqrts(1:ns)
         zcon(lb:ub,1) = zcon(lb:ub,0) * sqrts(1:ns)
#endif

      END DO
!$OMP END PARALLEL DO      



 100  CONTINUE


      END SUBROUTINE forces
