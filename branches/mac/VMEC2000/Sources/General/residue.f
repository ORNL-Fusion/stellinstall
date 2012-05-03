      SUBROUTINE residue (gcr, gcz, gcl, xcz)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc, zsc, 
     1                       meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE vsvd
      USE xstuff
      USE precon2d, ONLY: ictrl_prec2d, block_precond
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), 
     1  INTENT(inout) :: gcr, gcz, gcl
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), 
     1  INTENT(in) :: xcz
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: n0 = 0, m0 = 0, m1 = 1
      INTEGER, PARAMETER :: n3d = 0, nasym = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsfix, iflag, jedge, delIter
      REAL(rprec) :: r1, tnorm, bz1
C-----------------------------------------------

!
!       IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!       INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!       ( ZCS = RSS, ZCC = RSC ARE THE CORRECT POLAR RELATIONS; HERE WE USE
!         THE SIMPLER RELATION ZCS = 0 INSTEAD)
!
!
!     IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!     INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        RSS(n) = ZCS(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        RSC(n) = ZCC(n), all n 
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0  
!        3D, ASYM: GC(zcc) = 0
!
      IF (lthreed .and. fsqz.lt.1.0_dp)
     1   CALL constrain_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs),
     2                     xcz(:,:,1,zcs), 1)
      IF (lasym)
     1   CALL constrain_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc),
     2                     xcz(:,:,1,zcc), 2)
!
!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
!
      IF (ictrl_prec2d .ge. 2) RETURN
!
!     PUT FORCES INTO PHIFSAVE UNITS THAT PRECONDITIONERS, FNORM ARE IN
!
      IF (phifac .eq. zero) THEN
         STOP 'phifac = 0 in residue'
      ELSE
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      END IF


      IF (lrecon) THEN
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS

        r1 = SUM(gcr(:ns,n0,m0,1))
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        IF (iopt_raxis.gt.0 .and. iresidue.eq.2
     1     .and. fsq.lt.fopt_axis) iresidue = 3
        IF (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      ELSE
!
!     ADJUST PHIEDGE
!
         IF (imovephi .gt. 0) CALL movephi1 (gphifac)
      ENDIF
      gc(neqs1) = gphifac

!
!     CONSTRUCT INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0    
!SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
!ADD A V3FIT RELATED FLAG? ADD fsq criterion first
      delIter = iter2-iter1
!      IF (delIter.lt.50 .and. 
!     1   (fsqr+fsqz).lt.1.E-6_dp) jedge = 1

      IF (l_v3fit) THEN
!  Coding for when run by V3FIT. Needed for correct computation
!  of partial derivatives
         IF (iter2-iter1.lt.50) jedge = 1
      ELSE
!  Coding for VMEC2000 run stand-alone
         IF (delIter.lt.50 .and.
     1      (fsqr+fsqz).lt.1.E-6_dp) jedge = 1
      ENDIF

      CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

      bz1 = 2*hs*r1/((rbtor*tnorm)**2)
      fsql = bz1*SUM(gcl*gcl)
      fedge = r1*fnorm*SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      IF (ictrl_prec2d .eq. 1) THEN

         CALL block_precond(gc)

         IF (.not.lfreeb .and. ANY(gcr(ns,:,:,:) .ne. zero))
     1      STOP 'gcr(ns) != 0 for fixed boundary in residue'
         IF (.not.lfreeb .and. ANY(gcz(ns,:,:,:) .ne. zero))
     1      STOP 'gcz(ns) != 0 for fixed boundary in residue'
         IF (ANY(gcl(:,1:,0,zsc) .ne. zero))
     1      STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         IF (lthreed) THEN
            IF (ANY(gcl(:,n0,:,zcs) .ne. zero))
     1      STOP 'gcl(n=0,m,cs) != 0 in residue'
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE

         IF (lthreed) 
     1      CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
         IF (lasym)
     1      CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

         iflag = 0
         CALL scalfor (gcr, arm, brm, ard, brd, crd, iflag)
         iflag = 1
         CALL scalfor (gcz, azm, bzm, azd, bzd, crd, iflag)

         CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)

         gcl = faclam*gcl
         r1  = r1*SUM(phip(2:ns)**2)
         fsql1 = SUM(gcl*gcl)/r1

!
!     SOFT-STARTUP AVOIDS SOME "HANG-UPS"  (04/2002)
!
!         fac = one / (one+fsqr+fsqz+fsql)
!         gcr = fac*gcr
!         gcz = fac*gcz

      ENDIF

      END SUBROUTINE residue

      SUBROUTINE constrain_m1(gcr, gcz, z1, nType)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: ntmax
      USE precon2d, ONLY: ictrl_prec2d, z1_save
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nType
      REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
      REAL(rprec), DIMENSION(ns,0:ntor), INTENT(in)    :: z1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL, PARAMETER :: FThreshold = 1.E-6_dp
      REAL(rprec) :: temp(ns,0:ntor)
C-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      IF (lconm1) THEN
         temp = gcr
         gcr = p5*(gcr + gcz)
         gcz = p5*(temp- gcz)
      END IF

      IF (fsqz.lt.FThreshold .or. ictrl_prec2d.ne.0) gcz = 0

      RETURN

!OLD STYLE: NO LONGER NECESSARY, WE NOW ADD 1 TO DIAGONAL HESSIAN
      IF (ictrl_prec2d .gt. 0) THEN
         IF (ictrl_prec2d .eq. 2) z1_save(:,:,nType) = z1
!MUST DO IT THIS WAY TO AVOID SINGULAR HESSIAN
         gcz = z1 - z1_save(:,:,nType)         
!      ELSE IF (fsqz.lt.FThreshold .and. .not.lconm1) THEN
      ELSE IF (fsqz.lt.FThreshold) THEN
         gcz = 0
      END IF
 
      END SUBROUTINE constrain_m1

      SUBROUTINE scale_m1(gcr, gcz)
      USE vmec_main, p5 => cp5
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER :: n
      REAL(rprec) :: fac(ns)
C-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac = (ard(:,2)+brd(:,2))/
     1        (p5*(ard(:,2)+brd(:,2)+azd(:,2)+bzd(:,2)))
      DO n = 0, ntor
         gcr(:,n) = fac*gcr(:,n)
      END DO

      fac = (azd(:,2)+bzd(:,2))/
     1       (p5*(ard(:,2)+brd(:,2)+azd(:,2)+bzd(:,2)))
      DO n = 0, ntor
         gcz(:,n) = fac*gcz(:,n)
      END DO
  
      END SUBROUTINE scale_m1
