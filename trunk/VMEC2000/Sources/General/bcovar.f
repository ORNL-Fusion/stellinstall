!  Small equilibria resolutions can be slower when run in parallel. Allow V3FIT
!  to be compiled with OpenMP support while disabling VMEC OpenMP support.
#ifdef NO_VMEC_OMP
#define OMP NO_OMP
#endif

      SUBROUTINE bcovar (lu, lv, lmnsc00)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale
      USE realspace
      USE vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1             rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2             bsubu_e => clmn_e, bsubv_e => blmn_e, 
     3             bsubu_o => clmn_o, bsubv_o => blmn_o,
     4             bsq => bzmn_o, phipog => brmn_o
      USE vsvd, ONLY: phifac, phifsave, imovephi
      USE xstuff, ONLY: xc
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact,
     1                    ctor_prec2d
#ifdef _HBANGLE
      USE angle_constraints, ONLY: precondn_rho, ard2, arm2,            
     1                             azd2, azm2, brd2, brm2, bzd2, bzm2
#endif
      USE fbal
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt,0:1), INTENT(inout) :: lu, lv
      REAL(rprec), DIMENSION(ns), INTENT(inout) :: lmnsc00
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(rprec), PARAMETER :: c1p5 = (one + p5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, js, ndim, ub, lb
      REAL(rprec) :: r2, volume, curpol_temp
#ifndef _HBANGLE
      REAL(rprec), DIMENSION(2:ns-1):: maznorm, marnorm
      REAL(rprec) :: tcon_mul
#endif
      REAL(rprec), POINTER, DIMENSION(:) :: luu, luv, lvv, tau
      REAL(rprec), DIMENSION(:), POINTER :: bsupu, bsubuh, 
     1                                      bsupv, bsubvh, r12sq
      REAL(rprec), DIMENSION(ns):: temp1, temp2, temp3
      REAL(rprec) :: t1, t2
      LOGICAL :: lctor
!-----------------------------------------------
      ndim = 1+nrzt

!
!     POINTER ALIAS ASSIGNMENTS
!   
      tau => extra1(:,1);  luu => extra2(:,1);  
      luv => extra3(:,1);  lvv => extra4(:,1)

      bsupu => luu;  bsubuh => bsubu_o
      bsupv => luv;  bsubvh => bsubv_o
      r12sq => bsq

      vp=0
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub,temp1,temp2,temp3)
!$OMP& REDUCTION(+:vp)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l

         temp3 = sqrts(lb:ub)*sqrts(lb:ub)

         temp1= ru(lb:ub,0)*ru(lb:ub,0) 
     1        + zu(lb:ub,0)*zu(lb:ub,0) + temp3*
     2        ( ru(lb:ub,1)*ru(lb:ub,1)
     3        + zu(lb:ub,1)*zu(lb:ub,1))
         
         temp2= (ru(lb:ub,0)*ru(lb:ub,1) 
     1        +  zu(lb:ub,0)*zu(lb:ub,1))*2
         guu(lb)=0
         guu(lb+1:ub) = p5*(temp1(2:ns) + temp1(1:ns-1) + 
     1        shalf(lb+1:ub)*(temp2(2:ns) + temp2(1:ns-1)))

         IF (lthreed) THEN
         temp1= ru(lb:ub,0)*rv(lb:ub,0)
     1        + zu(lb:ub,0)*zv(lb:ub,0) + temp3*
     2        ( ru(lb:ub,1)*rv(lb:ub,1) 
     3        + zu(lb:ub,1)*zv(lb:ub,1) )
         temp2= ru(lb:ub,0)*rv(lb:ub,1) 
     1        + ru(lb:ub,1)*rv(lb:ub,0)
     2        + zu(lb:ub,0)*zv(lb:ub,1) 
     3        + zu(lb:ub,1)*zv(lb:ub,0)

         guv(lb)=0
         guv(lb+1:ub) = p5*(temp1(2:ns) + temp1(1:ns-1) + 
     1        shalf(lb+1:ub)*(temp2(2:ns) + temp2(1:ns-1)))
         ELSE
            guv(lb:ub) = 0
         END IF

         temp1= r1(lb:ub,0)*r1(lb:ub,0) 
     1        + temp3*r1(lb:ub,1)*r1(lb:ub,1)
         temp2 = 2*r1(lb:ub,0)*r1(lb:ub,1)
         IF (lthreed) THEN
         temp1= temp1
     1        + rv(lb:ub,0)*rv(lb:ub,0)
     3        + zv(lb:ub,0)*zv(lb:ub,0) 
     4        + temp3*( rv(lb:ub,1)*rv(lb:ub,1)
     5        + zv(lb:ub,1)*zv(lb:ub,1) )
          temp2 = temp2
     1        + (rv(lb:ub,0)*rv(lb:ub,1) 
     2        +  zv(lb:ub,0)*zv(lb:ub,1))*2
         END IF
         gvv(lb)=0
         gvv(lb+1:ub) = p5*(temp1(2:ns) + temp1(1:ns-1) + 
     1      shalf(lb+1:ub)*(temp2(2:ns) + temp2(1:ns-1)))
    
         tau(lb:ub) = gsqrt(lb:ub)
         gsqrt(lb:ub) = r12(lb:ub)*gsqrt(lb:ub)
         gsqrt(lb)=gsqrt(lb+1)

         phipog(lb)=0
         phipog(lb+1:ub)=one/gsqrt(lb+1:ub)
         
         vp(2:ns)=vp(2:ns)+signgs*gsqrt(lb+1:ub)*wint(lb+1:ub)
!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!     TO ACCOMODATE LRFP=T CASES, THE OVERALL PHIP FACTOR (PRIOR TO v8.46)
!     HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
!
!     NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE
!     WILL ADD CHIP IN CALL TO ADD_FLUXES. THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
!
!          BSUPU = PHIPOG*(chip - LAMV*LAMSCALE),
!          BSUPV = PHIPOG*(phip + LAMU*LAMSCALE)
!
         lu(lb:ub,:) = lu(lb:ub,:)*lamscale
         lv(lb:ub,:) = lv(lb:ub,:)*lamscale
         lu(lb:ub,0)=lu(lb:ub,0)+phipf(1:ns)

         bsupu(lb+1:ub)=p5*phipog(lb+1:ub)*(lv(lb+1:ub,0)+lv(lb:ub-1,0) 
     1              + shalf(2:ns)*(lv(lb+1:ub,1) + lv(lb:ub-1,1)))
         bsupv(lb+1:ub)=p5*phipog(lb+1:ub)*(lu(lb+1:ub,0)+lu(lb:ub-1,0) 
     1              + shalf(2:ns)*(lu(lb+1:ub,1) + lu(lb:ub-1,1)))
!v8.49: add ndim points
         bsupu(lb)=0
         bsupv(lb)=0
      END DO 
!$OMP END PARALLEL DO 

      IF (iter2 .eq. 1) voli = twopi*twopi*hs*SUM(vp(2:ns))

!
!     UPDATE IOTA EITHER OF TWO WAYS:
!     1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv 
!         FOR iotas  
!     2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv. 
!
!     NEED TO DO IT WAY (#2) TO EASILY COMPUTE THE HESSIAN ELEMENTS DUE TO LAMBDA-VARIATIONS. 
!     IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY [USE THIS COMPONENT SO IT 
!     WILL WORK EVEN FOR 2D PLASMA], ALTHOUGH ITS VARIATION IS LIKE THAT OF LV-CS(0,0), 
!     WITH N -> 1 IN THE HESSIAN CALCULATION ROUTINES (Compute_Hessian_Flam_lam, etc.)
!

      IF (ncurr .eq. 1) THEN
         IF (ictrl_prec2d .eq. 2) THEN
            lmnsc00(2:ns) = chips(2:ns)            !Initialize
!!          lmnsc00(2:ns) = iotas(2:ns)            !Initialize
         ELSE IF (ictrl_prec2d .ne. 0) THEN
            chips(2:ns) = lmnsc00(2:ns)            !Evolution
!!          iotas(2:ns) = lmnsc00(2:ns)            !Evolution
         END IF
      END IF

!     COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
      CALL add_fluxes(phipog, bsupu, bsupv, ictrl_prec2d.eq.0)

      wb = 0
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
!$OMP& REDUCTION(+:wb)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l         
!
!     COMPUTE LAMBDA FORCE KERNELS (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
         bsubuh(lb:ub)=guu(lb:ub)*bsupu(lb:ub)+guv(lb:ub)*bsupv(lb:ub)
         bsubvh(lb:ub)=guv(lb:ub)*bsupu(lb:ub)+gvv(lb:ub)*bsupv(lb:ub)

!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
         bsq(lb:ub) = p5*(bsupu(lb:ub)*bsubuh(lb:ub) 
     1              +     bsupv(lb:ub)*bsubvh(lb:ub))

         wb = wb+hs*ABS(SUM(wint(lb:ub)*gsqrt(lb:ub)*bsq(lb:ub)))
      END DO 
!$OMP END PARALLEL DO 

#ifdef _ANIMEC
!SPH: MAKE CALL HERE (bsubX_e are used for scratch arrays)
      CALL an_pressure(bsubu_e, bsubv_e)

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l         
         bsq(lb+1:ub) = bsq(lb+1:ub) + pperp(lb+1:ub)

!WAC-SPH: MODIFY EFFECTIVE CURRENT K = curl(sigma_an*B)
         phipog(lb:ub) = phipog(lb:ub)*sigma_an(lb:ub)
         bsubuh(lb:ub) = bsubuh(lb:ub)*sigma_an(lb:ub)
         bsubvh(lb:ub) = bsubvh(lb:ub)*sigma_an(lb:ub)
      END DO 
!$OMP END PARALLEL DO 

#else
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*SUM(vp(2:ns)*pres(2:ns))

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l         
         bsq(lb+1:ub) = bsq(lb+1:ub) + pres(2:ns)
      END DO 
!$OMP END PARALLEL DO 
#endif

!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
!SPH122407-MOVED HERE: COMPUTE LAMBDA FULL MESH FORCES
!     NOTE: bsubu_e is used here ONLY as a temporary array
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l         
         lvv(lb:ub) = phipog(lb:ub)*gvv(lb:ub)
         bsubv_e(lb:ub-1)=p5*(lvv(lb:ub-1)+lvv(lb+1:ub))*lu(lb:ub-1,0)
         bsubv_e(ub) = p5*lvv(ub)*lu(ub,0)

         lvv(lb:ub) = lvv(lb:ub)*shalf(1:ns)
         bsubu_e(lb:ub)=guv(lb:ub)*bsupu(lb:ub)*sigma_an(lb:ub)          !sigma_an=1 isotropic
         bsubv_e(lb:ub-1) = bsubv_e(lb:ub-1) 
     1            + p5*((lvv(lb:ub-1) + lvv(lb+1:ub))*lu(lb:ub-1,1)
     2            +      bsubu_e(lb:ub-1) + bsubu_e(lb+1:ub))
         bsubv_e(ub) = bsubv_e(ub) + p5*(lvv(ub)*lu(ub,1)+
     1                 bsubu_e(ub))
      END DO 
!$OMP END PARALLEL DO 

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
!WAC: UPDATE buco, bvco AFTER pressure called
#ifdef _ANIMEC
      IF (iequi .EQ. 1) papr = pmap*pres/vp
#endif
      CALL calc_fbal(bsubuh, bsubvh)
    
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
!
!     (SPH:08/19/04)
!     MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
!     CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
!     THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
!
!     ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
!     TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
!
!     TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
!     FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)
!
!

!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY

      IF (lHess_exact) THEN
         lctor = lfreeb .and. ictrl_prec2d.ne.0      !Yields correct hessian near edge
      ELSE
         lctor = lfreeb .and. ictrl_prec2d.gt.1      !Yields better accuracy in solution
      END IF
      IF (lctor) THEN       
         IF (ictrl_prec2d .eq. 2) 
     1       ctor_prec2d = signgs*twopi*p5*(buco(ns) - buco(ns1))
         ctor = ctor_prec2d + signgs*twopi*buco(ns)
      ELSE
         ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))
      END IF

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l         
         IF (bsubuh(lb) .ne. zero) STOP 'BSUBUH != 0 AT JS=1'
         IF (bsubvh(lb) .ne. zero) STOP 'BSUBVH != 0 AT JS=1'
         bsubu_e(lb:ub-1) = p5*(bsubuh(lb:ub-1) + bsubuh(lb+1:ub))
         bsubu_e(ub) = p5*bsubuh(ub)
         bsubv_e(lb:ub-1) = bsubv_e(lb:ub-1)*bdamp(1:ns-1) + 
     1     p5*(1-bdamp(1:ns-1))*(bsubvh(lb:ub-1) + bsubvh(lb+1:ub))
         bsubv_e(ub) = bsubv_e(ub)*bdamp(ns) + p5*(1-bdamp(ns))*
     1                 bsubvh(ub)
!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS: NOTE THAT lu=>czmn, lv=>crmn externally
!     SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
!
         IF (iequi .EQ. 1) THEN
            lu(lb:ub,0) = bsupv(lb:ub)
            lv(lb:ub,0) = bsupu(lb:ub)
         END IF
      END DO 
!$OMP END PARALLEL DO 
!
!     COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
!     NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
!
      IF ((MOD(iter2-iter1,ns4).eq.0 .and. iequi.eq.0) 
     1        .and. ictrl_prec2d.eq.0) THEN
         phifsave = phifac
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
         DO l=1,nznt
            lb=ns*(l-1)+1
            ub=ns*l         
            phipog(lb:ub) = phipog(lb:ub)*wint(lb:ub)
         END DO
!$OMP END PARALLEL DO 
         CALL lamcal(phipog, guu, guv, gvv)
         CALL precondn(bsupv,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,
#ifdef _HBANGLE
     2                 arm2, ard2, brm2, brd2,
#endif
     3                 crd,rzu_fac,cos01)
         CALL precondn(bsupv,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                 r1(1,1),azm,azd,bzm,bzd,
#ifdef _HBANGLE
     2                 azm2, azd2, bzm2, bzd2,
#endif
     3                 crd,rru_fac,sin01)

         rzu_fac(2:ns-1) = sqrts(2:ns-1)*rzu_fac(2:ns-1)
         rru_fac(2:ns-1) = sqrts(2:ns-1)*rru_fac(2:ns-1)
         frcc_fac(2:ns-1) = one/rzu_fac(2:ns-1);  rzu_fac = rzu_fac/2
         fzsc_fac(2:ns-1) =-one/rru_fac(2:ns-1);  rru_fac = rru_fac/2
#ifdef _HBANGLE
         CALL precondn_rho
#endif
         volume = hs*SUM(vp(2:ns))
         r2 = MAX(wb,wp)/volume

         t1=0; t2=0
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
!$OMP& REDUCTION(+:t1,t2)
         DO l=1,nznt
            lb=ns*(l-1)+2
            ub=ns*l          
            guu(lb:ub) = guu(lb:ub)*r12(lb:ub)**2
            t1 = t1 + SUM(guu(lb:ub)*wint(lb:ub))
            t2 = t2 + SUM((bsubuh(lb:ub)**2 + bsubvh(lb:ub)**2)
     1                    *wint(lb:ub))
         END DO
!$OMP END PARALLEL DO 

         fnorm = one/(t1*r2*r2)                                    !Norm, unpreconditioned R,Z forces
         fnorm1 = one/SUM(xc(1+ns:2*irzloff)**2)                   !Norm for preconditioned R,Z forces
         fnormL = one/(t2*lamscale**2)                             !Norm for unpreconditioned Lambda force

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
#ifndef _HBANGLE
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 from old-style files
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((4*r0scale**2)**2)                   !!Scaling of ard, azd (2*r0scale**2); 
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)

         marnorm=0
         maznorm=0
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
!$OMP& REDUCTION(+:marnorm,maznorm)
         DO l=1,nznt
            lb=ns*(l-1)+2
            ub=ns*l-1
            marnorm = marnorm + wint(lb:ub)*ru0(lb:ub)**2
            maznorm = maznorm + wint(lb:ub)*zu0(lb:ub)**2
         END DO
!$OMP END PARALLEL DO 
         IF (ANY(marnorm(2:ns-1).EQ.zero .OR. maznorm(2:ns-1).EQ.zero))
     1        STOP 'arnorm or aznorm=0 in bcovar'

         DO js = 2, ns-1
           tcon(js) = MIN(ABS(ard(js,1)/marnorm(js)),
     1                ABS(azd(js,1)/maznorm(js))) * tcon_mul*(32*hs)**2
         END DO

         tcon(ns) = p5*tcon(ns-1)
         IF (lasym) tcon = p5*tcon
#endif
      ENDIF

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      IF (iequi .eq. 1) THEN

!         IF (.FALSE.) THEN
         DO js = ns-1,2,-1
            DO l = js, nrzt, ns
               bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
            END DO
         END DO

!     ADJUST <bsubvh> AFTER MESH-BLENDING
         DO js = 2, ns
            curpol_temp = fpsi(js) 
     1                  - SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))
           DO l = js, nrzt, ns
               bsubvh(l) = bsubvh(l) + curpol_temp
            END DO
         END DO
!         END IF

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)

         RETURN

      END IF
!$OMP PARALLEL DO 
!$OMP& PRIVATE(l,lb,ub)
      DO l=1,nznt
         lb=ns*(l-1)+1
         ub=ns*l      

!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
         bsubu_e(lb:ub) = -lamscale*bsubu_e(lb:ub)
         bsubv_e(lb:ub) = -lamscale*bsubv_e(lb:ub)
         bsubu_o(lb:ub)  = sqrts(lb:ub)*bsubu_e(lb:ub)
         bsubv_o(lb:ub)  = sqrts(lb:ub)*bsubv_e(lb:ub)
!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
!WAC, SPH122407: sigma_an (=1 for isotropic case)
         lvv(lb:ub) = gsqrt(lb:ub)*sigma_an(lb:ub)
         guu(lb:ub) = bsupu(lb:ub)*bsupu(lb:ub)*lvv(lb:ub)
         guv(lb:ub) = bsupu(lb:ub)*bsupv(lb:ub)*lvv(lb:ub)
         gvv(lb:ub) = bsupv(lb:ub)*bsupv(lb:ub)*lvv(lb:ub)
         lv(lb:ub,0)= bsq(lb:ub)*tau(lb:ub)
         lu(lb:ub,0)= bsq(lb:ub)*r12(lb:ub)
      END DO
!$OMP END PARALLEL DO 

      END SUBROUTINE bcovar

#ifdef _ANIMEC
      SUBROUTINE an_pressure(gp, scratch1)
      USE stel_kinds, ONLY: rprec, dp
      USE realspace, ONLY: sigma_an, wint, pperp, ppar, onembc
      USE vforces, r12 => armn_o, gsqrt => azmn_o,
     &             bsq => bzmn_o
      USE vmec_main, ONLY: phot, tpotb, pppr, pmap, mass, pres, vp,
     &                     wp, gamma, wpar, wper, ns, nznt, nrzt, 
     &                     zero, one, nthreed, pperp_ns1=>dbsq,
     &                     bcrit, medge, phedg, hs, pperp_ns
      USE vmec_params, ONLY: signgs
      USE fbal
!
!     WAC (11/30/07): See description of anisotropic pressure below
!     SPH (12/24/07): Replaced "gp" with bsubu_e to avoid overwriting phipog
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt)               :: gp, scratch1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER     :: js   , lk   , l
      REAL        :: pparden, pres_pv, pppr_pv
!-----------------------------------------------
!
!********0*********0*********0*********0*********0*********0*********0**
!                                                                      *
!                 Anisotropic Pressure Model                           *
!                  specific to case where:                             *
!          a Bi-Maxwellian distribution is considered (by J. Graves)   *
!          p_parallel(s,B) = pth(1 + phot(s)*H(s,B))                     *
!      H(s,B)=(B/B_crit)/[1-(T_perp/T_par)(1-B/B_crit)] for B>B_crit   *
!      For B<B_crit,                                                   *
!      H(s,B)=H(s,B>B_crit){1-2[(T_perp/T_par)(1-B/B_crit)]^(5/2) /    *
!                                    [1+(T_perp/T_par)(1-B/B_crit)]}   *
!                                                                      *
!********0*********0*********0*********0*********0*********0*********0**
!
!********0*********0*********0*********0*********0*********0*********0**
!   1.  Compute Thermal Pressure and Hot Particle Pressure Amplitude.  *
!********0*********0*********0*********0*********0*********0*********0**
      DO js = 1,ns
         gp(js:nrzt:ns) = tpotb(js)
         scratch1(js:nrzt:ns) = phot(js)
      END DO

      bsq(1:nrzt:ns) = 0
      onembc = one - SQRT(2*bsq(1:nrzt))/bcrit
      WHERE (one .ne. gp*onembc) sigma_an = (one-onembc)/(one-gp*onembc)

      WHERE (onembc .gt. zero) sigma_an = sigma_an*
     1      (one-2*(gp*onembc)**2.5_dp/(one+gp*onembc))

      pperp = (1 + scratch1*sigma_an)*gsqrt(1:nrzt)

      DO js = 2,ns
         pmap(js) = DOT_PRODUCT(pperp(js:nrzt:ns), wint(js:nrzt:ns))
      END DO
        
      DO js = 2,ns
           pmap(js) = signgs*pmap(js)
           pres(js) = mass(js) / pmap(js)**gamma
           pppr(js) = pres(js) * phot(js)
      END DO
!
!********0*********0*********0*********0*********0*********0*********0**
!   3.  Compute P-Parallel, P-Perp.                                    *
!********0*********0*********0*********0*********0*********0*********0**

      DO js = 2,ns
         scratch1(js:nrzt:ns) = pppr(js)
      END DO

      ppar = scratch1*sigma_an

!FORTRAN 95 CONSTRUCT ALLOWS ELSEWHERE (TEST), F90 DOES NOT
#if defined(WIN32)
      WHERE (onembc .le. zero)
         pperp = gp*(one-onembc)/(one-gp*onembc)*ppar
      ELSEWHERE
         WHERE (onembc*gp .eq. one)
            ppar = 2*scratch1*(one - onembc)
            pperp= (7*ppar*(gp-one))/16
         ELSEWHERE
            pperp = (one - (5-(gp*onembc)**2)*
     &              (gp*onembc)**1.5_dp/(one+gp*onembc)**2)*scratch1
     &              *gp*(one-onembc)**2/(one-gp*onembc)**2
         ENDWHERE
      ENDWHERE

#else
      WHERE (onembc .le. zero)
 
         pperp = gp*(one-onembc)/(one-gp*onembc)*ppar

      ELSEWHERE (onembc*gp .eq. one)

         ppar = 2*scratch1*(one - onembc)
         pperp= (7*ppar*(gp-one))/16

      ELSEWHERE

         pperp = (one - (5-(gp*onembc)**2)*
     &           (gp*onembc)**1.5_dp/(one+gp*onembc)**2)*scratch1
     &           *gp*(one-onembc)**2/(one-gp*onembc)**2

      END WHERE
#endif

      DO js = 2,ns
         scratch1(js:nrzt:ns) = pres(js)
      END DO

      ppar = ppar + scratch1
      pperp= pperp+ scratch1

!
!********0*********0*********0*********0*********0*********0*********0**
!   4.  Compute P_perp at the plasma-vacuum interface.                 *
!********0*********0*********0*********0*********0*********0*********0**
!
      pparden = MAX(pppr(ns-1),1.e-30_dp)
      DO lk=1,nznt
         l = ns-1 + ns*(lk-1)
         pperp_ns1(lk) = (pperp(l)-pres(ns-1))/pparden
      END DO
      pparden = MAX(pppr(ns),1.e-30_dp)
      DO lk=1,nznt
         l = ns + ns*(lk-1)                             !!SPH12-27-12: l = ns, not ns-1
         pperp_ns(lk) = (pperp(l)-pres(ns))/pparden
      END DO
            
      pres_pv = medge / (1.5_dp*pmap(ns)-0.5_dp*pmap(ns-1))**gamma
      pppr_pv = pres_pv * phedg

      DO lk=1,nznt
         pperp_ns(lk)=(1.5_dp*pperp_ns(lk)-
     &                 0.5_dp*pperp_ns1(lk))*pppr_pv + pres_pv
      END DO
!
!********0*********0*********0*********0*********0*********0*********0**
!   5.  Compute Sigma_an. Determine Volume Averaged Pressures.            *
!********0*********0*********0*********0*********0*********0*********0**
!
      wper = 0
      gp = gsqrt(1:nrzt)*pperp
      sigma_an = one + (pperp-ppar)/(2*bsq(1:nrzt))

      pperp(1:nrzt:ns) = 0
      sigma_an(1:nrzt:ns) = 1

      IF (ALL(phot.eq.zero) .AND. ANY(sigma_an.ne.one)) 
     1   STOP 'SIGMA_AN != 1'

      pppr(1) = 0
      DO js = 2,ns
         pppr(js) = DOT_PRODUCT(gp(js:nrzt:ns),wint(js:nrzt:ns))
      END DO
      pppr(2:ns) = signgs*pppr(2:ns)/vp(2:ns)

      wp    = hs*DOT_PRODUCT(vp(2:ns),pres(2:ns))
      wpar  = hs*DOT_PRODUCT(pmap(2:ns),pres(2:ns))
!      whpar = wpar - wp
      wper  = hs*DOT_PRODUCT(vp(2:ns),pppr(2:ns))
!      whper = wper - wp


      END SUBROUTINE an_pressure
#endif
