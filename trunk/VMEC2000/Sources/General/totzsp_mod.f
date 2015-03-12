#define oldtotzsps
      MODULE totzsp_mod
      USE stel_kinds, ONLY: rprec
      USE vmec_main
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, n0=0
!-----------------------------------------------
#if !defined(oldtotzsps)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE:: mycosnv, mysinnv,
     1     mysinnvn,mycosnvn
      LOGICAL:: first=.true.
      INTEGER:: d3,d4
      REAL(rprec):: eins,nix
#endif

      CONTAINS

#if defined(oldtotzsps)
      SUBROUTINE totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, lb, ub
      INTEGER :: ioff, joff, mj, ni, nsz
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
!-----------------------------------------------
!
!     WORK1 Array of inverse transforms in toroidal angle (zeta), for all radial positions
!     NOTE: ORDERING OF LAST INDEX IS DIFFERENT HERE THAN IN PREVIOUS VMEC2000 VERSIONS
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(stored at rss) = .5(rss + zcs),
!     R-(stored at zcs) = .5(rss - zcs), TO EXTERNAL ("PHYSICAL") rss, zcs FORMS. NEED THIS EVEN 
!     WHEN COMPUTING HESSIAN FOR FREE BOUNDARY (rmnss, zmncs at JS=NS needed in vacuum call)
!
!     WHEN COMPUTING PRECONDITIONER, USE FASTER HESSIAN VERSION (totzsps_hess) INSTEAD.
      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)            !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)      !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)    !!COS(mu) SIN(nv)
#ifndef _HBANGLE
         CALL convert_sym (rmnss, zmncs)
#endif
      END IF

!v8.50: Norm for preconditioned R,Z forces: scale to boundary value only
!v8.51  Restore hs dependence (1:ns, not just ns)
!     fnorm1 = one/SUM(rzl_array(1:ns,:,m1:,1:2*ntmax)**2)    

!
!     ORIGIN EXTRAPOLATION (JS=1) FOR M=1 MODES
!     NOTE: PREVIOUS VERSIONS OF VMEC USED TWO-POINT EXTRAPOLATION 
!           FOR R,Z. HOWEVER,THIS CAN NOT BE USED TO COMPUTE THE 
!           TRI-DIAG 2D PRECONDITIONER
!
      rzl_array(1,:,m1,:)  = rzl_array(2,:,m1,:)
#ifdef _HBANGLE
!      rzl_array(1,:,m1+2,1:2*ntmax)= rzl_array(2,:,m1+2,1:2*ntmax)
#endif
      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)
      nsz = ns*nzeta
!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)
 
      IF (ictrl_prec2d .eq. 3) RETURN
     
      ALLOCATE (work1(nsz,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

!      r11 = 0;  ru1 = 0;  rv1 = 0;  rcn1 = 0
!      z11 = 0;  zu1 = 0;  zv1 = 0;  zcn1 = 0
!      lu1 = 0;  lv1 = 0
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED) 
!$OMP DO PRIVATE(i,l,lb,ub)
         DO i = 1, ntheta2
            l = nsz*(i-1)
            lb = l+1;  ub = nsz+l
            r11(lb:ub,:) = 0;  ru1(lb:ub,:) = 0
            rv1(lb:ub,:) = 0;  rcn1(lb:ub,:) = 0
            z11(lb:ub,:) = 0;  zu1(lb:ub,:) = 0
            zv1(lb:ub,:) = 0;  zcn1(lb:ub,:) = 0
            lu1(lb:ub,:) = 0;  lv1(lb:ub,:) = 0
         END DO
!$OMP END DO
!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
!         work1 = 0
!         j1 = jmin1(m)
#ifdef _HBANGLE
!         IF (m .EQ. 3) j1 = jmin1(1)
#endif
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
!$OMP DO PRIVATE(k,l,lb,ub,n,ni)
         DO k = 1, nzeta
            l = ns*(k-1)
            lb = l+1;  ub = ns+l
            work1(lb:ub,:) = 0
            DO n = 0, ntor
               ni = n+ioff
               work1(lb:ub,1) = work1(lb:ub,1) 
     1                        + rmncc(:,ni,mj)*cosnv(k,n)
               work1(lb:ub,6) = work1(lb:ub,6)
     1                        + zmnsc(:,ni,mj)*cosnv(k,n)
               work1(lb:ub,10) = work1(lb:ub,10) 
     1                        + lmnsc(:,ni,mj)*cosnv(k,n)
               IF (.NOT.lthreed) CYCLE
               
               work1(lb:ub,4) = work1(lb:ub,4) 
     1                        + rmnss(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,7) = work1(lb:ub,7) 
     1                        + zmncs(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,11) = work1(lb:ub,11)
     1                        + lmncs(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,2) = work1(lb:ub,2) 
     1                        + rmnss(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,5) = work1(lb:ub,5) 
     1                        + zmncs(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,9) = work1(lb:ub,9) 
     1                        + lmncs(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,3) = work1(lb:ub,3) 
     1                        + rmncc(:,ni,mj)*sinnvn(k,n)
               work1(lb:ub,8) = work1(lb:ub,8) 
     1                        + zmnsc(:,ni,mj)*sinnvn(k,n)
               work1(lb:ub,12) = work1(lb:ub,12)
     1                        + lmnsc(:,ni,mj)*sinnvn(k,n)
            END DO
         END DO
!$OMP END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!

!$OMP DO PRIVATE(i,l,lb,ub,cosmux,sinmux)
         DO i = 1, ntheta2
            l = nsz*(i-1)
            lb = l+1;  ub = nsz+l
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(lb:ub,mparity)  = r11(lb:ub,mparity)  
     1                          + work1(:,1)*cosmu(i,m)
            ru1(lb:ub,mparity)  = ru1(lb:ub,mparity)  
     1                          + work1(:,1)*sinmum(i,m)
#ifndef _HBANGLE
            rcn1(lb:ub,mparity) = rcn1(lb:ub,mparity) 
     1                          + work1(:,1)*cosmux
#endif
            z11(lb:ub,mparity)  = z11(lb:ub,mparity)  
     1                          + work1(:,6)*sinmu(i,m)

            zu1(lb:ub,mparity)  = zu1(lb:ub,mparity)  
     1                          + work1(:,6)*cosmum(i,m)
#ifndef _HBANGLE
            zcn1(lb:ub,mparity) = zcn1(lb:ub,mparity) 
     1                          + work1(:,6)*sinmux
#endif             
            lu1(lb:ub,mparity)  = lu1(lb:ub,mparity)  
     1                          + work1(:,10)*cosmum(i,m)

            IF (.NOT.lthreed) CYCLE

            r11(lb:ub,mparity)  = r11(lb:ub,mparity)  
     1                          + work1(:,2)*sinmu(i,m)
            ru1(lb:ub,mparity)  = ru1(lb:ub,mparity)  
     1                          + work1(:,2)*cosmum(i,m)
#ifndef _HBANGLE
            rcn1(lb:ub,mparity) = rcn1(lb:ub,mparity) 
     1                          + work1(:,2)*sinmux
#endif
            rv1(lb:ub,mparity)  = rv1(lb:ub,mparity)  
     1                          + work1(:,3)*cosmu(i,m) 
     1                          + work1(:,4)*sinmu(i,m)
            z11(lb:ub,mparity)  = z11(lb:ub,mparity)  
     1                          + work1(:,5)*cosmu(i,m)

            zu1(lb:ub,mparity)  = zu1(lb:ub,mparity)  
     1                          + work1(:,5)*sinmum(i,m)
#ifndef _HBANGLE
            zcn1(lb:ub,mparity) = zcn1(lb:ub,mparity) 
     1                          + work1(:,5)*cosmux
#endif
            zv1(lb:ub,mparity)  = zv1(lb:ub,mparity)  
     1                          + work1(:,7)*cosmu(i,m) 
     1                          + work1(:,8)*sinmu(i,m)

            lu1(lb:ub,mparity)  = lu1(lb:ub,mparity)  
     1                          + work1(:,9)*sinmum(i,m)
            lv1(lb:ub,mparity)  = lv1(lb:ub,mparity)  
     1                          - (work1(:,11)*cosmu(i,m) 
     1                          +  work1(:,12)*sinmu(i,m))
         END DO
!$OMP END DO
      END DO
!$OMP END PARALLEL
      DEALLOCATE (work1)

      z01(1:ns) = zmnsc(1:ns,n0+ioff,m1+joff)
      r01(1:ns) = rmncc(1:ns,n0+ioff,m1+joff)
      IF (r01(1) .EQ. zero) STOP 'r01(0) = 0 in totzsps'
      dkappa = z01(1)/r01(1)

      END SUBROUTINE totzsps

#else
      SUBROUTINE totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d, ONLY: ictrl_prec2d
      USE xstuff, ONLY: xc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns,nzeta,ntheta3,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, d2, s
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,6,ns):: myrzl
      REAL(rprec), DIMENSION(nzeta,ntheta3,0:1,ns):: myr11,myru1,
     1     myrv1, myz11, myzu1,  myzv1, mylu1, mylv1, myrcn1, myzcn1
      REAL(rprec), DIMENSION(nzeta,2):: mywork
      REAL(rprec), DIMENSION(nzeta,0:mpol1,2):: mywork2
      REAL(rprec), DIMENSION(nzeta) :: mytemp

C-----------------------------------------------
!
!     WORK1 Array of inverse transforms in toroidal angle (zeta), for all radial positions
!     NOTE: ORDERING OF LAST INDEX IS DIFFERENT HERE THAN IN PREVIOUS VMEC2000 VERSIONS
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(stored at rss) = .5(rss + zcs),
!     R-(stored at zcs) = .5(rss - zcs), TO EXTERNAL ("PHYSICAL") rss, zcs FORMS. NEED THIS EVEN 
!     WHEN COMPUTING HESSIAN FOR FREE BOUNDARY (rmnss, zmncs at JS=NS needed in vacuum call)
!
!     WHEN COMPUTING PRECONDITIONER, USE FASTER HESSIAN VERSION (totzsps_hess) INSTEAD.
      if (.not.lthreed) then
         print*,'lthreed=.false. is is not implemented in new totzps'
         print*,'please use old totzps implementation'
         stop
      endif
      if (lasym) then
         print*,'lasym=.true. is is not yet implemented in new totzps'
         print*,'please use old totzps implementation'
         stop
      endif

      if (first) then
         allocate(mycosnv(0:ntor,1:nzeta),mysinnv(0:ntor,1:nzeta),
     1        mycosnvn(0:ntor,1:nzeta),mysinnvn(0:ntor,1:nzeta))
         do n = 0, ntor
            do k = 1, nzeta
               mycosnv(n,k)=cosnv(k,n)
               mysinnv(n,k)=sinnv(k,n)
               mycosnvn(n,k)=cosnvn(k,n)
               mysinnvn(n,k)=sinnvn(k,n)
            end do
         end do
         first=.false.
         d3=ntor+1
         d4=mpol1+1
         eins=1.
         nix=0.
!         print*, nzeta,d3,(ntor+1),(mpol1+1)
      end if      

      CALL convert_sym (rzl_array(:,:,:,2), rzl_array(:,:,:,4))
 
!
!     ORIGIN EXTRAPOLATION (JS=1) FOR M=1 MODES
!     NOTE: PREVIOUS VERSIONS OF VMEC USED TWO-POINT EXTRAPOLATION 
!           FOR R,Z. HOWEVER,THIS CAN NOT BE USED TO COMPUTE THE 
!           TRI-DIAG 2D PRECONDITIONER
!
      rzl_array(1,:,m1,:)  = rzl_array(2,:,m1,:)

!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
!      call perfon('totzsps2')
      IF (lthreed .and. jlam(m0).gt.1) 
     1   rzl_array(1,:,m0,6) = rzl_array(2,:,m0,6)
 
      IF (ictrl_prec2d .eq. 3) RETURN

      z01(1:ns) = rzl_array(1:ns,n0,m1,3)
      r01(1:ns) = rzl_array(1:ns,n0,m1,1)
      IF (r01(1) .eq. zero) STOP 'r01(0) = 0 in totzsps'
      dkappa = z01(1)/r01(1)

      d2=(ntor+1)*(mpol1+1)*6

!      call perfon('trans')
!      call mkl_domatcopy('C','T',ns,d2,eins,rzl_array,ns,myrzl,d2)
      call mytrans(rzl_array,myrzl,d2,ns)

!      call perfoff

!      call perfon('totzsploop')

!$OMP PARALLEL DO NUM_THREADS(OMP_NUM_THREADS)
!$OMP&  DEFAULT(NONE)
!$OMP&  PRIVATE(s,m,i,mparity,mywork,mywork2,mytemp)
!$OMP&  SHARED(mpol1,nzeta,ns,ntheta2,ntor,mycosnv,mysinnv,
!$OMP&         mycosnvn,mysinnvn,cosmu,sinmu,cosmum,sinmum,xmpq,
!$OMP&         myrzl,myr11,myru1,myrv1,myz11,myzu1,myzv1,mylu1,
!$OMP&         mylv1,myrcn1,myzcn1,d3,eins,nix,d4)
      do s=2,ns
!         call perfon('innen')
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnv,d3,
     1        myrzl(0,0,2,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnv,d3,
     1        myrzl(0,0,1,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         myr11(:,:,:,s) = 0  
         myru1(:,:,:,s) = 0  
         myrcn1(:,:,:,s) = 0
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               mytemp = mywork2(:,m,1)*cosmu(i,m)
     1              + mywork2(:,m,2)*sinmu(i,m)
               myr11(:,i,mparity,s) = myr11(:,i,mparity,s) + mytemp
               myrcn1(:,i,mparity,s) = myrcn1(:,i,mparity,s) 
     1              + mytemp * xmpq(m,1)
            END DO
            DO i = 1, ntheta2
               myru1(:,i,mparity,s)  = myru1(:,i,mparity,s)
     1              + mywork2(:,m,1)*sinmum(i,m)
     2              + mywork2(:,m,2)*cosmum(i,m)
            END DO
         END DO    
!         call perfoff
         
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnvn,d3,
     1        myrzl(0,0,2,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnvn,d3,
     1        myrzl(0,0,1,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         myrv1(:,:,:,s) = 0.
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               myrv1(:,i,mparity,s)  = myrv1(:,i,mparity,s)
     1              + mywork2(:,m,1)*cosmu(i,m) 
     1              + mywork2(:,m,2)*sinmu(i,m)
            END DO
         END DO    


         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnv,d3,
     1        myrzl(0,0,4,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnv,d3,
     1        myrzl(0,0,3,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         myz11(:,:,:,s) = 0
         myzu1(:,:,:,s) = 0
         myzcn1(:,:,:,s) = 0
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               mytemp(:)  = mywork2(:,m,2)*sinmu(i,m)
     1              + mywork2(:,m,1)*cosmu(i,m)
               myz11(:,i,mparity,s)  = myz11(:,i,mparity,s) + mytemp
               myzcn1(:,i,mparity,s) = myzcn1(:,i,mparity,s) 
     1              + mytemp * xmpq(m,1)
            END DO
            DO i = 1, ntheta2
               myzu1(:,i,mparity,s)  = myzu1(:,i,mparity,s)
     1              + mywork2(:,m,2)*cosmum(i,m)
     2              + mywork2(:,m,1)*sinmum(i,m)
            END DO
         END DO            


         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnvn,d3,
     1        myrzl(0,0,4,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnvn,d3,
     1        myrzl(0,0,3,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         myzv1(:,:,:,s) = 0 
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               myzv1(:,i,mparity,s)  = myzv1(:,i,mparity,s)
     1              + mywork2(:,m,1)*cosmu(i,m) 
     1              + mywork2(:,m,2)*sinmu(i,m)
            END DO
         END DO            


         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnv,d3,
     1        myrzl(0,0,6,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnv,d3,
     1        myrzl(0,0,5,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         mylu1(:,:,:,s) = 0  
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               mylu1(:,i,mparity,s)  = mylu1(:,i,mparity,s)
     1              + mywork2(:,m,2)*cosmum(i,m)
     2              + mywork2(:,m,1)*sinmum(i,m)
            END DO
         END DO

         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mycosnvn,d3,
     1        myrzl(0,0,6,s),d3,
     1        nix,mywork2(1,0,1),nzeta)
         call dgemm('T','N',nzeta,d4,d3,eins,
     1        mysinnvn,d3,
     1        myrzl(0,0,5,s),d3,
     1        nix,mywork2(1,0,2),nzeta)
         mylv1(:,:,:,s) = 0  
         DO m = 0, mpol1
            mparity = MOD(m,2)
            DO i = 1, ntheta2
               mylv1(:,i,mparity,s)  = mylv1(:,i,mparity,s)
     1              - (mywork2(:,m,1)*cosmu(i,m) 
     1              +  mywork2(:,m,2)*sinmu(i,m))
            END DO
         END DO

      END DO
!$OMP END PARALLEL DO

!      call perfon ('snull')
      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mycosnv,d3,myrzl(0,m,1,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mysinnv,d3,myrzl(0,m,2,1),
     1        1,nix,mywork(1,2),1)
         DO i = 1, ntheta2
            mytemp = mywork(:,1)*cosmu(i,m)
     1           + mywork(:,2)*sinmu(i,m)
            myr11(:,i,mparity,1) = mytemp
            myrcn1(:,i,mparity,1) = mytemp * xmpq(m,1)
         END DO
         DO i = 1, ntheta2
            myru1(:,i,mparity,1) = mywork(:,1)*sinmum(i,m)
     2           + mywork(:,2)*cosmum(i,m)
         END DO
      END DO    

      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mysinnvn,d3,myrzl(0,m,1,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mycosnvn,d3,myrzl(0,m,2,1),
     1        1,nix,mywork(1,2),1)
         DO i = 1, ntheta2
            myrv1(:,i,mparity,1) = mywork(:,1)*cosmu(i,m) 
     1           + mywork(:,2)*sinmu(i,m)
         END DO
      END DO    
      
      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mysinnv,d3,myrzl(0,m,4,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mycosnv,d3,myrzl(0,m,3,1),
     1        1,nix,mywork(1,2),1)
         DO i = 1, ntheta2
            mytemp(:)  = mywork(:,2)*sinmu(i,m)
     1           + mywork(:,1)*cosmu(i,m)
            myz11(:,i,mparity,1)  = mytemp
            myzcn1(:,i,mparity,1) = mytemp * xmpq(m,1)
         END DO
         DO i = 1, ntheta2
            myzu1(:,i,mparity,1)  = mywork(:,2)*cosmum(i,m)
     2           + mywork(:,1)*sinmum(i,m)
         END DO
      END DO            
      
      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mycosnvn,d3,myrzl(0,m,4,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mysinnvn,d3,myrzl(0,m,3,1),
     1        1,nix,mywork(1,2),1)
         DO i = 1, ntheta2
            myzv1(:,i,mparity,1) = mywork(:,1)*cosmu(i,m) 
     1           + mywork(:,2)*sinmu(i,m)
         END DO
      END DO            
      
      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mysinnv,d3,myrzl(0,m,6,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mycosnv,d3,myrzl(0,m,5,1),
     1        1,nix,mywork(1,2),1)
         
         DO i = 1, ntheta2
            mylu1(:,i,mparity,1) = mywork(:,2)*cosmum(i,m)
     2           + mywork(:,1)*sinmum(i,m)
         END DO
      END DO
      
      DO m = 0, 1
         mparity = MOD(m,2)
         call dgemv('T',d3,nzeta,eins,mycosnvn,d3,myrzl(0,m,6,1),
     1        1,nix,mywork(1,1),1)
         call dgemv('T',d3,nzeta,eins,mysinnvn,d3,myrzl(0,m,5,1),
     1        1,nix,mywork(1,2),1)
         DO i = 1, ntheta2
            mylv1(:,i,mparity,1) = - mywork(:,1)*cosmu(i,m) 
     1           - mywork(:,2)*sinmu(i,m)
         END DO
      END DO
!      call perfoff
!      call perfoff

      d2=nzeta*ntheta3*2
!      call perfon('trans2')
!$OMP PARALLEL SECTIONS NUM_THREADS(OMP_NUM_THREADS)
!$OMP SECTION
      call mytrans(myr11,r11,ns,d2)
!$OMP SECTION
      call mytrans(myru1,ru1,ns,d2)
!$OMP SECTION
      call mytrans(myrcn1,rcn1,ns,d2)
!$OMP SECTION
      call mytrans(myrv1,rv1,ns,d2)
!$OMP SECTION
      call mytrans(myz11,z11,ns,d2)
!$OMP SECTION
      call mytrans(myzu1,zu1,ns,d2)
!$OMP SECTION
      call mytrans(myzcn1,zcn1,ns,d2)
!$OMP SECTION
      call mytrans(myzv1,zv1,ns,d2)
!$OMP SECTION
      call mytrans(mylu1,lu1,ns,d2)
!$OMP SECTION
      call mytrans(mylv1,lv1,ns,d2)
!$OMP END PARALLEL SECTIONS

!      call perfoff
!      call perfoff

      END SUBROUTINE totzsps

      subroutine mytrans(m_in,m_out,dim1,dim2)
      integer:: dim1,dim2
      REAL(rprec), DIMENSION(dim2,dim1):: m_in
      REAL(rprec), DIMENSION(dim1,dim2):: m_out
      REAL(rprec), PARAMETER:: one=1
#if 1
      m_out=transpose(m_in)
#else
      call mkl_domatcopy('C','T',dim2,dim1,one,m_in,dim2,m_out,dim1)
#endif      
      end subroutine mytrans
#endif

 
      SUBROUTINE totzsps_hess(rzl_array, r11, ru1, rv1, z11, zu1, 
     1                        zv1, lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
      LOGICAL :: logl, logr, logz
C-----------------------------------------------
!
!     SAME AS totzsps, BUT ONLY COMPUTES perturbation (rzl_array is perturbation) TO EXISTING R,Z FOR
!     A PARTICULAR (m,n,ntype) = (m_2d,n_2d,ntype_2d) VALUE
!
!     
!     STORE (FOR ictrl_prec2d=2) OR RESTORE (ictrl_prec2d=3) VARIABLES NOT RECOMPUTED 
!     IN TOTZSP DURING HESSIAN JOGS
!
      IF (ictrl_prec2d .eq. 2) THEN
         lus_save = lu1;  lvs_save = lv1;
         rus_save = ru1;  rvs_save = rv1; r1s_save = r11
         rcons_save = rcn1
         zus_save = zu1;  zvs_save = zv1; z1s_save = z11
         zcons_save = zcn1
         RETURN
      ELSE
         lu1 = lus_save;  lv1 = lvs_save
         ru1 = rus_save;  rv1 = rvs_save; r11 = r1s_save
         rcn1 = rcons_save
         zu1 = zus_save;  zv1 = zvs_save; z11 = z1s_save
         zcn1 = zcons_save
      END IF

      logr = ntype_2d .le. ntmax
      logz = (ntype_2d.gt.ntmax) .and. (ntype_2d.le.2*ntmax)
      logl = ntype_2d .gt. 2*ntmax

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)            !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)      !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)    !!COS(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

!
!     ENFORCE CONSTRAINT ON m=1 MODES FOR 3D (CONSISTENT WITH gcz(zcs) = 0 IN RESIDUE)
!     NOTE: Since r,z variations are coupled, must vary BOTH whenever ONE is varied
!

      IF (lthreed .and. m_2d.eq.1 .and. n_2d.ne.0 .and. 
     1   (ntype_2d.eq.rss .or. ntype_2d.eq.(zcs+ntmax))) THEN
!         zmncs(:,:,m1+joff) = 0
!         IF (logz) logz = .false.
!         IF (logr .and. lconm1) logz = .true.
         logr = .true.;  logz = .true.
         CALL convert_sym (rmnss, zmncs)
      END IF

!
!     EXTRAPOLATION AT JS=1 FOR M=1 MODES
!
      rzl_array(1,:,m1,:) = rzl_array(2,:,m1,:)   

!
!     EXTRAPOLATION AT JS=1 OF M=0 MODES FOR LAMBDA
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
         m = m_2d
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         n = n_2d
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               IF (logr)
     1         work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     2                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logz)
     1         work1(j1l:nsl,6) = work1(j1l:nsl,6)
     2                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logl)
     1         work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     2                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
     
               IF (logr) THEN          
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logz) THEN
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logl) THEN
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

            END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            cosmux
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            sinmux
            END IF
 
            IF (logl)
     1      lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     2            cosmum(i,m)

            IF (.not.lthreed) CYCLE

            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1            cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1            sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1            cosmu(i,m) + work1(:,4)*sinmu(i,m)
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1            sinmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1            cosmux
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1            cosmu(i,m) + work1(:,8)*sinmu(i,m)
            END IF

            IF (logl) THEN
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1            sinmum(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1            cosmu(i,m) + work1(:,12)*sinmu(i,m))
            END IF

         END DO

      DEALLOCATE (work1)

      END SUBROUTINE totzsps_hess

      SUBROUTINE convert_sym(rmnss, zmncs)
      USE vmec_main, p5 => cp5
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), INTENT(inout) :: 
     1                                           rmnss, zmncs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO "PHYSICAL" RMNSS, ZMNCS FOURIER FORM
!
      IF (.NOT.lconm1) RETURN
!$OMP PARALLEL WORKSHARE NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
      temp(:,:) = rmnss(:,:,1)
      rmnss(:,:,1) = temp(:,:) + zmncs(:,:,1)
      zmncs(:,:,1) = temp(:,:) - zmncs(:,:,1)
!$OMP END PARALLEL WORKSHARE

      END SUBROUTINE convert_sym


      SUBROUTINE totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1,
     1   lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1, lb, ub, nsz
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
C-----------------------------------------------
!     WHEN COMPUTING PRECONDITIONER, USE FASTER HESSIAN VERSION (totzsps_hess) INSTEAD.

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      CALL convert_asym (rmnsc, zmncc)

      IF (ictrl_prec2d .EQ. 3) RETURN

      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)
      nsz = ns*nzeta

      z00b = zmncc(ns,ioff,joff)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC totzspa'

!
!     INITIALIZATION BLOCK
!
!      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
!      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0
!$OMP PARALLEL NUM_THREADS(OMP_NUM_THREADS) DEFAULT(SHARED)
!$OMP DO PRIVATE(i,l,lb,ub)
         DO i = 1, ntheta2
            l = nsz*(i-1)
            lb = l+1;  ub = nsz+l
            r11(lb:ub,:) = 0;  ru1(lb:ub,:) = 0
            rv1(lb:ub,:) = 0;  rcn1(lb:ub,:) = 0
            z11(lb:ub,:) = 0;  zu1(lb:ub,:) = 0
            zv1(lb:ub,:) = 0;  zcn1(lb:ub,:) = 0
            lu1(lb:ub,:) = 0;  lv1(lb:ub,:) = 0
         END DO
!$OMP END DO

      IF (jlam(m0) .gt. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
!         work1 = 0
!         j1 = jmin1(m)
!$OMP DO PRIVATE(k,l,lb,ub,n,ni)
         DO k = 1, nzeta
            l = ns*(k-1)
            lb = l+1;  ub = ns+l
            work1(lb:ub,:) = 0

            DO n = 0, ntor
               ni = n+ioff
               work1(lb:ub,1) = work1(lb:ub,1)
     1                          + rmnsc(:,ni,mj)*cosnv(k,n)
               work1(lb:ub,6) = work1(lb:ub,6) 
     1                          + zmncc(:,ni,mj)*cosnv(k,n)
               work1(lb:ub,10) = work1(lb:ub,10)
     1                          + lmncc(:,ni,mj)*cosnv(k,n)

               IF (.NOT.lthreed) CYCLE

               work1(lb:ub,2) = work1(lb:ub,2)
     1                          + rmncs(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,3) = work1(lb:ub,3) 
     1                          + rmnsc(:,ni,mj)*sinnvn(k,n)
               work1(lb:ub,4) = work1(lb:ub,4) 
     1                          + rmncs(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,5) = work1(lb:ub,5) 
     1                          + zmnss(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,7) = work1(lb:ub,7)
     1                          + zmnss(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,8) = work1(lb:ub,8) 
     1                          + zmncc(:,ni,mj)*sinnvn(k,n)
               work1(lb:ub,9) = work1(lb:ub,9) 
     1                          + lmnss(:,ni,mj)*sinnv(k,n)
               work1(lb:ub,11) = work1(lb:ub,11)
     1                          + lmnss(:,ni,mj)*cosnvn(k,n)
               work1(lb:ub,12) = work1(lb:ub,12) 
     1                          + lmncc(:,ni,mj)*sinnvn(k,n)
            END DO
         END DO
!$OMP END DO
!
!        INVERSE TRANSFORM IN M-THETA
!
!$OMP DO PRIVATE(i,l,lb,ub,cosmux,sinmux)
         DO i = 1, ntheta2
            l = nsz*(i-1)
            lb = l+1;  ub = nsz+l
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            lb = l+1;  ub = nsz+l
            r11(lb:ub,mparity) = r11(lb:ub,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(lb:ub,mparity) = ru1(lb:ub,mparity) + work1(:,1)*
     1            cosmum(i,m)
            z11(lb:ub,mparity) = z11(lb:ub,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(lb:ub,mparity) = zu1(lb:ub,mparity) + work1(:,6)*
     1            sinmum(i,m)
            lu1(lb:ub,mparity) = lu1(lb:ub,mparity) + work1(:,10)*
     1            sinmum(i,m)
            rcn1(lb:ub,mparity) = rcn1(lb:ub,mparity) + work1(:,1)*
     1            sinmux
            zcn1(lb:ub,mparity) = zcn1(lb:ub,mparity) + work1(:,6)*
     1            cosmux

            IF (.NOT.lthreed) CYCLE
               
            r11(lb:ub,mparity) = r11(lb:ub,mparity) + work1(:,2)*
     1               cosmu(i,m)
            ru1(lb:ub,mparity) = ru1(lb:ub,mparity) + work1(:,2)*
     1               sinmum(i,m)
            z11(lb:ub,mparity) = z11(lb:ub,mparity) + work1(:,5)*
     1               sinmu(i,m)
            zu1(lb:ub,mparity) = zu1(lb:ub,mparity) + work1(:,5)*
     1               cosmum(i,m)
            lu1(lb:ub,mparity) = lu1(lb:ub,mparity) + work1(:,9)*
     1               cosmum(i,m)
            rcn1(lb:ub,mparity) = rcn1(lb:ub,mparity) + work1(:,2)*
     1               cosmux
            zcn1(lb:ub,mparity) = zcn1(lb:ub,mparity) + work1(:,5)*
     1               sinmux
            rv1(lb:ub,mparity) = rv1(lb:ub,mparity) + work1(:,3)*
     1               sinmu(i,m) + work1(:,4)*cosmu(i,m)
            zv1(lb:ub,mparity) = zv1(lb:ub,mparity) + work1(:,7)*
     1               sinmu(i,m) + work1(:,8)*cosmu(i,m)
            lv1(lb:ub,mparity) = lv1(lb:ub,mparity) - (work1(:,11)*
     1               sinmu(i,m)+work1(:,12)*cosmu(i,m))
         END DO
!$OMP END DO
      END DO
!$OMP END PARALLEL
      DEALLOCATE (work1)

      END SUBROUTINE totzspa

      SUBROUTINE totzspa_hess(rzl_array, r11, ru1, rv1, z11, zu1, 
     1                        zv1, lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
      LOGICAL :: logl, logr, logz
C-----------------------------------------------
!
!     SAME AS totzspa, BUT ONLY COMPUTES perturbation (rzl_array is perturbation) TO EXISTING R,Z, L FOR
!     A PARTICULAR (m,n,ntype) = (m_2d,n_2d,ntype_2d) VALUE
!
!     
!     STORE (FOR ictrl_prec2d=2) OR RESTORE (ictrl_prec2d=3) VARIABLES NOT RECOMPUTED 
!     IN TOTZSP DURING HESSIAN JOGS
!
      IF (ictrl_prec2d .eq. 2) THEN
         lua_save = lu1;  lva_save = lv1;
         rua_save = ru1;  rva_save = rv1; r1a_save = r11 
         rcona_save = rcn1
         zua_save = zu1;  zva_save = zv1; z1a_save = z11 
         zcona_save = zcn1
         RETURN
      ELSE
         lu1 = lua_save;  lv1 = lva_save
         ru1 = rua_save;  rv1 = rva_save; r11 = r1a_save 
         rcn1 = rcona_save
         zu1 = zua_save;  zv1 = zva_save; z11 = z1a_save
         zcn1 = zcona_save
      END IF

      logr = ntype_2d .le. ntmax
      logz = (ntype_2d.gt.ntmax) .and. (ntype_2d.le.2*ntmax)
      logl = ntype_2d .gt. 2*ntmax

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)

!
!     ENFORCE CONSTRAINT ON m=1 MODES FOR 3D (CONSISTENT WITH gcz(zcc) = 0 IN RESIDUE)
!     NOTE: Since r,z variations are coupled, must vary BOTH whenever ONE is varied
!
      IF (m_2d.eq.1 .and. 
     1   (ntype_2d.eq.rsc .or. ntype_2d.eq.(ntmax+zcc))) THEN
!         zmncc(:,:,m1+joff) = 0
!         IF (logz) logz = .false.
!         IF (logr .and. lconm1) logz = .true.
         logr = .true.;  logz = .true.
         CALL convert_asym (rmnsc, zmncc)
      END IF


      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'


!     EXTRAPOLATION OF M=0 MODES FOR LAMBDA
!
      IF (jlam(m0) .gt. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
         m = m_2d
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         n = n_2d
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               IF (logr)
     1         work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     2                          + rmnsc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logz)
     1         work1(j1l:nsl,6) = work1(j1l:nsl,6)
     2                          + zmncc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logl)
     1         work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     2                          + lmncc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
     
               IF (logr) THEN          
               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmncs(j1:ns,ni,mj)*cosnvn(k,n)
               END IF

               IF (logz) THEN
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logl) THEN
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

            END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux
            END IF
 
            IF (logl)
     1      lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     2            sinmum(i,m)

            IF (.not.lthreed) CYCLE

            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1            cosmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1            sinmu(i,m) + work1(:,4)*cosmu(i,m)
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1            sinmux
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1            sinmu(i,m) + work1(:,8)*cosmu(i,m)
            END IF

            IF (logl) THEN
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1            cosmum(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1            sinmu(i,m) + work1(:,12)*cosmu(i,m))
            END IF

         END DO
      
      DEALLOCATE (work1)

      END SUBROUTINE totzspa_hess


      SUBROUTINE convert_asym(rmnsc, zmncc)
      USE vmec_main, p5_2 => cp5
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), INTENT(inout) :: 
     1                                           rmnsc, zmncc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO RMNSC, ZMNCC FOURIER FORM
!
      IF (.not.lconm1) RETURN

      temp(:,:) = rmnsc(:,:,1)
      rmnsc(:,:,1) = temp(:,:) + zmncc(:,:,1)
      zmncc(:,:,1) = temp(:,:) - zmncc(:,:,1)

      END SUBROUTINE convert_asym

      END MODULE totzsp_mod
