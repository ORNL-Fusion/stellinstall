    MODULE pchelms
    USE stel_kinds
    USE metrics, ONLY: init_metric_elements
    USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i
    USE island_params, ntheta=>nu_i, nzeta=>nv_i, hs=>hs_i
    USE island_params, ohs => ohs_i, mnmax=>mnmax_i, nfp=>nfp_i
    USE hessian, ONLY: initHess, Compute_Hessian_Blocks, asym_index, levmarq_param
    USE blocktridiagonalsolver_s, ONLY: Initialize, GetRanks
    USE nscalingtools, ONLY: MyEnvVariables, PARSOLVER, disp, startglobrow, endglobrow, rcounts, MPI_ERR
    USE shared_data, ONLY: gc, xc, neqs, gc0, xc0, fsq_total, mblk_size, ndims
    USE shared_data, ONLY: asubsmnsf, asubumncf, asubvmncf,             &
      fsupsmnsf, fsupumncf, fsupvmncf,                                  &
      asubsmncf, asubumnsf, asubvmnsf,                                  &
      fsupsmncf, fsupumnsf, fsupvmnsf
    USE shared_data, ONLY: torflux, polflux, nsp, r1_i, z1_i, ru_i,     &
                           zu_i, rv_i, zv_i
    USE descriptor_mod, ONLY: iam, nprocs, lscalapack, INHESSIAN, SIESTA_COMM
    USE v3_utilities, ONLY: assert

    USE vmec_info, ONLY: jcurrumnc, jcurrvmnc, jcurrumns, jcurrvmns
    USE fourier
    USE utilities, ONLY: to_half_mesh, to_full_mesh
    USE metrics, ONLY: tolowerh
    USE quantities, ONLY: jacobh, jbsupsmnsh, jbsupumnch, jbsupvmnch,   &
                          jbsupsmnch, jbsupumnsh, jbsupvmnsh

    USE mpi_inc
    USE diagnostics_mod, ONLY: divb, divb_rms
    IMPLICIT NONE

    INTEGER, PRIVATE :: nsmin, nsmax, ns_match

    !Added SKS
    LOGICAL, PRIVATE            :: linit, lHessian
    LOGICAL, PARAMETER, PRIVATE :: l_AsEdge = .TRUE.                    !evolve (T) or fix (F) A_s at edge
    INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1, m0=0, m1=1, n0=0
    REAL(dp), PRIVATE           :: bnorm=-1, line_Bu

    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), PRIVATE ::                 &
      Bsupsijh, Bsupuijh, Bsupvijh,                                     &
      Bsubsijh, Bsubuijh, Bsubvijh, jacobmnch, jacobmnsh

    CONTAINS

      SUBROUTINE run_pchelms
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!     INTEGER, PARAMETER :: ns_in=101, mpol_in=12, ntor_in=6    !read in from file eventually
      INTEGER            :: istat, js
      CHARACTER(LEN=*), PARAMETER ::                                    &
         banner='----------------------------------------------'
!	  CHARACTER(LEN=*), PARAMETER :: wout_file="wout_akima.nc"
!     CHARACTER(LEN=*), PARAMETER :: wout_file="wout_w7x_mark.nc"
!     CHARACTER(LEN=*), PARAMETER :: wout_file="wout_w7x.00sfix.nc"
      REAL(dp) :: ton, toff
      LOGICAL  :: flag
!-----------------------------------------------
      IF (iam .EQ. 0) THEN
         PRINT *
         PRINT *,banner
         PRINT *,'STARTING EXTERNAL B-FIELD CALCULATION (PCHELMS)'
         PRINT *,banner
      END IF

      CALL init_pchelms

      CALL second0(ton)
!IN SIESTA_RUN    CALL initHess
      levmarq_param = 1.E-6_dp   !add diffusive term to preconditioner instead

!    CALL Compute_Hessian_Blocks (Compute_Forces_Lin, .TRUE.)         !diag prec
      CALL Compute_Hessian_Blocks (Compute_Forces_Lin, .FALSE.)           !full prec
      CALL second0(toff)
      IF (iam .EQ. 0) THEN
         PRINT 90,' ASYMMETRY INDEX: ', asym_index
         PRINT 100,' PCHELMS HESSIAN COMPUTATION TIME: ', toff-ton, 's'
    END IF

      CALL second0(ton)
!    CALL CONJUGATE_GRADIENT_PCHELMS
      CALL GMRES_PCHELMS
      CALL second0(toff)

      IF (iam .EQ. 0) PRINT 100,' PCHELMS SOLVER TIME: ', toff-ton, 's'

!update solution (initial guess + perturbation) and check final forces
      CALL BACKSLV

!dump output: MUST call backslv to get updated B's before dumping
      DO js = 1, ns
         CALL DUMP_A(js, 333)
         CALL DUMP_B(js, 666)
      END DO

      CALL divb(nsmin,nsmax)

      IF (iam .EQ. 0) PRINT '(/,a,1p,e12.3)',' |DivB| (normed) = ', divb_rms
    
      IF (iam .EQ. 0) THEN
         PRINT *
         PRINT *,banner
         PRINT *,'ENDING PCHELMS'
         PRINT *,banner
      END IF
    

      CALL CHECK_CURRENT

  90  FORMAT(a,1p,e12.4)
 100  FORMAT(/,a,f10.1,a)

      CALL CLEANUP
 
      END SUBROUTINE run_pchelms


!>    \brief Compute real-space components of contravariant jac*B^i on 
!>     half radial grid from curl(A)
      SUBROUTINE CURLA_PCHELMS(Asubsmnf,  Asubumnf,  Asubvmnf,          &
                              jbsupsmnh, jbsupumnh, jbsupvmnh, ipar)
      USE utilities, ONLY: CURL_FtoH
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar
      REAL(dp), POINTER, DIMENSION(:,:,:) ::                            &
                Asubsmnf, Asubumnf, Asubvmnf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                jbsupsmnh, jbsupumnh, jbsupvmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fcomb, fours, fouru, fourv, moff, nmin, nmax, nmatch
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fcomb = f_none; fours = f_sin; fouru = f_cos; fourv = f_cos
      ELSE
         fcomb = f_sum;  fours = f_cos; fouru = f_sin; fourv = f_sin
      END IF

      moff = LBOUND(asubsmnf,1)
      nmin = nsmin;  nmax = nsmax
      
      IF (nsmin .EQ. 1) THEN
          asubsmnf(:,:,1) = 0 
          asubumnf(:,:,1) = 0;  asubvmnf(moff+m1:,:,1) = 0
      END IF

!Dirichelet bc (fixed) for A_u, A_v (A_s floats) consistent with Poynting flux = 0
      IF (nmax.EQ.ns .AND. lHessian .AND. .NOT.lInit) THEN
         IF (.NOT.l_AsEdge) asubsmnf(:,:,ns) = 0 
         asubumnf(:,:,ns) = 0; asubvmnf(:,:,ns) = 0
      END IF

!need this to keep Hessian symmetric      
      IF (nmin.LE.ns_match .AND. lHessian .AND. .NOT.lInit) THEN
          nmatch = MIN(ns_match, nmax)
          asubumnf(m1+moff,:,nmin:nmatch) = 0
          asubvmnf(m1+moff,:,nmin:nmatch) = 0
      END IF


      CALL CURL_FtoH(Asubsmnf,  Asubumnf,  Asubvmnf,                    &
                    jbsupsmnh, jbsupumnh, jbsupvmnh, ipar, nmin, nmax)
      nmin = startglobrow

!Never use this except to make B_X finite at s=1
      IF (nmin .EQ. 1) THEN
         jbsupsmnh(:,:,1) = 0; jbsupumnh(:,:,1) = 0; jbsupvmnh(:,:,1) = 0
      END IF

!
!     convert to real space
!
      CALL toijsp_par(jbsupsmnh(:,:,nmin:nmax), Bsupsijh, fcomb, fours)
      CALL toijsp_par(jbsupumnh(:,:,nmin:nmax), Bsupuijh, fcomb, fouru)
      CALL toijsp_par(jbsupvmnh(:,:,nmin:nmax), Bsupvijh, fcomb, fourv)

      
      END SUBROUTINE CURLA_PCHELMS
      

!>     \brief Compute real-space components of contravariant jac*J^i == K^i 
!>      on full radial grid from curl(B)
      SUBROUTINE CURLB_PCHELMS(ksupsmnf, ksupumnf, ksupvmnf, asubsmnf,  &
                               bsupsmnh, ipar)
      USE utilities, ONLY: CURL_HtoF, GradientHalf
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar
      REAL(dp), POINTER, DIMENSION(:,:,:) :: asubsmnf
      REAL(dp), POINTER, DIMENSION(:,:,:) ::                            &
                ksupsmnf, ksupumnf, ksupvmnf
      REAL(dp), ALLOCATABLE, INTENT(IN) :: bsupsmnh(:,:,:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: fours, fouru, fourv, moff, noff, n1, n2,              &
                  istat, nmin, nmax, mp, np, m, n, sparity
!----------------------------------------------- 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
         db_uds, db_vds, bsubsmnf, bsubumnf, bsubvmnf, dA_sds
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                bsubsmnh, bsubumnh, bsubvmnh
      REAL(dp) :: Diff, t1(nsp)
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fours = f_sin; fouru = f_cos; fourv = f_cos
         sparity = 1
      ELSE
         fours = f_cos; fouru = f_sin; fourv = f_sin
         sparity = -1
      END IF

      nmin = startglobrow;  nmax = MIN(ns, endglobrow+1)
      moff = LBOUND(ksupsmnf,1); noff = LBOUND(ksupsmnf,2) + ntor

      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,nmin:nmax),                   & 
               bsubumnh(0:mpol,-ntor:ntor,nmin:nmax),                   &
               bsubvmnh(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)   
      CALL ASSERT(istat.EQ.0,'Allocation1 failed in CURLB_PCHELMS')
      
      CALL tomnsp_par(bsubsijh, bsubsmnh(:,:,nmin:nmax), fours)          ! Compute BSUBXMNH    
      CALL tomnsp_par(bsubuijh, bsubumnh(:,:,nmin:nmax), fouru)
      CALL tomnsp_par(bsubvijh, bsubvmnh(:,:,nmin:nmax), fourv)
      
      IF (nmax.EQ.ns .AND. ipar.EQ.isym) line_Bu = bsubumnh(0,0,ns)

      CALL CURL_HtoF(Bsubsmnh, Bsubumnh, Bsubvmnh,                      &
                     Ksupsmnf, Ksupumnf, Ksupvmnf, ipar, nmin, nmax)
      nmax = endglobrow
      DEALLOCATE (bsubsmnh, bsubumnh, bsubvmnh, stat=istat)   
      CALL ASSERT(istat.EQ.0,'Deallocation failed in CURLB_PCHELMS')

 
!Add smoothing diffusion for A_s in Hessian: db_uds stores dA_s/ds on half grid
      IF (lHessian) THEN
         np = MAX(1,nmin-1);  mp = MIN(ns, nmax+1)
         ALLOCATE (dA_sds(0:mpol,-ntor:ntor,np:mp))
         CALL GradientHalf(dA_sds, asubsmnf)
         Diff = 1.E-4_dp*ohs
         ksupsmnf(:,:,nmin:mp-1) = ksupsmnf(:,:,nmin:mp-1) +  Diff*     &
                      (dA_sds(:,:,nmin+1:mp) - dA_sds(:,:,nmin:mp-1))
         IF (nmax .EQ. ns) THEN
            ksupsmnf(:,:,ns) = ksupsmnf(:,:,ns) - 2*Diff*dA_sds(:,:,ns)
         END IF
         IF (ipar .EQ. isym) ksupsmnf(m0+moff,n0+noff,nmin:nmax) = 0
         DEALLOCATE(dA_sds)
      END IF

      n2 = MIN(nmax, nsp)
      IF (n2 .GE. nmin) THEN
!Force b^s -> 0 for js <= nsp (add (B^s)**2 to energy)

! bsupsmnh(m,n,nmin:nmax) = -mp*asubvmnh(m,n,nmin:nmax) 
!                         +  np*asubumnh(m,n,nmin:nmax)
      n1 = nmin
      Diff = 0.1_dp
      DO m = 0,mpol
         mp = m*sparity
         DO n = -ntor,ntor
            np = n*sparity*nfp
            t1(n1:n2) = Diff*(bsupsmnh(m,n,n1:n2)                       &
                      +       bsupsmnh(m,n,n1+1:n2+1))                  
            ksupumnf(m+moff,n+noff,n1:n2) =                             &
            ksupumnf(m+moff,n+noff,n1:n2) - np*t1(n1:n2)    
            ksupvmnf(m+moff,n+noff,n1:n2) =                             & 
            ksupvmnf(m+moff,n+noff,n1:n2) + mp*t1(n1:n2)    
         END DO
      END DO
      END IF

!ORIGIN BOUNDARY CONDITIONS: IF asubs(m0,:,1) = 0 IN CURLA, MUST
!                            SET ksupsmnf(m0,:,1) = 0 HERE FOR HESSIAN SYMMETRY
!NEED FACTOR OF 1/2 FOR HESSIAN SYMMETRY AT JS=1
      IF (nmin .EQ. 1) THEN
         ksupsmnf(m0+moff,:,1) = 0*0.5_dp*ksupsmnf(m0+moff,:,1)
         ksupvmnf(m0+moff,:,1) = 0.5_dp*ksupvmnf(m0+moff,:,1)
         ksupsmnf(m1+moff:,:,1) = 0
         ksupumnf(:,:,1) = 0
         ksupvmnf(m1+moff:,:,1) = 0
      END IF
      IF (nmax .EQ. ns) THEN
         ksupsmnf(:,:,ns) = 0.5_dp*ksupsmnf(:,:,ns)
         ksupumnf(:,:,ns) = 0.5_dp*ksupumnf(:,:,ns)
         ksupvmnf(:,:,ns) = 0.5_dp*ksupvmnf(:,:,ns)
      END IF
      

!FLIP SIGN TO GET CORRECT DIAGONAL SIGN (>0)
      ksupsmnf(:,:,nmin:nmax) = -ksupsmnf(:,:,nmin:nmax)
      ksupumnf(:,:,nmin:nmax) = -ksupumnf(:,:,nmin:nmax)
      ksupvmnf(:,:,nmin:nmax) = -ksupvmnf(:,:,nmin:nmax)
     

      END SUBROUTINE CURLB_PCHELMS


      SUBROUTINE CYL2VMEC_A(A_r,A_p,A_z,cA_s,cA_u,cA_v)

      REAL(dp), DIMENSION(ntheta,nzeta), INTENT(IN)    ::  A_r, A_p, A_z
      REAL(dp), DIMENSION(ntheta,nzeta), INTENT(OUT)   ::  cA_s, cA_u, cA_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   ::  rs, zs
      INTEGER :: u, v, js

!      DO u=1, ntheta
!        DO v=1, nzeta
!          DO js=1, ns
!            WRITE(400,*) u, v, js, r1_i(u,v,js), z1_i(u,v,js)
!            CALL FLUSH(400)
!          END DO
!        END DO
!      END DO

      !print *, "NS in CYC2VMEC: ", ns, nsp
      ALLOCATE(rs(ntheta, nzeta))
      ALLOCATE(zs(ntheta, nzeta))
      rs=0
      zs=0
      rs = (r1_i(:,:,ns) - r1_i(:,:,ns-1))*ohs
      zs = (z1_i(:,:,ns) - z1_i(:,:,ns-1))*ohs
!      DO u=1, ntheta
!        DO v=1, nzeta
!          WRITE(300,*) u, v, rs(u,v), zs(u,v)
!          CALL FLUSH(300)
!        END DO
!      END DO
      cA_s(:,:)=A_r(:,:)*rs(:,:) + A_z(:,:)*zs(:,:)
      cA_u(:,:)=A_r(:,:)*ru_i(:,:,ns) + A_z(:,:)*zu_i(:,:,ns)
      cA_v(:,:)=A_r(:,:)*rv_i(:,:,ns) + A_z(:,:)*zv_i(:,:,ns) + A_p(:,:)*r1_i(:,:,ns)
      DEALLOCATE(rs,zs)

      END SUBROUTINE CYL2VMEC_A


      SUBROUTINE INIT_A (A_s, A_u, A_v, lmn, ipar)
      USE metrics, ONLY: phipf_i, chipf_i
      USE bmw_run, ONLY: bmw_exec
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ipar
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns)   ::  A_s, A_u, A_v
      REAL(dp), DIMENSION(ntheta,nzeta,1) ::  A_r, A_p, A_z
      REAL(dp), DIMENSION(ntheta,nzeta,1) ::  cA_s, cA_u, cA_v
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns+1), INTENT(IN) ::  lmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, istat, m, n, u, v, nmatch
      INTEGER   :: fours, fouru, fourv
      REAL(dp) :: xm, xn, du, dv, t1, p1, q1
      REAL(dp) :: dphi, s_edge
      CHARACTER(LEN=256) :: LABEL
!-----------------------------------------------
!Note: lmns is on half grid here (splined to ns+1 in metrics)
!Only do this once, so no need to parallelize

      A_s = 0; A_u = 0;  A_v = 0
      
#if defined(NOSKS)
    !XXX: This is just a test...
      OPEN(unit=336, file="A_sA_uA_v.txt", status='old', iostat=istat)
      IF (istat .EQ. 0) THEN

         READ(336, *) LABEL
         IF (iam .EQ. 0) PRINT *,'Read in A_X from Biot Savart boundary file'

         DO WHILE (istat .EQ. 0) 
            READ(336, *, END=100, ERR=100,IOSTAT=ISTAT) xm, xn, ds, du, dv
            m = INT(xm);  n = INT(xn)
            IF (m.GT.mpol .OR. ABS(n).GT.ntor) CYCLE
            A_s(m,n,ns) = ds/orthonorm(m,n)
            A_u(m,n,ns) = du/orthonorm(m,n)
            A_v(m,n,ns) = dv/orthonorm(m,n)
         END DO
 100     CONTINUE

         CLOSE (unit=336)
         DO js = 1, ns
           t1 = (REAL(js-1)/(ns-1))**2
           A_s(:,:,js) = A_s(:,:,ns)*t1
           A_u(:,:,js) = A_u(:,:,ns)*t1
           A_v(:,:,js) = A_v(:,:,ns)*t1
         END DO

      END IF
      RETURN
      !XXX: End of test
#endif

        dphi=(2*pi)/(nfp*nzeta)
        CALL bmw_exec ("mgrid_w7x.nc","wout_w7x_mark.nc",r1_i,z1_i,dphi,&
                        ns,A_r,A_p,A_z)

!        ns=nsp !XXX: SKS - REMEMEBR TO DELETE THIS LINE

        CALL CYL2VMEC_A(A_r,A_p,A_z,cA_s,cA_u,cA_v)

        IF (ipar.EQ.0) THEN            !e_s (sin), e_u (cos), e_v (cos)
          fours=f_sin; fouru=f_cos; fourv=f_cos;
        ELSE                           !e_s (cos), e_u (sin), e_v (sin)
          fours=f_cos; fouru=f_sin; fourv=f_sin;
        END IF

        CALL tomnsp_par(cA_s,A_s(:,:,ns:ns),fours)
        CALL tomnsp_par(cA_u,A_u(:,:,ns:ns),fouru)
        CALL tomnsp_par(cA_v,A_v(:,:,ns:ns),fourv)

#if defined(NOSKS)
!     INSIDE js = ns_match < nsp (VMEC surface)
      nmatch = MIN(ns_match, nsp)
      IF (nmatch .GT. 0) THEN
      A_s(:,:,1:nmatch) = 0; A_u(:,:,1:nmatch) = 0; A_v(:,:,1:nmatch) = 0
      DO js = 1, nmatch
         A_s(:,:,js) = phipf_i(js)*(lmn(:,:,js)+lmn(:,:,js+1))/2
         IF (ipar .EQ. iasym) CYCLE
         A_u(m0, n0, js) =-torflux(js)/orthonorm(m0,n0)
         A_v(m0, n0, js) = polflux(js)/orthonorm(m0,n0)
      END DO
      
!      A_v(m0,n0,ns) = A_v(m0,n0,nmatch) * (ns / nmatch)**2          !scales like polflux
      END IF

      q1 = 1
      DO js = nmatch+1, ns
        t1 = REAL(js-nmatch)/(ns-nmatch); t1 = t1*t1
#else    
      nmatch = 0
      t1 = 1 
      DO js = 1, ns
#endif           
        p1 = REAL(js-1,dp)/(ns-1); p1 = p1*p1; q1 = 1                   !p1 ~ flux
        DO m= 0, mpol
        q1 = q1*p1           !p1**m
        DO n = -ntor, ntor
           A_s(m,n,js) = A_s(m,n,ns)*q1
           A_u(m,n,js) = A_u(m,n,ns)*q1 
           A_v(m,n,js) = A_v(m,n,ns)*q1 
           IF (nmatch .GT. 0) THEN
              A_s(m,n,js) = t1*A_s(m,n,js) + (1-t1)*A_s(m,n,nmatch)
              A_u(m,n,js) = t1*A_u(m,n,js) + (1-t1)*A_u(m,n,nmatch)
              A_v(m,n,js) = t1*A_v(m,n,js) + (1-t1)*A_v(m,n,nmatch)
           END IF
        END DO
        END DO
      END DO
       
      END SUBROUTINE INIT_A

      
      SUBROUTINE INIT_F (Fsupsmn, Fsupumn, Fsupvmn, jcurrumn,           &
                        jcurrvmn, ipar)
      USE metrics, ONLY: phipf_i, chipf_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns)  :: Fsupsmn, Fsupumn, Fsupvmn
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns)  :: jcurrumn, jcurrvmn
      INTEGER :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: js
      REAL(dp) :: t1 
!-----------------------------------------------
!NOTE: mu0 factored into j's in metric; (note: F -> -F in CURLB_PCHELMS,
!to get correct sign of hessian diagonal elements)
      
      nsmin = startglobrow;  nsmax = endglobrow
      
      CALL ASSERT(ALL(jcurrumn(:,:,nsp+1:).EQ.0._dp), 'jcurumn != 0 in vacuum!')
      CALL ASSERT(ALL(jcurrvmn(:,:,nsp+1:).EQ.0._dp), 'jcurvmn != 0 in vacuum!')
      
   
      Fsupsmn(:,:,nsmin:nsmax) = 0
      Fsupumn(:,:,nsmin:nsmax) = jcurrumn(:,:,nsmin:nsmax)
      Fsupvmn(:,:,nsmin:nsmax) = jcurrvmn(:,:,nsmin:nsmax)
      
!SMOOTH EDGE CURRENT GRADIENT      
      DO js = nsp-10, nsp
         t1 = 1 - REAL(js-1,dp)/(nsp-1)
         Fsupumn(:,:,js) = t1*Fsupumn(:,:,js)
         Fsupvmn(:,:,js) = t1*Fsupvmn(:,:,js)
      END DO

!IMPROVES CONVERGENCE BY SUPPRESSING jac*B^u(m=1) NEAR AXIS (AT LEAST FOR W7X CASE)
      IF (ns_match .GT. 0) RETURN
      Fsupumn(m1,:,1:nsp/4) = 0
      Fsupvmn(m1,:,1:nsp/4) = 0

      END SUBROUTINE INIT_F


      SUBROUTINE COMPARE_CURRENT(Fsupumn, Fsupvmn, jcurrumn, jcurrvmn)         
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:,:) :: Fsupumn, Fsupvmn,                   &
                                    jcurrumn, jcurrvmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ns_mid, m, n, moff, noff
      REAL(dp) :: t1 
!-----------------------------------------------
!NOTE: mu0 factored into j's in metric; minus sign because r = b - Ax
!Only do this once, so no need to parallelize (note: F -> -F to get correct
!sign of hessian diagonal elements)
      RETURN

      nsmin = startglobrow;  nsmax = endglobrow
      ns_mid = (nsmin+nsmax)/2

      CALL ASSERT(LBOUND(fsupumn,1) .EQ. LBOUND(jcurrumn,1),' LBOUND1 FAILED')
      CALL ASSERT(LBOUND(fsupumn,2) .EQ. LBOUND(jcurrumn,2),' LBOUND2 FAILED')

      WRITE (1000+iam, *) ' JS = ', ns_mid
      WRITE (1000+iam, *) '    M     N      J^u         F^u         J^v         F^v'
      DO n = -ntor, ntor
         noff = n + ntor + LBOUND(fsupumn,2)
         DO m = 0, mpol
            moff = m + LBOUND(fsupumn,1)
            t1 = orthonorm(m, n)
            WRITE (1000+iam, 100) m, n,                                 &
            jcurrumn(moff,noff,ns_mid)*t1, -Fsupumn(moff,noff,ns_mid)*t1,&
            jcurrvmn(moff,noff,ns_mid)*t1, -Fsupvmn(moff,noff,ns_mid)*t1
         END DO
      END DO

 100  FORMAT (2i6, 1p,4e12.4)

      END SUBROUTINE COMPARE_CURRENT


      SUBROUTINE DUMP_A(js, iunit)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: js, iunit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, moff, noff, mp, np, nb, is
!-----------------------------------------------
!     xc = xc+xc0 ON ENTRY
      
      nsmin = LBOUND(Asubsmnsf,3); nsmax = UBOUND(Asubsmnsf,3)
      moff  = LBOUND(Asubsmnsf,1); noff  = LBOUND(Asubsmnsf,2)

      IF (js.LT.nsmin .OR. js.GT.nsmax) RETURN
      nb = MIN(2, ntor)

      IF (js .EQ. 1) THEN
           WRITE(iunit+10,*)' RAD     A_s         A_u         A_v'
         DO n = -nb, nb
     np = n+ntor+noff
     DO m = 0, 1
     mp = m+moff
     IF (n.LT.0 .AND. m.EQ.0) CYCLE
     WRITE(iunit+10,*)' N=',n,' M=',m
     DO is = 1, ns
     WRITE(iunit+10,110)  is, Asubsmnsf(mp,np,is)*orthonorm(m,n),       &
                              Asubumncf(mp,np,is)*orthonorm(m,n),       &
                              Asubvmncf(mp,np,is)*orthonorm(m,n)
     END DO
     WRITE(iunit+10,*)
     END DO
     END DO
  110    FORMAT(i4,3(1p,e12.4))
      END IF


      WRITE(iunit, *) ' RADIUS(JS): ', js, ' NEQS: ', neqs  

      WRITE(iunit, '(2a)') '       M       N      A_s         ',       &
           'A_u         A_v       jcurru      jcurrv'
      DO n=-ntor, ntor
         np = n+ntor+noff
      DO m= 0, mpol
         mp = m+moff
         WRITE(iunit, 100) m, n, Asubsmnsf(mp,np,js)*orthonorm(m,n),   &
                       Asubumncf(mp,np,js)*orthonorm(m,n),         & 
               Asubvmncf(mp,np,js)*orthonorm(m,n),         &
               jcurrumnc(m,n,js)*orthonorm(m,n),           &
               jcurrvmnc(m,n,js)*orthonorm(m,n)
         IF (.NOT.lasym) CYCLE
         WRITE(iunit+100, 100) m, n, Asubsmncf(mp,np,js)*orthonorm(m,n), &
                       Asubumnsf(mp,np,js)*orthonorm(m,n),              &
               Asubvmnsf(mp,np,js)*orthonorm(m,n),                      &
               jcurrumns(m,n,js)*orthonorm(m,n),                        &
               jcurrvmns(m,n,js)*orthonorm(m,n)
      END DO
      END DO

      WRITE(iunit,*)

 100  FORMAT(2i8, 5(1p,e12.4))

      END SUBROUTINE DUMP_A


      SUBROUTINE DUMP_B(js, iunit)
      USE metrics, ONLY: r1_i, chipf_i, phipf_i
      USE hessian, ONLY: gather_array
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: js, iunit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, is, i, j, n1, n2
      REAL(dp) :: rjac, rbsupv, rho, t1
!-----------------------------------------------
!     xc = xc+xc0 ON ENTRY
!     Here, bsupX = jac*bsupX
      n1 = startglobrow;  n2 = endglobrow

      IF (js .EQ. 1) THEN
      
         ALLOCATE (jacobmnch(0:mpol,-ntor:ntor,ns))
         CALL tomnsp_par(jacobh, jacobmnch, f_cos)
         IF (lasym) THEN
            ALLOCATE (jacobmnsh(0:mpol,-ntor:ntor,ns))
            CALL tomnsp_par(jacobh, jacobmnsh, f_sin)
         END IF
 
         CALL Gather_Array(jbsupsmnsh(:,:,n1:n2), jbsupsmnsh)
         CALL Gather_Array(jbsupumnch(:,:,n1:n2), jbsupumnch)
         CALL Gather_Array(jbsupvmnch(:,:,n1:n2), jbsupvmnch)
      
         WRITE(iunit+10,*)' EDGE VALUES OF jB^X FOR ALL M,N'
         WRITE(iunit+10,'(a)')'   M   N   jB^s      jB^u      jB^v'
                           
      IF (endglobrow .EQ. NS) THEN                      
      is = ns
      DO m = 0, mpol
      DO n = -ntor,ntor
      IF (n.LT.0 .AND. m.EQ.0) CYCLE
      t1 = -orthonorm(m,n)
      WRITE(iunit+10,50)  m, n,                                         &
       (1.5_dp*jbsupsmnsh(m,n,is) - 0.5_dp*jbsupsmnsh(m,n,is-1))*t1,    &
       (1.5_dp*jbsupumnch(m,n,is) - 0.5_dp*jbsupumnch(m,n,is-1))*t1,    &
       (1.5_dp*jbsupvmnch(m,n,is) - 0.5_dp*jbsupvmnch(m,n,is-1))*t1
      END DO
      END DO
 50   FORMAT(2i4,3(1p,e10.2))
      WRITE(iunit+10,*)
      WRITE(iunit+10,'(a)')' RAD    jB^s       jB^u       jB^v' 
      DO n = -1,1
      DO m = 0, 1
      IF (n.LT.0 .AND. m.EQ.0) CYCLE
      t1 = -orthonorm(m,n)
      WRITE(iunit+10,*)' N=',n,' M=',m
      DO is = 2, ns
      WRITE(iunit+10,110)  is,                                          &
                          jbsupsmnsh(m,n,is)*t1,                        &
                          jbsupumnch(m,n,is)*t1,                        &
                          jbsupvmnch(m,n,is)*t1
      END DO
      WRITE(iunit+10,*)
      END DO
      END DO
  110 FORMAT(i4,3(1p,e11.3))
      CALL FLUSH(iunit+10)
      END IF   !endglobrow=ns

      RETURN
      END IF   !js=1

      IF (iam .EQ. 0) THEN
      WRITE(iunit, *) 'HALF RADIUS (JS): ', js, ' NEQS: ', neqs  

      WRITE(iunit, '(2a)') '       M       N   jB^s(sin)   ',           &
           'jB^u(cos)   jB^v(cos)  jacob(cos)'
      IF (lasym) THEN
      WRITE(iunit+100, '(2a)') '       M       N      jB^s(cos)       ',&
           'jB^u(sin)       jB^v(sin)      jacob(sin)'
      END IF

      DO n=-ntor, ntor
      DO m= 0, mpol
         t1 = -orthonorm(m,n)
         WRITE(iunit, 100) m, n, jbsupsmnsh(m,n,js)*t1,   &
                       jbsupumnch(m,n,js)*t1,             & 
                       jbsupvmnch(m,n,js)*t1,             &
                       jacobmnch(m,n,js)*t1
         IF (.NOT.lasym) CYCLE
         WRITE(iunit+100, 100) m, n, jbsupsmnch(m,n,js)*t1, &
                       jbsupumnsh(m,n,js)*t1,         &
                       jbsupvmnsh(m,n,js)*t1,         &
                       jacobmnsh(m,n,js)*t1
      END DO
      END DO

      WRITE(iunit,*)

 100  FORMAT(2i8, 5(1p,e12.4))
      CALL FLUSH(iunit)
      END IF !iam=0

      IF (js .NE. ns) RETURN
      
!ESTIMATE IOTA PROFILE ON EXTENDED MESH (NOTE bsupX is jac * bsupX here and minus sign for sign of jacobian)
      DO is = nsp+1, ns
         chipf_i(is) = -chipf_i(is-1)-2*jbsupumnch(0,0,is)
         phipf_i(is) = -phipf_i(is-1)-2*jbsupvmnch(0,0,is)
      END DO
      
!TEST : CHECK EDGE VALUE
      IF (iam .EQ. 0) THEN

      rjac = 0; rbsupv = 0
      DO n=-ntor, ntor
      DO m= 0, mpol
         rjac = rjac + jacobmnch(m,n,ns)*orthonorm(m,n)
         rbsupv = rbsupv + jbsupvmnch(m,n,ns)*orthonorm(m,n)
      END DO
      END DO

      rbsupv = rbsupv*r1_i(1,1,ns)/rjac

      PRINT 150,' jac(u=0,v=0,ns): ', rjac,' jacobh: ', jacobh(1,1,ns)
 150  FORMAT(2(a,1p,e12.3))
      PRINT 155,' R*B^v(u=0,v=0,ns) : ', rbsupv
 155  FORMAT(a,1p,e12.3)

!ESTIMATE IOTA PROFILE (NOTE bsupX is jac * bsupX here and minus sign for sign of jacobian)
      PRINT *
      PRINT *,'  RAD     CHIP(WOUT)     CHIP(PCH)    PHIP(WOUT)     PHIP(PCH)'
      DO is = 2, ns, 10
          rho = (is-1.5_dp)*hs
          PRINT 200,rho, 0.5_dp*(chipf_i(is)+chipf_i(is-1)),            &
            -jbsupumnch(0,0,is),                                        &
            0.5_dp*(phipf_i(is)+phipf_i(is-1)), -jbsupvmnch(0,0,is)
      END DO
 200  FORMAT(f7.3,1p,4e14.3)
 
      END IF   !iam = 0
!END TEST

      END SUBROUTINE DUMP_B

      SUBROUTINE GETFSQ_PCHELMS
      USE shared_data, ONLY: fsupsmnsf, fsupumncf, fsupvmncf,           &
                             fsupsmncf, fsupumnsf, fsupvmnsf

#if defined(MPI_OPT)
      USE mpi_inc
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: f1
!-----------------------------------------------
      nsmin = startglobrow;  nsmax = endglobrow
      
      f1 = SUM(Fsupsmnsf(:,:,nsmin:nsmax)**2 +                          &
               Fsupumncf(:,:,nsmin:nsmax)**2 +                          &
               Fsupvmncf(:,:,nsmin:nsmax)**2)
      IF (lasym) f1 = f1 +                                              &
           SUM(Fsupsmncf(:,:,nsmin:nsmax)**2 +                          &
               Fsupumnsf(:,:,nsmin:nsmax)**2 +                          &
               Fsupvmnsf(:,:,nsmin:nsmax)**2)
#if defined(MPI_OPT)
      CALL MPI_ALLREDUCE(f1, fsq_total, 1, MPI_REAL8, MPI_SUM,          & 
                         SIESTA_COMM, MPI_ERR)
#else
      fsq_total = f1
#endif
      IF (bnorm .GT. zero) fsq_total = fsq_total/bnorm
      
      END SUBROUTINE GETFSQ_PCHELMS


      SUBROUTINE CLEANUP

      DEALLOCATE(xc0, gc0)
      DEALLOCATE(jacobmnch)

      IF (lasym) THEN
         DEALLOCATE(jacobmnsh)
      END IF

      DEALLOCATE(Bsupsijh, Bsupuijh, Bsupvijh)
      DEALLOCATE(Bsubsijh, Bsubuijh, Bsubvijh)
      
      xc = 0;  gc = 0

      END SUBROUTINE CLEANUP

      
      SUBROUTINE CHECK_CURRENT
      USE metrics, ONLY: jsupvdotA
!Check that the line integral of Bsubuij equals the integrated toroidal current
!-----------------------------------------------
      INTEGER  :: i
      REAL(dp) :: currv_int
!-----------------------------------------------

!Average over toroidal angle too (should not dependent on toroidal angle)
      IF (endglobrow .NE. ns) RETURN
      
      currv_int = hs*(SUM(jcurrvmnc(0,0,1:ns-1)) + jcurrvmnc(0,0,ns)/2)

      PRINT 100,' Line Integral B_u: ', line_bu,                        &
                ' Surface-integrated current: ', currv_int,             &
                ' VMEC Surface-integrated current: ', jsupvdotA

 100  FORMAT(a,1p,e12.4,4x,a,1pe12.4)

      END SUBROUTINE CHECK_CURRENT

      SUBROUTINE COMPUTE_FORCES_LIN

      lHessian = .TRUE.
      CALL COMPUTE_FORCES(.FALSE.)
      lHessian = .FALSE.

      END SUBROUTINE COMPUTE_FORCES_LIN

      SUBROUTINE BACKSLV
      USE metrics, ONLY: r1_i
      USE hessian, ONLY: gather_array
      USE vmec_info, ONLY: lmns=>lmns_i, lmnc=>lmnc_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, js
!-----------------------------------------------
!update solution (already gathered to all processors)
      xc = xc + xc0

      CALL INIT_F(Fsupsmnsf, Fsupumncf, Fsupvmncf,                      &
                  jcurrumnc, jcurrvmnc, isym)
      IF (lasym)              &
      CALL INIT_F(Fsupsmncf, Fsupumnsf, Fsupvmnsf,                      &
                  jcurrumns, jcurrvmns, iasym)
      gc0 = gc

      CALL COMPUTE_FORCES(.TRUE.)
      CALL gather_array(gc)

      IF (iam .EQ. 0) PRINT 200,' Back Solve Check', fsq_total,         &
                      MAXVAL(ABS(gc)), SQRT(SUM(xc*xc))

 200  FORMAT(/,a,/,' |F|^2: ',1p,e10.3,' MAX |F|: ', 1pe10.3, ' |X|: ',1pe10.3)      

      END SUBROUTINE BACKSLV

      SUBROUTINE init_pchelms
      USE blocktridiagonalsolver_s, ONLY: Initialize
      USE descriptor_mod, ONLY: iam, nprocs
      USE nscalingtools, ONLY: initRemap
      USE shared_functions, ONLY: init_ptrs
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n1, istat
!-----------------------------------------------
      IF (lasym) THEN
         ndims = 6
      ELSE
         ndims = 3
      END IF

      mblk_size=ndims*mnmax
      
#if defined(SKS)
      CALL Initialize(ns, mblk_size) 
#else
          startglobrow = 1; endglobrow = ns
#endif

          n1 = ns*mnmax
          neqs = ndims*n1

!          ALLOCATE(col_scale(0:mpol,-ntor:ntor,ndims,ns))
!          col_scale = 1

!          ALLOCATE(xc(neqs), stat=istat)
!          CALL ASSERT(istat.eq.0,'Allocate xc failed!')
!          xc = 0
          CALL init_ptrs (xc, asubsmnsf, asubumncf, asubvmncf,          &
                          asubsmncf, asubumnsf, asubvmnsf)

!          ALLOCATE(gc(neqs), stat=istat)
!          CALL ASSERT(istat.eq.0,'Allocate gc failed!')
!          gc = 0
          CALL init_ptrs (gc, fsupsmnsf, fsupumncf, fsupvmncf,          &
                          fsupsmncf, fsupumnsf, fsupvmnsf)

!Initialize forces
#if defined(SKS)
          CALL initRemap(mpol, ntor, ns, nprocs, iam)
#endif

          CALL Compute_Forces(.TRUE.)
          xc = 0

      END SUBROUTINE init_pchelms

!>     \brief Compute Fourier components of curl curl (A) - mu0J operator
!>      on full radial grid
      SUBROUTINE COMPUTE_FORCES (lNLinear)
      USE vmec_info, ONLY:   lmns=>lmns_i, lmnc=>lmnc_i
      USE shared_data, ONLY: asubsmnsf, asubumncf, asubvmncf,           &
                             fsupsmnsf, fsupumncf, fsupvmncf,           &
                             asubsmncf, asubumnsf, asubvmnsf,           &
                             fsupsmncf, fsupumnsf, fsupvmnsf, neqs
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN) :: lNLinear
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat, moff, noff
!-----------------------------------------------
      ! XXX : Which are nsmin and nsmax bounds not applied to all the
      ! arrays below?
      
      nsmin = startglobrow; nsmax = MIN(ns,endglobrow+1)
      ns_match = 0  !nsp /4

      IF (bnorm .LT. 0) THEN
         linit = .TRUE.
         IF (.NOT. ALLOCATED(xc0)) THEN
         ALLOCATE(xc0(neqs), gc0(neqs))
         CALL ASSERT(ALLOCATED(jbsupsmnsh),'B^S not allocated!')
         IF (lasym)                                                     &  
         CALL ASSERT(ALLOCATED(jbsupsmnch),'B^S not allocated!')
         ALLOCATE(Bsupsijh(ntheta,nzeta,nsmin:nsmax),                   &
                  Bsupuijh(ntheta,nzeta,nsmin:nsmax),                   &
                  Bsupvijh(ntheta,nzeta,nsmin:nsmax),                   &
                  Bsubsijh(ntheta,nzeta,nsmin:nsmax),                   &
                  Bsubuijh(ntheta,nzeta,nsmin:nsmax),                   &
                  Bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
         CALL ASSERT(istat.EQ.0,'Allocation failed in COMPUTE_FORCES')
         END IF
    
         CALL INIT_F (Fsupsmnsf, Fsupumncf, Fsupvmncf,                  &
                      jcurrumnc, jcurrvmnc, isym)
         IF (lasym)                                                     &
         CALL INIT_F (Fsupsmncf, Fsupumnsf, Fsupvmnsf,                  &
                      jcurrumns, jcurrvmns, iasym)
         gc0 = gc
         
         CALL GetFSQ_PCHELMS
         bnorm = fsq_total
         
         CALL INIT_A(asubsmnsf, asubumncf, asubvmncf, lmns, isym)
         IF (lasym)                                                     &
         CALL INIT_A(asubsmncf, asubumnsf, asubvmnsf, lmnc, iasym)
         xc0 = xc
      END IF

      gc = 0

!     1. Get sqrt(g)*CURLA dot grad(s,u,v) [sqrt(g) B^i] on half grid
!        in real space, for input AsubXmn on full grid

      nsmin = MAX(1, startglobrow-1); nsmax = MIN(ns, endglobrow+1)
      CALL CURLA_PCHELMS(Asubsmnsf,  Asubumncf,  Asubvmncf,             &
                        jbsupsmnsh, jbsupumnch, jbsupvmnch, isym)
                        
      IF (lasym)                                                        &                        
      CALL CURLA_PCHELMS(Asubsmncf,  Asubumnsf,  Asubvmnsf,             &
                         jbsupsmnch, jbsupumnsh, jbsupvmnsh, iasym)

!Add 1/jac factor
      nsmin = startglobrow

      Bsupsijh = Bsupsijh/jacobh(:,:,nsmin:nsmax)
      Bsupuijh = Bsupuijh/jacobh(:,:,nsmin:nsmax)
      Bsupvijh = Bsupvijh/jacobh(:,:,nsmin:nsmax)
      
!     2. Convert to covariant component of B_i on half grid
      CALL tolowerh(Bsupsijh, Bsupuijh, Bsupvijh,                       &
                    Bsubsijh, Bsubuijh, Bsubvijh, nsmin, nsmax)
                    

!     3. Get -sqrt(g)*CURLB dot grad(s,u,v) Fourier coefficients 
!        [sqrt(g)*J^i] on full grid and store in FsupXmnf  (sign is flipped
!        to make positive contributions to Hessian)
      
      nsmin = startglobrow;  nsmax = MIN(ns, endglobrow+1)
      CALL CURLB_PCHELMS(Fsupsmnsf, Fsupumncf, Fsupvmncf, Asubsmnsf,    &
                         jbsupsmnsh, isym)
      IF (lasym)                                                        &
      CALL CURLB_PCHELMS(Fsupsmncf, Fsupumnsf, Fsupvmnsf, Asubsmncf,    &
                         jbsupsmnch, iasym)
      nsmax = endglobrow
      
!DUMP DEBUG INFO      
      IF (linit) THEN
         CALL Compare_Current(Fsupumncf, Fsupvmncf, jcurrumnc, jcurrvmnc)
         IF (lasym)  &
         CALL Compare_Current(Fsupumnsf, Fsupvmnsf, jcurrumns, jcurrvmns)
      END IF

      IF (linit .OR. lNlinear) gc = gc + gc0
           
!Compute ONLY -Ax, NOT residue= -Ax + b (for preconditioner and GMRES)
!     4. Subtract [sqrt(g)*J]mn (input from VMEC) to get net force (residual=b-Ax)
      CALL BoundaryConditions(Fsupsmnsf, Fsupumncf, Fsupvmncf)
      IF (lasym)                                                        &
      CALL BoundaryConditions(Fsupsmncf, Fsupumnsf, Fsupvmnsf)

      
      
!     5. Subtract A*x0 since A_u, A_v are fixed at the edge
!        b' = b + A*x0 (new right-hand side)

      IF (linit) gc0 = gc
      
      IF (lNLinear .OR. lInit) THEN 
         linit = .FALSE.
         CALL GETFSQ_PCHELMS
      END IF

      END SUBROUTINE COMPUTE_FORCES
      
      
      SUBROUTINE BoundaryConditions(Fsupsmn, Fsupumn, Fsupvmn)
      USE hessian, ONLY: l_Compute_Hessian
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns) :: Fsupsmn, Fsupumn, Fsupvmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nmatch
!-----------------------------------------------
      IF (endglobrow .EQ. ns) THEN
        IF (.NOT.l_AsEdge) fsupsmn(:,:,ns) = 0    !Nat bc -> this can vary on s=s_edge surface
        fsupumn(:,:,ns) = 0; fsupvmn(:,:,ns) = 0
      END IF
      
     
!PRESERVE SURFACES FOR js < ns_match
      IF (ns_match .GE. nsmin) THEN
         nmatch = MIN(ns_match, nsmax)
         fsupumn(m1,:,nsmin:nmatch) = 0
         fsupvmn(m1,:,nsmin:nmatch) = 0
      END IF
      
      END SUBROUTINE BoundaryConditions
      
      
      SUBROUTINE GMRES_PCHELMS
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE gmres_lib, ONLY: gmres_info, gmres_par, gmres_ser
      USE hessian, ONLY: Apply_Precond
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE(gmres_info)   :: gi
      INTEGER, PARAMETER   :: noPrec=0, leftPrec=1, rightPrec=2, dblePrec=3
      INTEGER :: n, m, js
      REAL(dp), PARAMETER :: one = 1, ftol = 1.E-18_dp
	  REAL(dp), ALLOCATABLE :: b(:), x0(:)
!-----------------------------------------------
!
!     GMRES SOLVES r == b + Ax = 0  (it appears the sign is +)
!
      ALLOCATE(b(neqs), x0(neqs))

!      CALL TriDiag_Test (xc, gc)

      xc = 0

      CALL Compute_Forces(.TRUE.)
      IF (iam .EQ. 0) THEN
         PRINT *
         PRINT 100,' ITER: ', 1, '|F|^2: ', fsq_total, ' |X| = ', SQRT(SUM(xc0**2))
         WRITE (33, 100)' NEQS : ', neqs
         WRITE (33, 100)' ITER: ',1, '|F|^2: ', fsq_total, ' |X| = ', SQRT(SUM(xc0**2))
      END IF

         b = gc0
!USE THIS FOR FULL PRECONDITIONER ONLY
!	  x0 = -gc0
!	  CALL Apply_Precond(x0)
!OTHERWISE USE THIS
         x0 = xc
         x0 = 0
 
 100  FORMAT(a,i6,2(3x,a,1pE12.3))

!---------------
      CALL init_dgmres(gi%icntl,gi%cntl)

!*************************
!*  Tune some parameters
!*************************
! Tolerance
       gi%cntl(1) = SQRT(ftol)
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      gi%icntl(2) = 21
! Save the convergence history in file fort.20
      gi%icntl(3) = 20
! Preconditioning INPUT flags (note: different than revcom flags in driver)
!     gi%icntl(4) = leftPrec
      gi%icntl(4) = rightPrec 
!     gi%icntl(4) = dblePrec
!      gi%icntl(4) = noPrec
!! ICGS orthogonalization
      gi%icntl(5) = 3
! Initial guess (use it if = 1, otherwise ignore it)
      gi%icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
      gi%icntl(7) = 20

      gi%icntl(8) = 1             !Default
      gi%icntl(9) = 1             !Steps for peek at progress during rev com loop
      gi%l_nonlinear = .FALSE.    !linear solver (output on screen)

!*********************************
! Choose the restart parameter
!*********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      n   = neqs
      gi%m = 200
!---------------
      gi%iam = iam; gi%nprocs=nprocs; gi%ngmres_type=1
      gi%ftol = ftol
      IF (PARSOLVER) THEN
         gi%endglobrow=endglobrow; gi%startglobrow=startglobrow
         gi%mblk_size=mblk_size; gi%rcounts=>rcounts; gi%disp=>disp

         gi%my_comm = SIESTA_COMM
         gi%my_comm_world = SIESTA_COMM
!         gi%lactive = lactive

#if defined(SKS)
         CALL gmres_par (n, gi, matvec_par_PCHELMS, apply_precond,      &
                         GetNLForce_PCHELMS, x0, b)
#endif                         
      ELSE
         CALL gmres_ser (n, gi, matvec_PCHELMS, apply_precond,          &
                         GetNLForce_PCHELMS, x0, b)
      END IF
      
      xc = x0
      
      DEALLOCATE(x0, b)

      END SUBROUTINE GMRES_PCHELMS


      SUBROUTINE MatVec_PCHELMS (p, Ap, ndim)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: ndim
      REAL(dp), INTENT(IN)  :: p(ndim)
      REAL(dp), INTENT(OUT) :: Ap(ndim)
!-----------------------------------------------
      CALL ASSERT(neqs.EQ.ndim,' neqs != ndim')

      xc = p
      CALL Compute_Forces(.FALSE.)
      Ap = -gc

      END SUBROUTINE MatVec_PCHELMS

      SUBROUTINE MatVec_par_PCHELMS (ploc, Ap, nloc)
      USE hessian, ONLY: eps_factor
      USE blocktridiagonalsolver_s, ONLY: ParMatVec, PadSides 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                    :: nloc
      REAL(dp), INTENT(IN), DIMENSION(nloc)  :: ploc
      REAL(dp), INTENT(OUT), DIMENSION(nloc) :: Ap
!-----------------------------------------------
      LOGICAL, PARAMETER    :: LPARMATVEC=.FALSE.
      INTEGER               :: myrowstart, myrowend, istat
      REAL(dp), ALLOCATABLE :: p(:)
!-----------------------------------------------
!SPH NOTE 032713: EVENTUALLY REMOVE CALL TO ParMatVec IF THE FUNCT_ISLAND IS FASTER
!
!     NOTE: DO NOT CALL ParMatVec WHEN nprecon_type==PREDIAG
!           SINCE IN HESSIAN, WE SET ALL THE OFF-DIAGONAL COMPONENTS OF 
!           THE A-MATRIX  (USED IN ParMatVec) TO ZERO
!
!      CALL ASSERT(.NOT.l_init_state,'l_init_state = T in matvec_par')
      istat = (endglobrow-startglobrow+1)*mblk_size
      CALL ASSERT(istat.EQ.nloc, 'nloc wrong in matvec_par')

      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1

      ALLOCATE (p(ns*mblk_size), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in matvec_par')

      p(myrowstart:myrowend) = ploc
      CALL PadSides(p, mblk_size, 1, 1)


      IF (.NOT.LPARMATVEC) THEN
         xc = p
         CALL Compute_Forces(.FALSE.)
         Ap = -gc(myrowstart:myrowend)

!      ELSE 
!         CALL ParMatVec(p,Ap,nloc)
      END IF

      DEALLOCATE (p)

      END SUBROUTINE MatVec_par_PCHELMS

      SUBROUTINE GetNLForce_PCHELMS(xcstate, fsq_nl, bnorm)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------

	  PRINT *,' SHOULD NOT GET HERE'
	  fsq_nl = 1
	  RETURN

	  END SUBROUTINE GetNLForce_PCHELMS


      SUBROUTINE Conjugate_Gradient_PCHELMS
      USE descriptor_mod, ONLY: iam
      USE hessian, ONLY: Apply_Precond, gather_array
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: one = 1, ftol = 1.E-16_dp
      INTEGER :: niter,  istat, js, moff, noff
      REAL(dp), ALLOCATABLE :: xcdot(:), xsave(:), gsave(:), zsave(:)
	  REAL(dp) :: beta, alpha, fsq, fsq1, fsq0, xsq
	  LOGICAL :: lRestart=.TRUE.
!-----------------------------------------------
!     COMPUTE A*xc + b' = 0, WHERE b' = b+A*xc0
!     FINAL SOLUTION IS xc = xc+xc0

      ALLOCATE (xcdot(neqs), xsave(neqs), gsave(neqs),   &
	            zsave(neqs), stat=istat)
	  xcdot = 0
	  CALL ASSERT(istat.EQ.0,'Conj Grad Allocation error')

      beta = 0

!      CALL Check_Symmetry (xc, gc)
!      CALL TriDiag_Test (xc, gc)

      xc = 0

 122  CONTINUE

      CALL Compute_Forces(.TRUE.)
  	  fsq   = fsq_total
	  IF (iam.EQ.0 .AND. lRestart) THEN
	     PRINT 100,' ITER: ', 1, '|F|^2: ', fsq
	     WRITE (33, 100)' NEQS : ', neqs
         WRITE (33, 100)' ITER: ',1, '|F|^2: ', fsq
	  END IF

 100  FORMAT(a,i6,2(3x,a,1pE12.3))


!Initialize zsave = preconditioned forces      
	  zsave = gc
      CALL Apply_Precond(zsave)
      fsq = -SUM(gc*zsave)

      DO niter = 2, 100

         xcdot = zsave + beta*xcdot
         xsave = xc
		 gsave = gc
         
		 xc = xcdot
!	     Compute and store gc=A*xcdot
   	     CALL Compute_Forces(.FALSE.)
		 
         alpha = fsq/SUM(xcdot*gc)

		 xc = xsave + alpha*xcdot
		 gc = gsave + alpha*gc

         zsave = gc
	     CALL Apply_Precond(zsave)
         fsq1 = -SUM(gc*zsave)
		 beta = fsq1/fsq
		 fsq  = fsq1
		 IF (beta .GT. 1._dp) THEN
		    xc = xsave
			IF (.NOT.lRestart) EXIT
			lRestart=.FALSE.
			GOTO 122
		 END IF

	     IF (iam .EQ. 0) THEN
		    xsq = SQRT(SUM(xc*xc))
            fsq0 = SUM(gc*gc)
		    PRINT 100,' ITER: ', niter, '|F|^2: ', fsq0, ' |x|: ', xsq
		    WRITE(33,100)' ITER: ', niter, '|F|^2: ', fsq0, ' |x|: ', xsq
		 END IF
		 IF (fsq .LE. ftol) EXIT

	  END DO

      CALL gather_array(xc)

      DEALLOCATE (xcdot, xsave, zsave, gsave, stat=istat)

	  END SUBROUTINE Conjugate_Gradient_PCHELMS

    END MODULE pchelms
