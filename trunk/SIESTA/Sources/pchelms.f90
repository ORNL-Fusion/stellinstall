!*******************************************************************************
!>  @file pchelms.f90
!>  @brief Contains module @ref pchelms.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file solves the helmholtz equation to set inital fields that match vmec
!>  and vacuum currents from the vector potential. Initial vector potential on
!>  the edge is supplied by the BMW code.
!*******************************************************************************
      MODULE pchelms
      USE stel_kinds
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i
      USE island_params, ntheta=>nu_i, nzeta=>nv_i, hs=>hs_i
      USE island_params, ohs => ohs_i, mnmax=>mnmax_i, nfp=>nfp_i
      USE hessian, ONLY: initHess, Compute_Hessian_Blocks, asym_index, levmarq_param
      USE blocktridiagonalsolver_s, ONLY: Initialize, GetRanks
      USE nscalingtools, ONLY: MyEnvVariables, PARSOLVER, disp, startglobrow, endglobrow, rcounts, MPI_ERR
      USE shared_data, ONLY: gc, xc, neqs, gc0, xc0, fsq_total, mblk_size, ndims
      USE shared_data, ONLY: asubsmnsf, asubumncf, asubvmncf,                  &
                             fsupsmnsf, fsupumncf, fsupvmncf,                  &
                             asubsmncf, asubumnsf, asubvmnsf,                  &
                             fsupsmncf, fsupumnsf, fsupvmnsf
      USE shared_data, ONLY: torflux, polflux, nsp, r1_i, z1_i, ru_i,          &
                             zu_i, rv_i, zv_i, lcolscale
      USE descriptor_mod, ONLY: iam, nprocs, lscalapack, INHESSIAN, SIESTA_COMM
      USE v3_utilities, ONLY: assert

      USE vmec_info, ONLY: jcurrumnc, jcurrvmnc, jcurrumns, jcurrvmns
      USE fourier
      USE metrics, ONLY: tolowerh
      USE quantities, ONLY: jacobh, jbsupsmnsh, jbsupumnch, jbsupvmnch,        &
                            jbsupsmnch, jbsupumnsh, jbsupvmnsh

      USE mpi_inc

      IMPLICIT NONE

!*******************************************************************************
!  pchelms module parameters
!*******************************************************************************
!>  Controls if the edge Aubs values are evolved or fixed.
!>    * True  Evolve the edge Asubs values
!>    * False Fix edge Asubs values.
      LOGICAL, PARAMETER, PRIVATE :: l_AsEdge = .TRUE.

!  FIXME Remove these. They should either be passed in as arguments.
      INTEGER, PRIVATE            :: nsmin
      INTEGER, PRIVATE            :: nsmax
      INTEGER, PRIVATE            :: ns_match

    !Added SKS
      LOGICAL, PRIVATE            :: linit
      LOGICAL, PRIVATE            :: lHessian
      REAL (dp), PRIVATE          :: bnorm = -1
      REAL (dp), PRIVATE          :: line_Bu

      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsubsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsubuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: Bsubvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: jacobmnch
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: jacobmnsh

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Run the pchelms code to solve for inital jbsup values.
!-------------------------------------------------------------------------------
      SUBROUTINE run_pchelms
      USE shared_data, ONLY: lverbose
      USE diagnostics_mod, ONLY: divb, divb_rms

      IMPLICIT NONE

!  local variables
      INTEGER                     :: istat
      INTEGER                     :: js
      REAL (dp)                   :: ton
      REAL (dp)                   :: toff
      LOGICAL                     :: flag
      LOGICAL                     :: lcolscale_save

!  local parameters
      REAL (dp), PARAMETER        :: levmarq_param_inital = 1.E-6_dp

!  Start of executable code
!  Turn of column scaling. Save the value of the column scaling flag so it can
!  be reset afterwards.
      lcolscale_save = lcolscale
      lcolscale = .false.

      IF (lverbose .and. iam .EQ. 0) THEN
         WRITE (*,1000)
      END IF

      CALL init_pchelms

      CALL second0(ton)

!  Add initial diffusive term to preconditioner instead
      levmarq_param = levmarq_param_inital

!  Note initHess must have already been called prior to this point.
      CALL Compute_Hessian_Blocks(Compute_Forces_Lin, .FALSE.)
      CALL second0(toff)
      IF (lverbose .and. iam .EQ. 0) THEN
         WRITE (*,1001) asym_index, toff - ton
      END IF

      CALL second0(ton)
!  Call conjugate gradient pchelms
      CALL GMRES_PCHELMS
      CALL second0(toff)

      IF (lverbose .and. iam .EQ. 0) THEN
         WRITE (*,1002) toff - ton
      END IF

!  Update solution (initial guess + perturbation) and check final forces
      CALL BACKSLV

!  Dump output: MUST call backslv to get updated B's before dumping
      DO js = 1, ns
         CALL DUMP_A(js, 333)
         CALL DUMP_B(js, 666)
      END DO

      nsmin = startglobrow
      nsmax = MIN(ns, endglobrow + 1)
      CALL divb(nsmin, nsmax)
    
      IF (lverbose .and. iam .EQ. 0) THEN
         WRITE (*,1003) divb_rms
      END IF

      CALL CHECK_CURRENT

      CALL CLEANUP

      lcolscale = lcolscale_save

1000  FORMAT(/,'-----------------------------------------------',              &
             /,'STARTING EXTERNAL B-FIELD CALCULATION (PCHELMS)',              &
             /,'-----------------------------------------------')
1001  FORMAT(' ASYMMETRY INDEX: ',1p,e12.4,/,                                  &
             ' PCHELMS HESSIAN COMPUTATION TIME: ',f10.1,'s')
1002  FORMAT(' PCHELMS HESSIAN SOLVER TIME: ',f10.1,'s')
1003  FORMAT(/,' |DivB| (normed) = ',1p,e12.3,/,                               &
             /,'-----------------------------------------------',              &
             /,'ENDING PCHELMS',                                               &
             /,'-----------------------------------------------')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute real-space components of contravariant jac*B^i on half radial
!>  grid from curl(A).
!>
!>  @param[in]    Asubsmnf  Covariant vector potential in the s direction.
!>  @param[in]    Asubumnf  Covariant vector potential in the u direction.
!>  @param[in]    Asubvmnf  Covariant vector potential in the v direction.
!>  @param[inout] jbsupsmnh Contravariant B field and jacobian in the s
!>                          direction.
!>  @param[inout] jbsupumnh Contravariant B field and jacobian in the u
!>                          direction.
!>  @param[inout] jbsupvmnh Contravariant B field and jacobian in the v
!>                          direction.
!>  @param[in]    iparity   Parity flag to indicate stellaratory symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE CURLA_PCHELMS(Asubsmnf,  Asubumnf,  Asubvmnf,                 &
                               jbsupsmnh, jbsupumnh, jbsupvmnh, iparity)
      USE utilities, ONLY: CURL_FtoH, set_bndy_fouier_m0, set_bndy_full_origin

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: Asubsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: Asubumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: Asubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupvmnh
      INTEGER, INTENT(in)                                     :: iparity

!  local variables
      INTEGER                                                 :: fcomb
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv
      INTEGER                                                 :: moff
      INTEGER                                                 :: noff
      INTEGER                                                 :: nmin
      INTEGER                                                 :: nmatch

      INTEGER :: m,n

!  Start of executable code
      moff = LBOUND(asubsmnf, 1)
      noff = ntor + LBOUND(asubsmnf, 2)
      nmin = nsmin

      IF (iparity .EQ. f_sin) THEN
         fcomb = f_none
         fours = f_sin
         fouruv = f_cos
      ELSE
         fcomb = f_sum
         fours = f_cos
         fouruv = f_sin
      END IF

      CALL set_bndy_fouier_m0(asubsmnf, asubumnf, asubvmnf, iparity)
      IF (nsmin .EQ. 1) THEN
         CALL set_bndy_full_origin(asubsmnf, asubumnf, asubvmnf)
      END IF

!  Dirichelet bc (fixed) for A_u, A_v (A_s floats) consistent with Poynting
!  flux = 0.
!      IF (nsmax .eq. ns .and.                                                  &
!          lHessian     .and.                                                   &
!          .not.lInit) THEN
!         IF (.not.l_AsEdge) THEN
!            asubsmnf(:,:,ns) = 0
!         END IF
!         asubumnf(:,:,ns) = 0
!         asubvmnf(:,:,ns) = 0
!      END IF

!  Need this to keep Hessian symmetric
      IF (nmin .le. ns_match .and. lHessian .and. .not.lInit) THEN
          nmatch = MIN(ns_match, nsmax)
          asubumnf(m1 + moff,:,nmin:nmatch) = 0
          asubvmnf(m1 + moff,:,nmin:nmatch) = 0
      END IF

      CALL CURL_FtoH(asubsmnf,  asubumnf,  asubvmnf,                           &
                     jbsupsmnh, jbsupumnh, jbsupvmnh, iparity, nmin, nsmax)

      nmin = startglobrow

!  Convert to real space.
      CALL toijsp(jbsupsmnh(:,:,nmin:nsmax), Bsupsijh, fcomb, fours)
      CALL toijsp(jbsupumnh(:,:,nmin:nsmax), Bsupuijh, fcomb, fouruv)
      CALL toijsp(jbsupvmnh(:,:,nmin:nsmax), Bsupvijh, fcomb, fouruv)

      END SUBROUTINE
      
!-------------------------------------------------------------------------------
!>  @brief Compute real-space components of contravariant jac*J^i on full radial
!>  grid from curl(B).
!>
!>  @param[inout] ksupsmnf Contravariant current and jacobian in the s direction
!>                         on the full grid.
!>  @param[inout] ksupumnf Contravariant current and jacobian in the u direction
!>                         on the full grid.
!>  @param[inout] ksupvmnf Contravariant current and jacobian in the v direction
!>                         on the full grid.
!>  @param[inout] asubsmnf Covariant vector potential in the s direction on the
!>                         full grid.
!>  @param[in]    bsupsmnh Contravariant current and jacobian in the v direction
!>                         on the half grid.
!>  @param[in]    iparity  Parity flag to indicate stellaratory symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE CURLB_PCHELMS(ksupsmnf, ksupumnf, ksupvmnf, asubsmnf,         &
                               bsupsmnh, iparity)
      USE utilities, ONLY: CURL_HtoF, GradientHalf, set_bndy_fouier_m0
!  Declare Arguments
      REAL(dp), POINTER, DIMENSION(:,:,:)                 :: ksupsmnf
      REAL(dp), POINTER, DIMENSION(:,:,:)                 :: ksupumnf
      REAL(dp), POINTER, DIMENSION(:,:,:)                 :: ksupvmnf
      REAL(dp), POINTER, DIMENSION(:,:,:)                 :: asubsmnf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN) :: bsupsmnh
      INTEGER, INTENT(IN)                                 :: iparity

!  Local variables
      INTEGER                                             :: fours
      INTEGER                                             :: fouruv
      INTEGER                                             :: moff
      INTEGER                                             :: noff
      INTEGER                                             :: nl
      INTEGER                                             :: nh
      INTEGER                                             :: istat
      INTEGER                                             :: nmin
      INTEGER                                             :: nmax
      INTEGER                                             :: mp
      INTEGER                                             :: np
      INTEGER                                             :: m
      INTEGER                                             :: n
      INTEGER                                             :: sparity
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE            :: dA_sds
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE            :: bsubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE            :: bsubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE            :: bsubvmnh
      REAL (dp)                                           :: Diff
      REAL (dp), DIMENSION(nsp)                           :: t1(nsp)

!  Start of executable code
      IF (iparity .eq. f_sin) THEN
         fours = f_sin
         fouruv = f_cos
         sparity = 1
      ELSE
         fours = f_cos
         fouruv = f_sin
         sparity = -1
      END IF

      nmin = startglobrow
      nmax = MIN(ns, endglobrow+1)
      moff = LBOUND(ksupsmnf,1)
      noff = LBOUND(ksupsmnf,2) + ntor

      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,nmin:nmax),                          &
               bsubumnh(0:mpol,-ntor:ntor,nmin:nmax),                          &
               bsubvmnh(0:mpol,-ntor:ntor,nmin:nmax), stat=istat)   
      CALL ASSERT(istat .eq. 0,'Allocation1 failed in CURLB_PCHELMS')

      CALL tomnsp(bsubsijh, bsubsmnh(:,:,nmin:nmax), fours)          ! Compute BSUBXMNH
      CALL tomnsp(bsubuijh, bsubumnh(:,:,nmin:nmax), fouruv)
      CALL tomnsp(bsubvijh, bsubvmnh(:,:,nmin:nmax), fouruv)

      CALL set_bndy_fouier_m0(bsubsmnh, bsubumnh, bsubvmnh, fours)

!  Sub s term has sin parity so the subu term has cos parity.
      IF (nmax .eq. ns .and. iparity .eq. f_sin) THEN
         line_Bu = bsubumnh(0,0,ns)
      END IF

      CALL CURL_HtoF(bsubsmnh, bsubumnh, bsubvmnh,                             &
                     Ksupsmnf, Ksupumnf, Ksupvmnf,                             &
                     iparity, nmin, nmax, .False.)

      nmax = endglobrow
      DEALLOCATE (bsubsmnh, bsubumnh, bsubvmnh, stat=istat)   
      CALL assert_eq(0, istat, 'Deallocation failed in CURLB_PCHELMS')
#if 0
!  Add smoothing diffusion for A_s in Hessian:
      IF (lHessian) THEN
         np = MAX(1, nmin - 1)
         mp = MIN(ns, nmax + 1)
         ALLOCATE (dA_sds(0:mpol,-ntor:ntor,np:mp))
         CALL GradientHalf(dA_sds, asubsmnf)
         Diff = 1.E-4_dp*ohs
         ksupsmnf(:,:,nmin:mp-1) = ksupsmnf(:,:,nmin:mp-1)                     &
                                 + Diff*(dA_sds(:,:,nmin+1:mp) -               &
                                         dA_sds(:,:,nmin:mp-1))
         IF (nmax .eq. ns) THEN
            ksupsmnf(:,:,ns) = ksupsmnf(:,:,ns) - 2*Diff*dA_sds(:,:,ns)
         END IF
         IF (iparity .eq. isym) THEN
            ksupsmnf(m0+moff,n0+noff,nmin:nmax) = 0
         END IF
         DEALLOCATE(dA_sds)
      END IF

      nh = MIN(nmax, nsp)
      IF (nh .GE. nmin) THEN
!Force b^s -> 0 for js <= nsp (add (B^s)**2 to energy)
         nl = nmin
         Diff = 0.1_dp
         DO m = 0,mpol
            mp = m*sparity
            DO n = -ntor,ntor
               np = n*sparity*nfp
               t1(nl:nh) = Diff*(bsupsmnh(m,n,nl:nh)                           &
                         +       bsupsmnh(m,n,nl+1:nh+1))
               ksupumnf(m+moff,n+noff,nl:nh) = ksupumnf(m+moff,n+noff,nl:nh)   &
                                             - np*t1(nl:nh)
               ksupvmnf(m+moff,n+noff,nl:nh) = ksupvmnf(m+moff,n+noff,nl:nh)   &
                                             + mp*t1(nl:nh)
            END DO
         END DO
      END IF
#endif

!  Origin boundary conditions: If asubs(0,:,1) = 0 in curla, must set
!                              ksupsmnf(0,:,1) = 0 here for hessian symmetry.
!  Need factor of 1/2 for hessian symmetry at js=1
!      IF (nmin .eq. 1) THEN
!         ksupsmnf(:,:,1) = 0
!         ksupumnf(:,:,1) = 0
!         ksupvmnf(m0+moff,:,1) = 0.5*ksupvmnf(m0+moff,:,1)
!         ksupvmnf(m1+moff,:,1) = 0
!      END IF
!      IF (nmax .eq. ns) THEN
!         ksupsmnf(:,:,ns) = 0.5*ksupsmnf(:,:,ns)
!         ksupumnf(:,:,ns) = 0.5*ksupumnf(:,:,ns)
!         ksupvmnf(:,:,ns) = 0.5*ksupvmnf(:,:,ns)
!      END IF

!  Flip sign to get correct diagonal sign (>0).
      ksupsmnf(:,:,nmin:nmax) = -ksupsmnf(:,:,nmin:nmax)
      ksupumnf(:,:,nmin:nmax) = -ksupumnf(:,:,nmin:nmax)
      ksupvmnf(:,:,nmin:nmax) = -ksupvmnf(:,:,nmin:nmax)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert cylindical vector potential to contravariant components.
!>
!>  @param[in] A_r Cylindical r component of the vector potential at the last
!>                 surface.
!>  @param[in] A_p Cylindical p component of the vector potential at the last
!>                 surface.
!>  @param[in] A_z Cylindical z component of the vector potential at the last
!>                 surface.
!>  @param[out] cA_s Contravariant vector potential in the s direction on the
!>                   last surface grid.
!>  @param[out] cA_u Contravariant vector potential in the u direction on the
!>                   last surface grid.
!>  @param[out] cA_v Contravariant vector potential in the v direction on the
!>                   last surface grid.
!-------------------------------------------------------------------------------
      SUBROUTINE CYL2VMEC_A(A_r, A_p, A_z, cA_s, cA_u, cA_v)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), INTENT(IN)  :: A_r
      REAL (dp), DIMENSION(:,:), INTENT(IN)  :: A_p
      REAL (dp), DIMENSION(:,:), INTENT(IN)  :: A_z
      REAL (dp), DIMENSION(:,:), INTENT(OUT) :: cA_s
      REAL (dp), DIMENSION(:,:), INTENT(OUT) :: cA_u
      REAL (dp), DIMENSION(:,:), INTENT(OUT) :: cA_v

!  Local variables
      REAL (dp), DIMENSION(:,:), ALLOCATABLE :: rs
      REAL (dp), DIMENSION(:,:), ALLOCATABLE :: zs
      INTEGER                                :: u
      INTEGER                                :: v
      INTEGER                                :: js

!  Start of executable code
      ALLOCATE(rs(ntheta, nzeta))
      ALLOCATE(zs(ntheta, nzeta))
      rs = (r1_i(:,:,ns) - r1_i(:,:,ns - 1))*ohs
      zs = (z1_i(:,:,ns) - z1_i(:,:,ns - 1))*ohs

      cA_s(:,:) = A_r(:,:)*rs(:,:) + A_z(:,:)*zs(:,:)
      cA_u(:,:) = A_r(:,:)*ru_i(:,:,ns) + A_z(:,:)*zu_i(:,:,ns)
      cA_v(:,:) = A_r(:,:)*rv_i(:,:,ns) + A_z(:,:)*zv_i(:,:,ns)                &
                + A_p(:,:)*r1_i(:,:,ns)

      DEALLOCATE(rs,zs)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize vector potential.
!>
!>  This subroutine is only called once so no need to parallelize. Linearly
!>  extrapolat to the center or match surface.
!>
!>  @param[inout] cA_s    Contravariant vector potential in the s direction on
!>                        the full grid.
!>  @param[inout] cA_u    Contravariant vector potential in the u direction on
!>                        the full grid.
!>  @param[inout] cA_v    Contravariant vector potential in the v direction on
!>                        the full grid.
!>  @param[out]   A_s     Contravariant fourier vector potential in the s
!>                        direction on the full grid.
!>  @param[out]   A_u     Contravariant fourier vector potential in the u
!>                        direction on the full grid.
!>  @param[out]   A_v     Contravariant fourier vector potential in the v
!>                        direction on the full grid.
!>  @param[in]    iparity Parity flag to indicate stellaratory symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE INIT_A(cA_s, cA_u, cA_v, A_s, A_u, A_v, iparity)
      USE metrics, ONLY: phipf_i, chipf_i
      USE utilities, ONLY: set_bndy_fouier_m0, set_bndy_full_origin

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: cA_s
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: cA_u
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: cA_v
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: A_s
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: A_u
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: A_v
      INTEGER, INTENT(IN)                                     :: iparity

!  Local variables
      INTEGER                                                 :: js
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: nmatch
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv
      REAL (dp)                                               :: t1
      REAL (dp)                                               :: p1
      REAL (dp)                                               :: q1

!  Start of executable code
!  Only do this once, so no need to parallelize
      A_s = 0
      A_u = 0
      A_v = 0

      IF (iparity .EQ. f_sin) THEN
!  e_s (sin), e_u (cos), e_v (cos)
         fours = f_sin
         fouruv = f_cos
      ELSE
!  e_s (cos), e_u (sin), e_v (sin)
         fours = f_cos
         fouruv = f_sin
      END IF

      CALL tomnsp(cA_s, A_s(:,:,ns:ns), fours)
      CALL tomnsp(cA_u, A_u(:,:,ns:ns), fouruv)
      CALL tomnsp(cA_v, A_v(:,:,ns:ns), fouruv)

!  Extapolate the vector potential from the edge to the center.
      nmatch = 0
      t1 = 1 
      DO js = 1, ns
         p1 = REAL(js-1,dp)/(ns-1)
         p1 = p1*p1
         q1 = 1                   !p1 ~ flux
         DO m = 0, mpol
            q1 = q1*p1           !p1**m
            DO n = -ntor, ntor
               A_s(m,n,js) = A_s(m,n,ns)*q1
               A_u(m,n,js) = A_u(m,n,ns)*q1
               A_v(m,n,js) = A_v(m,n,ns)*q1
               IF (nmatch .gt. 0) THEN
                  A_s(m,n,js) = t1*A_s(m,n,js) + (1 - t1)*A_s(m,n,nmatch)
                  A_u(m,n,js) = t1*A_u(m,n,js) + (1 - t1)*A_u(m,n,nmatch)
                  A_v(m,n,js) = t1*A_v(m,n,js) + (1 - t1)*A_v(m,n,nmatch)
               END IF
            END DO
         END DO
      END DO

      CALL set_bndy_fouier_m0(A_s, A_u, A_v, fours)
      CALL set_bndy_full_origin(A_s, A_u, A_v)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize expected currents.
!>
!>  Up to the vmec last closed flux surface, initialize with the VMEC currents.
!>  outside, initalize to zero.
!>
!>  @param[out] Fsupsmn  Contravariant current in the s direction on the full
!>                       grid.
!>  @param[out] Fsupumn  Contravariant current in the u direction on the full
!>                       grid.
!>  @param[out] Fsupvmn  Contravariant current in the v direction on the full
!>                       grid.
!>  @param[in]  jcurrumn Contravariant VMEC current in the u direction on the
!>                       full grid.
!>  @param[in]  jcurrvmn Contravariant VMEC current in the v direction on the
!>                       full grid.
!-------------------------------------------------------------------------------
      SUBROUTINE INIT_F(Fsupsmn, Fsupumn, Fsupvmn, jcurrumn, jcurrvmn)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: Fsupsmn
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: Fsupumn
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(out) :: Fsupvmn
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(in)  :: jcurrumn
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(in)  :: jcurrvmn

!  Local variables
      INTEGER                                                 :: js
      REAL (dp)                                               :: t1

!  Start of executable code
!NOTE: mu0 factored into j's in metric; (note: F -> -F in CURLB_PCHELMS,
!to get correct sign of hessian diagonal elements)
      nsmin = startglobrow
      nsmax = endglobrow
      
      CALL ASSERT(ALL(jcurrumn(:,:,nsp+1:).EQ.0._dp), 'jcurumn != 0 in vacuum!')
      CALL ASSERT(ALL(jcurrvmn(:,:,nsp+1:).EQ.0._dp), 'jcurvmn != 0 in vacuum!')

      Fsupsmn(:,:,nsmin:nsmax) = 0
      Fsupumn(:,:,nsmin:nsmax) = jcurrumn(:,:,nsmin:nsmax)
      Fsupvmn(:,:,nsmin:nsmax) = jcurrvmn(:,:,nsmin:nsmax)

!  Smooth edge current gradient.
      DO js = nsp - 10, nsp
         t1 = 1 - REAL(js-1,dp)/(nsp-1)
         Fsupumn(:,:,js) = t1*Fsupumn(:,:,js)
         Fsupvmn(:,:,js) = t1*Fsupvmn(:,:,js)
      END DO

!      Fsupumn(:,:,1) = 0
!      IF (iparity .eq. 0) THEN
!         Fsupvmn(m1:,:,1) = 0
!      ELSE
!         Fsupvmn(m1:,:,1) = 0
!      END IF

!IMPROVES CONVERGENCE BY SUPPRESSING jac*B^u(m=1) NEAR AXIS (AT LEAST FOR W7X CASE)
!      IF (ns_match .GT. 0) RETURN
!      Fsupumn(m1,:,1:nsp/4) = 0
!      Fsupvmn(m1,:,1:nsp/4) = 0

      END SUBROUTINE INIT_F

!-------------------------------------------------------------------------------
!>  @brief Compare Curl(Curl(A)) with the expected VMEC and vacuum currents.
!>
!>  @param[out] Fsupsmn  Contravariant current in the s direction on the full
!>                       grid.
!>  @param[out] Fsupumn  Contravariant current in the u direction on the full
!>                       grid.
!>  @param[out] Fsupvmn  Contravariant current in the v direction on the full
!>                       grid.
!>  @param[in]  jcurrumn Contravariant VMEC current in the u direction on the
!>                       full grid.
!>  @param[in]  jcurrvmn Contravariant VMEC current in the v direction on the
!>                       full grid.
!>  @param[in]  iparity  Parity flag to indicate stellaratory symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE COMPARE_CURRENT(Fsupumn, Fsupvmn, jcurrumn, jcurrvmn)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in) :: Fsupumn
      REAL (dp), DIMENSION(:,:,:), INTENT(in) :: Fsupvmn
      REAL (dp), DIMENSION(:,:,:), INTENT(in) :: jcurrumn
      REAL (dp), DIMENSION(:,:,:), INTENT(in) :: jcurrvmn

!  Local variables
      INTEGER                                 :: ns_mid
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: moff
      INTEGER                                 :: noff
      REAL (dp)                               :: t1

!  Start of executable code
!  Note: mu0 factored into j's in metric; minus sign because r = b - Ax
!        Only do this once, so no need to parallelize (note: F -> -F to get
!        correct sign of hessian diagonal elements)
      nsmin = startglobrow
      nsmax = endglobrow
      ns_mid = (nsmin + nsmax)/2

      CALL assert_eq(LBOUND(fsupumn,1), LBOUND(jcurrumn,1), ' LBOUND1 FAILED')
      CALL assert_eq(LBOUND(fsupumn,2), LBOUND(jcurrumn,2), ' LBOUND2 FAILED')

      WRITE (1000 + iam, 1000) ns_mid
      DO n = -ntor, ntor
         noff = n + ntor + LBOUND(fsupumn,2)
         DO m = 0, mpol
            moff = m + LBOUND(fsupumn,1)
            t1 = orthonorm(m, n)
            WRITE (1000 + iam, 1001) m, n,                                     &
                                     jcurrumn(moff,noff,ns_mid)*t1,            &
                                     -Fsupumn(moff,noff,ns_mid)*t1,            &
                                     jcurrvmn(moff,noff,ns_mid)*t1,            &
                                     -Fsupvmn(moff,noff,ns_mid)*t1
         END DO
      END DO

1000  FORMAT(' JS = ',i4,/                                                     &
             '    M     N      J^u         F^u         J^v         F^v')
1001  FORMAT(2i6, 1p,4e12.4)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out vector potential to file.
!>
!>  @param[in] js    Radial index.
!>  @param[in] iunit File unit to write to.
!-------------------------------------------------------------------------------
      SUBROUTINE DUMP_A(js, iunit)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: js
      INTEGER, INTENT(in) :: iunit

!  Local variables
      INTEGER             :: m
      INTEGER             :: n
      INTEGER             :: moff
      INTEGER             :: noff
      INTEGER             :: mp
      INTEGER             :: np
      INTEGER             :: nb
      INTEGER             :: is

!  Start of executable code
!  xc = xc+xc0 ON ENTRY
      nsmin = LBOUND(Asubsmnsf,3)
      nsmax = UBOUND(Asubsmnsf,3)
      moff  = LBOUND(Asubsmnsf,1)
      noff  = LBOUND(Asubsmnsf,2)

      IF (js .lt. nsmin .or. js .gt. nsmax) THEN
         RETURN
      END IF
      nb = MIN(2, ntor)

      IF (js .EQ. 1) THEN
         WRITE(iunit + 10, 1000)
         DO n = -nb, nb
            np = n + ntor + noff

            DO m = 0, 1
               mp = m + moff
               IF (n .lt. 0 .and. m .eq. 0) CYCLE
               WRITE(iunit + 10, 1001) n, m
               DO is = 1, ns
                  WRITE(iunit + 10,1002) is,                                   &
                                         Asubsmnsf(mp,np,is)*orthonorm(m,n),   &
                                         Asubumncf(mp,np,is)*orthonorm(m,n),   &
                                         Asubvmncf(mp,np,is)*orthonorm(m,n)
               END DO
            END DO
            WRITE(iunit + 10,*)
         END DO
      END IF

      WRITE(iunit, 1003) js, neqs
      DO n = -ntor, ntor
         np = n + ntor + noff
         DO m = 0, mpol
            mp = m + moff
            WRITE(iunit, 1004) m, n, Asubsmnsf(mp,np,js)*orthonorm(m,n),       &
                                     Asubumncf(mp,np,js)*orthonorm(m,n),       &
                                     Asubvmncf(mp,np,js)*orthonorm(m,n),       &
                                     jcurrumnc(m,n,js)*orthonorm(m,n),         &
                                     jcurrvmnc(m,n,js)*orthonorm(m,n)
            IF (.NOT.lasym) CYCLE
            WRITE(iunit + 100, 1004) m, n, Asubsmncf(mp,np,js)*orthonorm(m,n), &
                                           Asubumnsf(mp,np,js)*orthonorm(m,n), &
                                           Asubvmnsf(mp,np,js)*orthonorm(m,n), &
                                           jcurrumns(m,n,js)*orthonorm(m,n),   &
                                           jcurrvmns(m,n,js)*orthonorm(m,n)
         END DO
      END DO

      WRITE(iunit,*)

1000  FORMAT(' RAD     A_s         A_u         A_v')
1001  FORMAT(' N=',i2,' M=',i2)
1002  FORMAT(i4,3(1p,e12.4))
1003  FORMAT(' RADIUS(JS): ',i4,' NEQS: ',i6,                                  &
             /,'       M       N      A_s         A_u         A_v       ',     &
               'jcurru      jcurrv')
1004  FORMAT(2i8, 5(1p,e12.4))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out magnetic field to file.
!>
!>  @param[in] js    Radial index.
!>  @param[in] iunit File unit to write to.
!-------------------------------------------------------------------------------
      SUBROUTINE DUMP_B(js, iunit)
      USE metrics, ONLY: r1_i, chipf_i, phipf_i
      USE hessian, ONLY: gather_array
      USE siesta_namelist, ONLY: nsin

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN) :: js
      INTEGER, INTENT(IN) :: iunit

!  Local variables
      INTEGER             :: m
      INTEGER             :: n
      INTEGER             :: is
      INTEGER             :: i
      INTEGER             :: j
      INTEGER             :: n1
      INTEGER             :: n2
      REAL (dp)           :: rjac
      REAL (dp)           :: rbsupv
      REAL (dp)           :: rho
      REAL (dp)           :: t1

!  Start of executable code
!  xc = xc+xc0 ON ENTRY
!  Here, bsupX = jac*bsupX
      n1 = startglobrow
      n2 = endglobrow

      IF (js .EQ. 1) THEN
         ALLOCATE(jacobmnch(0:mpol,-ntor:ntor,ns))
         CALL tomnsp(jacobh, jacobmnch, f_cos)
         IF (lasym) THEN
            ALLOCATE(jacobmnsh(0:mpol,-ntor:ntor,ns))
            CALL tomnsp(jacobh, jacobmnsh, f_sin)
         END IF
 
         CALL Gather_Array(jbsupsmnsh)
         CALL Gather_Array(jbsupumnch)
         CALL Gather_Array(jbsupvmnch)
         IF (lasym) THEN
            CALL Gather_Array(jbsupsmnch)
            CALL Gather_Array(jbsupumnsh)
            CALL Gather_Array(jbsupvmnsh)
         END IF

         WRITE(iunit + 10,1000)
                           
         IF (endglobrow .EQ. NS) THEN
            is = ns
            DO m = 0, mpol
               DO n = -ntor, ntor
                  IF (n.lt.0 .AND. m.eq.0) CYCLE
                  t1 = -orthonorm(m,n)
                  WRITE(iunit + 10, 1001) m, n,                                &
                                          (1.5_dp*jbsupsmnsh(m,n,is) -         &
                                           0.5_dp*jbsupsmnsh(m,n,is - 1))*t1,  &
                                          (1.5_dp*jbsupumnch(m,n,is) -         &
                                           0.5_dp*jbsupumnch(m,n,is - 1))*t1,  &
                                          (1.5_dp*jbsupvmnch(m,n,is) -         &
                                           0.5_dp*jbsupvmnch(m,n,is - 1))*t1
               END DO
            END DO

            WRITE(iunit + 10,1002)
            DO n = -1,1
               DO m = 0, 1
                  IF (n.LT.0 .AND. m.EQ.0) CYCLE
                  t1 = -orthonorm(m,n)
                  WRITE(iunit + 10,1003) n, m
                  DO is = 2, ns
                     WRITE(iunit + 10, 1004)  is, jbsupsmnsh(m,n,is)*t1,       &
                                                  jbsupumnch(m,n,is)*t1,       &
                                                  jbsupvmnch(m,n,is)*t1
                  END DO
                  WRITE(iunit+10,*)
               END DO
            END DO

            CALL FLUSH(iunit + 10)
         END IF !endglobrow=ns

         RETURN
      END IF   !js=1

      IF (iam .EQ. 0) THEN
         WRITE(iunit, 1005) js, neqs

         WRITE(iunit, 1006)
         IF (lasym) THEN
            WRITE(iunit + 100, 1007)
         END IF

         DO n = -ntor, ntor
            DO m = 0, mpol
               t1 = -orthonorm(m,n)
               WRITE(iunit, 1008) m, n, jbsupsmnsh(m,n,js)*t1,                 &
                                        jbsupumnch(m,n,js)*t1,                 &
                                        jbsupvmnch(m,n,js)*t1,                 &
                                        jacobmnch(m,n,js)*t1
               IF (lasym) THEN
                  WRITE(iunit+100, 1008) m, n, jbsupsmnch(m,n,js)*t1,          &
                                               jbsupumnsh(m,n,js)*t1,          &
                                               jbsupvmnsh(m,n,js)*t1,          &
                                               jacobmnsh(m,n,js)*t1
               END IF
            END DO
         END DO

         WRITE(iunit,*)
         CALL FLUSH(iunit)
      END IF !iam=0

      IF (js .NE. ns) THEN
         RETURN
      END IF

!  Estimate iota profile on extended mesh. Note bsupx is jac*bupx here and minus
!  sign for sign of jacobian. FIXME: Sign of the jacobian is hard coded here.
!  This should be read from the wout file.
!      chipf_i(1) = 0.0
!      phipf_i(1) = 0.0
!      DO is = 2, ns - 1
      DO is = nsin, ns - 1
         chipf_i(is) = (-jbsupumnch(0,0,is + 1) - jbsupumnch(0,0,is))*0.5
         phipf_i(is) = (-jbsupvmnch(0,0,is + 1) - jbsupvmnch(0,0,is))*0.5
      END DO
      chipf_i(ns) = -2.5*jbsupumnch(0,0,ns) + 1.5*jbsupumnch(0,0,ns - 1)
      phipf_i(ns) = -2.5*jbsupvmnch(0,0,ns) + 1.5*jbsupvmnch(0,0,ns - 1)

!  Check Edge value
      IF (iam .EQ. 0) THEN

         rjac = 0
         rbsupv = 0
         DO n = -ntor, ntor
            DO m = 0, mpol
               rjac = rjac + jacobmnch(m,n,ns)*orthonorm(m,n)
               rbsupv = rbsupv + jbsupvmnch(m,n,ns)*orthonorm(m,n)
            END DO
         END DO

         rbsupv = rbsupv*r1_i(1,1,ns)/rjac

         WRITE (*,1009) rjac, jacobh(1,1,ns), rbsupv

!  Estimate iota profile. Note bsupx is jac*bsupx here and minus sign for sign
!  of jacobian. FIXME: Sign of the jacobian is hard coded here. This should be
!  read from the wout file.
         WRITE (*, 1010)
         DO is = 2, ns, 10 !Skip Every 10
!         DO is = 2, ns
             rho = (is-1.5_dp)*hs
             WRITE (*,1011) rho, 0.5_dp*(chipf_i(is) + chipf_i(is - 1)),       &
                                 -jbsupumnch(0,0,is),                          &
                                 -(asubvmncf(1,ntor + 1,is - 1) -              &
                                   asubvmncf(1,ntor + 1,is))*ohs,              &
                                 0.5_dp*(phipf_i(is) + phipf_i(is - 1)),       &
                                 -jbsupvmnch(0,0,is),                          &
                                 -(asubumncf(1,ntor + 1,is) -                  &
                                   asubumncf(1,ntor + 1,is-1))*ohs
         END DO
      END IF !iam = 0

      chipf_i(1) = 0.0
      phipf_i(1) = 0.0
      DO is = 2,nsin - 1
         chipf_i(is) = (-jbsupumnch(0,0,is + 1) - jbsupumnch(0,0,is))*0.5
         phipf_i(is) = (-jbsupvmnch(0,0,is + 1) - jbsupvmnch(0,0,is))*0.5
      END DO

1000  FORMAT(' EDGE VALUES OF jB^X FOR ALL M,N',/,                             &
             '   M   N   jB^s      jB^u      jB^v')
1001  FORMAT(2i4,3(1p,e10.2))
1002  FORMAT(/,' RAD    jB^s       jB^u       jB^v')
1003  FORMAT(' N=',i2,' M=',i2)
1004  FORMAT(i4,3(1p,e11.3))
1005  FORMAT('HALF RADIUS (JS): ',i4, ' NEQS: ',i6)
1006  FORMAT('       M       N   jB^s(sin)   jB^u(cos)   jB^v(cos)  jacob(cos)')
1007  FORMAT('       M       N   jB^s(cos)   jB^u(sin)   jB^v(sin)  jacob(sin)')
1008  FORMAT(2i8, 5(1p,e12.4))
1009  FORMAT(' jac(u=0,v=0,ns): ',1p,e12.3,' jacobh: ',1p,e12.3,/              &
             ' R*B^v(u=0,v=0,ns) : ',1p,e12.3)
1010  FORMAT(/,'  RAD      CHIP(WOUT)    CHIP(PCH)     CHIP(VP)'               &
               '      PHIP(WOUT)    PHIP(PCH)     PHIP(VP)')
!1011  FORMAT(f7.3,1p,4e14.3)
1011  FORMAT(f7.3,1p,6e14.3)

      END SUBROUTINE DUMP_B

!-------------------------------------------------------------------------------
!>  @brief The the current residual.
!-------------------------------------------------------------------------------
      SUBROUTINE GETFSQ_PCHELMS
      USE shared_data, ONLY: fsupsmnsf, fsupumncf, fsupvmncf,                  &
                             fsupsmncf, fsupumnsf, fsupvmnsf
      USE mpi_inc

      IMPLICIT NONE

!  Start of executable code
      nsmin = startglobrow
      nsmax = endglobrow
      
      fsq_total = SUM(Fsupsmnsf(:,:,nsmin:nsmax)**2 +                          &
                      Fsupumncf(:,:,nsmin:nsmax)**2 +                          &
                      Fsupvmncf(:,:,nsmin:nsmax)**2)
      IF (lasym) THEN
         fsq_total = fsq_total                                                 &
                   + SUM(Fsupsmncf(:,:,nsmin:nsmax)**2 +                       &
                         Fsupumnsf(:,:,nsmin:nsmax)**2 +                       &
                         Fsupvmnsf(:,:,nsmin:nsmax)**2)
      END IF
#if defined(MPI_OPT)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, fsq_total, 1, MPI_REAL8, MPI_SUM,       &
                         SIESTA_COMM, MPI_ERR)
#endif
      IF (bnorm .GT. zero) fsq_total = fsq_total/bnorm
      
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocate arrays.
!-------------------------------------------------------------------------------
      SUBROUTINE CLEANUP

      IMPLICIT NONE

!  Start of executable code
      DEALLOCATE(xc0, gc0)

      IF (ALLOCATED(jacobmnch)) THEN
         DEALLOCATE(jacobmnch)
      END IF
      IF (lasym .and. ALLOCATED(jacobmnsh)) THEN
         DEALLOCATE(jacobmnsh)
      END IF

      DEALLOCATE(Bsupsijh, Bsupuijh, Bsupvijh)
      DEALLOCATE(Bsubsijh, Bsubuijh, Bsubvijh)
      
      xc = 0
      gc = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Check that the line integral of Bsubuij equals the integrated
!>         toroidal current
!-------------------------------------------------------------------------------
      SUBROUTINE CHECK_CURRENT
      USE metrics, ONLY: jsupvdotA
      USE vmec_info, ONLY: Itor
      USE quantities, ONLY: b_factor

      IMPLICIT NONE

!  Local variables
      INTEGER   :: i
      REAL (dp) :: currv_int

!  Start of executable code
!  Average over toroidal angle too (should not dependent on toroidal angle)
      IF (endglobrow .NE. ns) RETURN
      
      currv_int = hs*(SUM(jcurrvmnc(0,0,1:ns - 1)) + jcurrvmnc(0,0,ns)/2)

      WRITE (*,1000) -twopi*line_Bu*orthonorm(0,0)/mu0, Itor
      WRITE (*,1001) line_bu, currv_int, jsupvdotA

1000  FORMAT(' Toroidal Current (SIESTA): '1p,e12.4,4x,1p,'(VMEC): ',1p,e12.4)
1001  FORMAT(' Line Integral B_u: ',1p,e12.4,4x,                               &
             ' Surface-integrated current: ',1p,e12.4,4x,                      &
             ' VMEC Surface-integrated current: ',1p,e12.4)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute linear forces.
!-------------------------------------------------------------------------------
      SUBROUTINE COMPUTE_FORCES_LIN

!  Start of executable code
      lHessian = .TRUE.
      CALL COMPUTE_FORCES(.FALSE.)
      lHessian = .FALSE.

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Update solution.
!-------------------------------------------------------------------------------
      SUBROUTINE BACKSLV
      USE metrics, ONLY: r1_i
      USE hessian, ONLY: gather_array

      IMPLICIT NONE

!  Start of executable code
!  update solution (already gathered to all processors)
      xc = xc + xc0

      CALL INIT_F(Fsupsmnsf, Fsupumncf, Fsupvmncf,                             &
                  jcurrumnc, jcurrvmnc)
      IF (lasym) THEN
         CALL INIT_F(Fsupsmncf, Fsupumnsf, Fsupvmnsf,                          &
                     jcurrumns, jcurrvmns)
      END IF
      gc0 = gc

      CALL COMPUTE_FORCES(.TRUE.)
      CALL gather_array(gc)

      IF (iam .EQ. 0) THEN
         WRITE(*,1000) fsq_total, MAXVAL(ABS(gc)), SQRT(SUM(xc*xc))
      END IF

1000  FORMAT(/,' Back Solve Check',                                            &
             /,' |F|^2: ',1p,e10.3,' MAX |F|: ', 1pe10.3, ' |X|: ',1pe10.3)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialize helmholtz solver.
!-------------------------------------------------------------------------------
      SUBROUTINE init_pchelms
      USE blocktridiagonalsolver_s, ONLY: Initialize
      USE descriptor_mod, ONLY: iam, nprocs
      USE nscalingtools, ONLY: initRemap
      USE shared_functions, ONLY: init_ptrs

      IMPLICIT NONE

!  Local variables
      INTEGER :: n1
      INTEGER :: istat

!  Start of executable code
      IF (lasym) THEN
         ndims = 6
      ELSE
         ndims = 3
      END IF

      mblk_size = ndims*mnmax
      CALL Initialize(ns, mblk_size)

      n1 = ns*mnmax
      neqs = ndims*n1

      CALL init_ptrs(xc, asubsmnsf, asubumncf, asubvmncf,                      &
                         asubsmncf, asubumnsf, asubvmnsf)

      CALL init_ptrs(gc, fsupsmnsf, fsupumncf, fsupvmncf,                      &
                         fsupsmncf, fsupumnsf, fsupvmnsf)

!  Initialize forces
      CALL initRemap(mpol, ntor, ns, nprocs, iam)

      CALL Compute_Forces(.TRUE.)
      xc = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute Fourier components of Curl(Curl(A)) - mu0J operator on full
!>         radial grid.
!>
!>  @param(in) lNLinear Use Non linear equations.
!-------------------------------------------------------------------------------
      SUBROUTINE COMPUTE_FORCES(lNLinear)
      USE shared_data, ONLY: asubsmnsf, asubumncf, asubvmncf,           &
                             fsupsmnsf, fsupumncf, fsupvmncf,           &
                             asubsmncf, asubumnsf, asubvmnsf,           &
                             fsupsmncf, fsupumnsf, fsupvmnsf, neqs
      USE bmw_run, ONLY: bmw_exec
      USE siesta_namelist, ONLY: wout_file, mgrid_file
      USE read_wout_mod, ONLY: Itor

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL, INTENT(IN) :: lNLinear

!  Local variables
      INTEGER                                  :: istat
      INTEGER                                  :: moff
      INTEGER                                  :: noff
      REAL (dp)                                :: dphi
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: A_r
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: A_p
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: A_z
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: cA_s
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: cA_u
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: cA_v

!  Start of executable code
      nsmin = startglobrow
      nsmax = MIN(ns,endglobrow + 1)
      ns_match = 0  !nsp /4

!  FIXME: All this is initalizion it should be somewhere else.
      IF (bnorm .LT. 0) THEN
         linit = .TRUE.
         IF (.NOT. ALLOCATED(xc0)) THEN
            ALLOCATE(xc0(neqs), gc0(neqs))
            CALL assert(ALLOCATED(jbsupsmnsh), 'B^S not allocated!')
            IF (lasym) THEN
               CALL assert(ALLOCATED(jbsupsmnch), 'B^S not allocated!')
            END IF
            ALLOCATE(Bsupsijh(ntheta,nzeta,nsmin:nsmax),                       &
                     Bsupuijh(ntheta,nzeta,nsmin:nsmax),                       &
                     Bsupvijh(ntheta,nzeta,nsmin:nsmax),                       &
                     Bsubsijh(ntheta,nzeta,nsmin:nsmax),                       &
                     Bsubuijh(ntheta,nzeta,nsmin:nsmax),                       &
                     Bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
            CALL assert_eq(0, istat, 'Allocation failed in COMPUTE_FORCES')
         END IF
    
         CALL INIT_F(Fsupsmnsf, Fsupumncf, Fsupvmncf,                          &
                     jcurrumnc, jcurrvmnc)
         IF (lasym) THEN
            CALL INIT_F(Fsupsmncf, Fsupumnsf, Fsupvmnsf,                       &
                        jcurrumns, jcurrvmns)
         END IF
         gc0 = gc
         
         CALL GetFSQ_PCHELMS
         bnorm = fsq_total

         dphi = (2*pi)/(nfp*nzeta)
         ALLOCATE(A_r(ntheta,nzeta,1),                                         &
                  A_p(ntheta,nzeta,1),                                         &
                  A_z(ntheta,nzeta,1))

         CALL bmw_exec(mgrid_file, wout_file, r1_i, z1_i, dphi, ns,            &
                       A_r, A_p, A_z)

         ALLOCATE(cA_s(ntheta,nzeta,1),                                        &
                  cA_u(ntheta,nzeta,1),                                        &
                  cA_v(ntheta,nzeta,1))
         CALL CYL2VMEC_A(A_r(:,:,1), A_p(:,:,1), A_z(:,:,1),                   &
                         cA_s(:,:,1), cA_u(:,:,1), cA_v(:,:,1))
         DEALLOCATE(A_r, A_p, A_z)

         CALL INIT_A(cA_s, cA_u, cA_v, asubsmnsf, asubumncf, asubvmncf, f_sin)
         IF (lasym) THEN
            CALL INIT_A(cA_s, cA_u, cA_v, asubsmncf, asubumnsf, asubvmnsf,     &
                        f_cos)
         END IF
         DEALLOCATE(cA_s, cA_u, cA_v)

         xc0 = xc
      END IF

      gc = 0

!  Get sqrt(g)*CURLA.grad(s,u,v) [sqrt(g) B^i] on half grid in real space, for
!  input AsubXmn on full grid.
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 1)
      CALL CURLA_PCHELMS(Asubsmnsf,  Asubumncf,  Asubvmncf,                    &
                         jbsupsmnsh, jbsupumnch, jbsupvmnch, f_sin)
      IF (lasym) THEN
         CALL CURLA_PCHELMS(Asubsmncf,  Asubumnsf,  Asubvmnsf,                 &
                            jbsupsmnch, jbsupumnsh, jbsupvmnsh, f_cos)
      END IF

!  Add 1/jac factor
      nsmin = startglobrow

      Bsupsijh = Bsupsijh/jacobh(:,:,nsmin:nsmax)
      Bsupuijh = Bsupuijh/jacobh(:,:,nsmin:nsmax)
      Bsupvijh = Bsupvijh/jacobh(:,:,nsmin:nsmax)
      
!  Convert to covariant component of B_i on half grid.
      CALL tolowerh(Bsupsijh, Bsupuijh, Bsupvijh,                              &
                    Bsubsijh, Bsubuijh, Bsubvijh, nsmin, nsmax)

!  Get -sqrt(g)*CURLB dot grad(s,u,v) Fourier coefficients [sqrt(g)*J^i] on full
!  grid and store in FsupXmnf (sign is flipped to make positive contributions to
!  Hessian)
      nsmin = startglobrow
      nsmax = MIN(ns, endglobrow + 1)
      CALL CURLB_PCHELMS(Fsupsmnsf, Fsupumncf, Fsupvmncf, Asubsmnsf,           &
                         jbsupsmnsh, f_sin)
      IF (lasym) THEN
         CALL CURLB_PCHELMS(Fsupsmncf, Fsupumnsf, Fsupvmnsf, Asubsmncf,        &
                            jbsupsmnch, f_cos)
      END IF
      nsmax = endglobrow
      
!  Dump debugging info.
!  FIXME: This should only get written out if compiled in debug mode. Check
!  linit to see if this is only written once.
      IF (linit) THEN
         CALL Compare_Current(Fsupumncf, Fsupvmncf, jcurrumnc, jcurrvmnc)
         IF (lasym) THEN
            CALL Compare_Current(Fsupumnsf, Fsupvmnsf, jcurrumns, jcurrvmns)
         END IF
      END IF

      IF (linit .OR. lNlinear) THEN
         gc = gc + gc0
      END IF

!  Compute ONLY -Ax, NOT residue= -Ax + b (for preconditioner and GMRES)
!  Subtract [sqrt(g)*J]mn (input from VMEC) to get net force (residual=b-Ax)
      CALL BoundaryConditions(Fsupsmnsf, Fsupumncf, Fsupvmncf, f_sin)
      IF (lasym) THEN
         CALL BoundaryConditions(Fsupsmncf, Fsupumnsf, Fsupvmnsf, f_cos)
      END IF
      
!  Subtract A*x0 since A_u, A_v are fixed at the edge b' = b + A*x0 (new
!  right-hand side)
      IF (linit) THEN
         gc0 = gc
      END IF

      IF (lNLinear .OR. lInit) THEN 
         linit = .FALSE.
         CALL GETFSQ_PCHELMS
      END IF

      END SUBROUTINE
      
!-------------------------------------------------------------------------------
!>  @brief Apply boundary conditions.
!>
!>  @param(inout) Fsupsmn Contravariant current in the s direction.
!>  @param(inout) Fsupumn Contravariant current in the u direction.
!>  @param(inout) Fsupvmn Contravariant current in the v direction.
!>  @param(in)    iparity Fourier parity of the residuals.
!-------------------------------------------------------------------------------
      SUBROUTINE BoundaryConditions(Fsupsmn, Fsupumn, Fsupvmn, iparity)
      USE hessian, ONLY: l_Compute_Hessian
      USE utilities, ONLY: set_bndy_fouier_m0, set_bndy_full_origin

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: Fsupsmn
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: Fsupumn
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: Fsupvmn
      INTEGER, INTENT(in)                                      :: iparity

!  Local variables
      INTEGER                                                  :: nmatch

!  Start of executable code
      IF (endglobrow .EQ. ns) THEN
         IF (.not.l_AsEdge) THEN
            fsupsmn(:,:,ns) = 0
         END IF
         fsupumn(:,:,ns) = 0
         fsupvmn(:,:,ns) = 0
      END IF

      CALL set_bndy_fouier_m0(Fsupsmn, Fsupumn, Fsupvmn, iparity)
      IF (startglobrow .eq. 1) THEN
         CALL set_bndy_full_origin(Fsupsmn, Fsupumn, Fsupvmn)
      END IF

!  Preserve surfaces for js < ns_match
      IF (ns_match .GE. nsmin) THEN
         nmatch = MIN(ns_match, nsmax)
         fsupumn(m1,:,nsmin:nmatch) = 0
         fsupvmn(m1,:,nsmin:nmatch) = 0
      END IF
      
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Setup and run GMRES solver for the helmholtz problem.
!>
!>  GMRES solves r == b + Ax = 0
!-------------------------------------------------------------------------------
      SUBROUTINE GMRES_PCHELMS
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE gmres_lib, ONLY: gmres_info, gmres_par, gmres_ser
      USE hessian, ONLY: Apply_Precond
      USE shared_data, ONLY: lverbose, unit_out, gmres_no_peak, iortho

      IMPLICIT NONE

!  Local variables
      TYPE (gmres_info)                    :: gi
      INTEGER                              :: n
      INTEGER                              :: m
      INTEGER                              :: js
      REAL (dp), DIMENSION(:), ALLOCATABLE :: b
      REAL (dp), DIMENSION(:), ALLOCATABLE :: x0

!  Local Parameters
      INTEGER, PARAMETER                   :: noPrec = 0
      INTEGER, PARAMETER                   :: leftPrec = 1
      INTEGER, PARAMETER                   :: rightPrec = 2
      INTEGER, PARAMETER                   :: dblePrec = 3
      REAL (dp), PARAMETER                 :: one = 1
      REAL (dp), PARAMETER                 :: ftol = 1.E-18_dp

!  Start of executable code
      ALLOCATE(b(neqs), x0(neqs))

      xc = 0

      CALL Compute_Forces(.TRUE.)
      IF (iam .EQ. 0) THEN
         IF (lverbose) THEN
            WRITE (*,*)
            WRITE (*,1001) 1, fsq_total, SQRT(SUM(xc0**2))
         END IF
         WRITE (unit_out, 1000) neqs
         WRITE (unit_out, 1001) 1, fsq_total, SQRT(SUM(xc0**2))
      END IF

      b = gc0
      x0 = xc
      x0 = 0

      CALL init_dgmres(gi%icntl, gi%cntl)

!  Tune GMRES parameters
!  Tolerance
       gi%cntl(1) = SQRT(ftol)
!  Write warnings to fort.21
       gi%icntl(2) = 21
!  Save the convergence history in file fort.20
       gi%icntl(3) = 20
!  Preconditioning input flags (note: different than revcom flags in driver)
      gi%icntl(4) = rightPrec ! leftPrec, dblePrec, noPrec
!  Orthogonalization type.
      gi%icntl(5) = iortho
!  Initial guess (use it if = 1, otherwise ignore it)
      gi%icntl(6) = 1
!  Maximum number of iterations at each step (~ns/5)
      gi%icntl(7) = 20
!  Default
      gi%icntl(8) = 1
!  Steps for peek at progress during rev com loop.
      gi%icntl(9) = 1
!  linear solver (output on screen).
      gi%l_nonlinear = .FALSE.


      n   = neqs
      gi%m = 200

!  SIESTA process mpi rank.
      gi%iam = iam
!  Number of siesta processes.
      gi%nprocs=nprocs
      gi%ngmres_type = gmres_no_peak
      gi%ftol = ftol
      IF (PARSOLVER) THEN
         gi%endglobrow = endglobrow
         gi%startglobrow = startglobrow
         gi%mblk_size = mblk_size
         gi%rcounts => rcounts
         gi%disp => disp

         gi%my_comm = SIESTA_COMM
         gi%my_comm_world = SIESTA_COMM

         CALL gmres_par(n, gi, matvec_par_PCHELMS, apply_precond,              &
                        GetNLForce_PCHELMS, x0, b)
      ELSE
         CALL gmres_ser(n, gi, matvec_PCHELMS, apply_precond,                  &
                        GetNLForce_PCHELMS, x0, b)
      END IF
      
      xc = x0
      
      DEALLOCATE(x0, b)

1000  FORMAT(' NEQS : ',i6)
1001  FORMAT(' ITER: ',i6,3x,'|F|^2: ',1pE12.3,3x,' |X| = ',1pE12.3)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Serial callback function for MatVec GMRES routine.
!>
!>  @param[in]  p    Displacement parameters.
!>  @param[out] Ap   Matrix element.
!>  @param[in]  ndim Number of dimensiona.
!-------------------------------------------------------------------------------
      SUBROUTINE MatVec_PCHELMS (p, Ap, ndim)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(IN)  :: p(ndim)
      REAL(dp), INTENT(OUT) :: Ap(ndim)
      INTEGER, INTENT(IN)   :: ndim

!  Start of executable code
      CALL assert_eq(neqs, ndim, ' neqs != ndim')

      xc = p
      CALL Compute_Forces(.FALSE.)
      Ap = -gc

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Parallel callback function for MatVec GMRES routine.
!>
!>  @param[in]  ploc Local displacement parameters.
!>  @param[out] Ap   Local Matrix element.
!>  @param[in]  nloc Local number of dimensiona.
!-------------------------------------------------------------------------------
      SUBROUTINE MatVec_par_PCHELMS (ploc, Ap, nloc)
      USE hessian, ONLY: eps_factor
      USE blocktridiagonalsolver_s, ONLY: ParMatVec, PadSides 

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nloc), INTENT(IN)  :: ploc
      REAL(dp), DIMENSION(nloc), INTENT(OUT) :: Ap
      INTEGER, INTENT(IN)                    :: nloc

!  local variables
      INTEGER                                :: myrowstart
      INTEGER                                :: myrowend
      INTEGER                                :: istat
      REAL (dp), DIMENSION(:), ALLOCATABLE   :: p

!  Start of executable code
      istat = (endglobrow - startglobrow + 1)*mblk_size
      CALL assert_eq(istat, nloc, 'nloc wrong in matvec_par')

      myrowstart = (startglobrow - 1)*mblk_size + 1
      myrowend = myrowstart + nloc - 1

      ALLOCATE (p(ns*mblk_size), stat=istat)
      CALL assert_eq(0, istat, 'Allocation error in matvec_par')

      p(myrowstart:myrowend) = ploc
      CALL PadSides(p, mblk_size, 1, 1)

      xc = p
      CALL Compute_Forces(.false.)
      Ap = -gc(myrowstart:myrowend)

      DEALLOCATE (p)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Non linear force callback.
!>
!>  @param[in]  xcstate Displacement parameters.
!>  @param[out] fsq_nl  |F^2| non-linear residual.
!>  @param[in]  bnorm   Internal GMRES normalization
!-------------------------------------------------------------------------------
      SUBROUTINE GetNLForce_PCHELMS(xcstate, fsq_nl, bnorm)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(neqs), INTENT(IN) :: xcstate
      REAL (dp), INTENT(OUT)                 :: fsq_nl
      REAL (dp), INTENT(IN)                  :: bnorm

!  Start of executable code
	  WRITE (*,1000)
	  fsq_nl = 1
	  RETURN

1000  FORMAT('PCHELMS should not get here.')

	  END SUBROUTINE

      END MODULE pchelms
