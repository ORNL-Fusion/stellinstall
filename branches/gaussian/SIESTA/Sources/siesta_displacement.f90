!*******************************************************************************
!>  \brief Module contained subroutines for computing the real-space 
!>  contravariant (full mesh) "displacements" as part of the SIESTA project.
!>  \author S. P. Hirshman and R. Sanchez
!>  \date Aug 18, 2006
!>
!>  Displacements are defined by 
!>  
!>     J*v*dt                                                                (1)
!>
!>  where J is the coordinate system Jacobian and v is the dispalcement 
!>  velocity. The velocity harmonics are obtained from solving the Fourier 
!>  components of the MHD force. Fourier harmonics of the displacement velocity 
!>  are converted to real space.
!*******************************************************************************
      MODULE siesta_displacement
      USE stel_kinds
      USE stel_constants        
      USE fourier 
      USE quantities
      USE shared_data, ONLY: xc, l_push_s, l_push_u, l_push_v,        &
                             l_push_edge, l_natural, col_scale
      USE timer_mod, ONLY: time_update_upperv
      USE descriptor_mod, ONLY: iam  

      CONTAINS
#if defined(SKS)
!*******************************************************************************
!>  \brief Parallel subroutine for initializing updated displacement
!>  components in real space.
!*******************************************************************************
      SUBROUTINE update_upperv_par
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE siesta_state, ONLY: clear_field_perts
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0,m1=1,m2=2,ifull=0,pcos=0,psin=1
      INTEGER ::  js, nsmin, nsmax, ns_span, istat
      REAL(dp) :: ton, toff, r0
      INTEGER :: i,j,k
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
#if defined(MRC)
                gvsupsmnsf, gvsupumncf, gvsupvmncf, work,              &
#endif
                gvsupsmncf, gvsupumnsf, gvsupvmnsf

!-----------------------------------------------
      CALL second0(ton)

      nsmin=MAX(1,startglobrow-2); nsmax=MIN(endglobrow+2,ns)
      ns_span = nsmax-nsmin+1

      CALL Clear_Field_Perts

!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components
       ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                 gvsupumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                 gvsupvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE (gvsupsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                    gvsupumncf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                    gvsupvmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          ALLOCATE (work(SIZE(jvsupsijf, 1), SIZE(jvsupsijf, 2), SIZE(jvsupsijf, 3)), stat=istat)
       END IF
#endif
       IF (istat .ne. 0) STOP 'Allocation error # 1in UPDATE_UPPERV'

!Apply scaling in real space (keep xc components unscaled)
       gvsupsmncf = jvsupsmncf(:,:,nsmin:nsmax)*col_scale(:,:,1,nsmin:nsmax)
       gvsupumnsf = jvsupumnsf(:,:,nsmin:nsmax)*col_scale(:,:,2,nsmin:nsmax)
       gvsupvmnsf = jvsupvmnsf(:,:,nsmin:nsmax)*col_scale(:,:,3,nsmin:nsmax)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
       IF (nsmin .EQ. 1) THEN
          gvsupsmncf(m0,:,1) = 0;  gvsupsmncf(m2:,:,1) = 0
          IF (.NOT.l_push_s) gvsupsmncf(m1,:,1) = 0

          gvsupumnsf(m0,:,1) = 0;  gvsupumnsf(m2:,:,1) = 0
          IF (.NOT.l_push_u) gvsupumnsf(m1,:,1) = 0

          gvsupvmnsf(m1:,:,1) = 0
          IF (.NOT.l_push_v) gvsupvmnsf(m0,:,1) = 0

          IF (.NOT. l_natural) THEN
             r0 = hs_i/2
             gvsupsmncf(m1,:,1) = -r0*gvsupumnsf(m1,:,1)
          END IF
       END IF   !nsmin == 1

!m=0, n<0 (need so Hessian catch-all loop gets this...)
       gvsupsmncf(m0,-ntor:-1,nsmin:nsmax) = 0
       gvsupumnsf(m0,-ntor:0,nsmin:nsmax) = 0
       gvsupvmnsf(m0,-ntor:0,nsmin:nsmax) = 0

!Edge boundary conditions
       IF (nsmax .EQ. ns) THEN
          gvsupsmncf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumnsf(:,:,ns) = 0
             gvsupvmnsf(:,:,ns) = 0
          END IF
       END IF  !!nsmax == ns

       CALL toijsp_par(gvsupsmncf, jvsupsijf(:,:,nsmin:nsmax), 0, 0, pcos, ifull, 1, ns_span)    ! Calculate covariant velocities in real space (full mesh)
       CALL toijsp_par(gvsupumnsf, jvsupuijf(:,:,nsmin:nsmax), 0, 0, psin, ifull, 1, ns_span)    ! NOTE: s-component is COS, remaining are SIN
       CALL toijsp_par(gvsupvmnsf, jvsupvijf(:,:,nsmin:nsmax), 0, 0, psin, ifull, 1, ns_span)

       DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf,stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          gvsupsmnsf = jvsupsmnsf(:,:,nsmin:nsmax)*col_scale(:,:,4,nsmin:nsmax)
          gvsupumncf = jvsupumncf(:,:,nsmin:nsmax)*col_scale(:,:,5,nsmin:nsmax)
          gvsupvmncf = jvsupvmncf(:,:,nsmin:nsmax)*col_scale(:,:,6,nsmin:nsmax)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
          IF (nsmin .EQ. 1) THEN
             gvsupsmnsf(m0,:,1) = 0;  gvsupsmnsf(m2:,:,1) = 0
             IF (.NOT.l_push_s) gvsupsmnsf(m1,:,1) = 0

             gvsupumncf(m0,:,1) = 0;  gvsupumncf(m2:,:,1) = 0
             IF (.NOT.l_push_u) gvsupumncf(m1,:,1) = 0

             gvsupvmncf(m1:,:,1) = 0
             IF (.NOT.l_push_v) gvsupvmncf(m0,:,1) = 0

             IF (.NOT. l_natural) THEN
                r0 = hs_i/2
                gvsupsmnsf(m1,:,1) = -r0*gvsupumncf(m1,:,1)
             END IF
          END IF   !nsmin == 1

!m=0, n<0 (need so Hessian catch-all loop gets this...)
          gvsupsmnsf(m0,-ntor:0,nsmin:nsmax) = 0
          gvsupumncf(m0,-ntor:-1,nsmin:nsmax) = 0
          gvsupvmncf(m0,-ntor:-1,nsmin:nsmax) = 0

!Edge boundary conditions
          IF (nsmax .EQ. ns) THEN
             gvsupsmnsf(:,:,ns) = 0
             IF (.NOT.l_push_edge) THEN
                gvsupumncf(:,:,ns) = 0
                gvsupvmncf(:,:,ns) = 0
             END IF
          END IF  !!nsmax == ns

          CALL toijsp_par(gvsupsmnsf, work(:,:,nsmin:nsmax), fourier_none, fourier_none, &
                          fourier_parity_odd, fourier_full_grid, 1, ns_span)    ! Calculate covariant velocities in real space (full mesh)
          jvsupsijf = jvsupsijf + work

          CALL toijsp_par(gvsupumncf, work(:,:,nsmin:nsmax), fourier_none, fourier_none, &
                          fourier_parity_even, fourier_full_grid, 1, ns_span)    ! NOTE: s-component is SIN, remaining are COSINE
          jvsupuijf = jvsupuijf + work

          CALL toijsp_par(gvsupvmncf, work(:,:,nsmin:nsmax), fourier_none, fourier_none, &
                          fourier_parity_even, fourier_full_grid, 1, ns_span)
          jvsupvijf = jvsupvijf + work

          DEALLOCATE (gvsupsmnsf, gvsupumncf, gvsupvmncf,stat=istat)
       END IF
#endif

       CALL second0(toff)
       time_update_upperv = time_update_upperv + (toff-ton)

       END SUBROUTINE update_upperv_par
#endif

!*******************************************************************************
!>  \brief Serial subroutine for initializing updated displacement components
!>  in real space.
!*******************************************************************************
      SUBROUTINE update_upperv
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0,m1=1,m2=2,ifull=0,pcos=0,psin=1
      REAL(dp) :: ton, toff, r0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
#if defined(MRC)
                          gvsupsmnsf, gvsupumncf, gvsupvmncf, work,    &
#endif
                          gvsupsmncf, gvsupumnsf, gvsupvmnsf
      INTEGER :: i, j, k, nsmin, nsmax, istat
!-----------------------------------------------
!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components
       CALL second0(ton)

       ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,ns),                      &
                 gvsupumnsf(0:mpol,-ntor:ntor,ns),                      &
                 gvsupvmnsf(0:mpol,-ntor:ntor,ns), stat=istat)

#if defined(MRC)
       IF (l_force_asym) THEN
          ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,ns),                      &
                    gvsupumnsf(0:mpol,-ntor:ntor,ns),                      &
                    gvsupvmnsf(0:mpol,-ntor:ntor,ns), stat=istat)
          ALLOCATE (work(SIZE(jvsupsijf, 1), SIZE(jvsupsijf, 2), SIZE(jvsupsijf, 3)), stat=istat)
       END IF
#endif
       IF (istat .ne. 0) STOP 'Allocation error # 1in UPDATE_UPPERV'

!Apply scaling in real space (keep xc components unscaled)
       gvsupsmncf = jvsupsmncf*col_scale(:,:,1,:)
       gvsupumnsf = jvsupumnsf*col_scale(:,:,2,:)
       gvsupvmnsf = jvsupvmnsf*col_scale(:,:,3,:)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
       gvsupsmncf(m0,:,1) = 0;  gvsupsmncf(m2:,:,1) = 0
       IF (.NOT.l_push_s) gvsupsmncf(m1,:,1) = 0

       gvsupumnsf(m0,:,1) = 0;  gvsupumnsf(m2:,:,1) = 0
       IF (.NOT.l_push_u) gvsupumnsf(m1,:,1) = 0

       gvsupvmnsf(m1:,:,1) = 0
       IF (.NOT.l_push_v) gvsupvmnsf(m0,:,1) = 0

       IF (.NOT. l_natural) THEN
          r0 = hs_i/2
          gvsupsmncf(m1,:,1) = -r0*gvsupumnsf(m1,:,1)
       END IF

!m=0, n<0 (need so Hessian catch-all loop gets this...)
       gvsupsmncf(m0,-ntor:-1,:) = 0       
       gvsupumnsf(m0,-ntor:0,:) = 0       
       gvsupvmnsf(m0,-ntor:0,:) = 0       

!Edge boundary conditions
       gvsupsmncf(:,:,ns) = 0
       IF (.NOT.l_push_edge) THEN
          gvsupumnsf(:,:,ns) = 0
          gvsupvmnsf(:,:,ns) = 0
       END IF

       CALL toijsp(gvsupsmncf, jvsupsijf, 0, 0, pcos, ifull)    ! Calculate upper velocities in real space (full mesh)
       CALL toijsp(gvsupumnsf, jvsupuijf, 0, 0, psin, ifull)    ! NOTE: s-component is COS, remaining are SIN
       CALL toijsp(gvsupvmnsf, jvsupvijf, 0, 0, psin, ifull)

       DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf)

#if defined(MRC)
       IF (l_force_asym) THEN
!Apply scaling in real space (keep xc components unscaled)
          gvsupsmnsf = jvsupsmncf*col_scale(:,:,1,:)
          gvsupumncf = jvsupumnsf*col_scale(:,:,2,:)
          gvsupvmncf = jvsupvmnsf*col_scale(:,:,3,:)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
          gvsupsmnsf(m0,:,1) = 0;  gvsupsmncf(m2:,:,1) = 0
          IF (.NOT.l_push_s) gvsupsmncf(m1,:,1) = 0

          gvsupumncf(m0,:,1) = 0;  gvsupumncf(m2:,:,1) = 0
          IF (.NOT.l_push_u) gvsupumncf(m1,:,1) = 0

          gvsupvmncf(m1:,:,1) = 0
          IF (.NOT.l_push_v) gvsupvmncf(m0,:,1) = 0

          IF (.NOT. l_natural) THEN
             r0 = hs_i/2
             gvsupsmnsf(m1,:,1) = -r0*gvsupumncf(m1,:,1)
          END IF

!m=0, n<0 (need so Hessian catch-all loop gets this...)
          gvsupsmnsf(m0,-ntor:0,:) = 0
          gvsupumncf(m0,-ntor:-1,:) = 0
          gvsupvmncf(m0,-ntor:-1,:) = 0

!Edge boundary conditions
          gvsupsmnsf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumncf(:,:,ns) = 0
             gvsupvmncf(:,:,ns) = 0
          END IF

          CALL toijsp(gvsupsmnsf, work, fourier_none, fourier_none, fourier_parity_even, fourier_full_grid)    ! Calculate upper velocities in real space (full mesh)
          jvsupsijf = jvsupsijf + work

          CALL toijsp(gvsupumncf, work, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid)    ! NOTE: s-component is SIN, remaining are COS
          jvsupuijf = jvsupuijf + work

          CALL toijsp(gvsupvmncf, work, fourier_none, fourier_none, fourier_parity_odd, fourier_full_grid)
          jvsupvijf = jvsupvijf + work

          DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf)
      END IF
#endif

       CALL second0(toff)
       time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv

      END MODULE siesta_displacement
