!*******************************************************************************
!>  \brief Module contained subroutines for computing the real-space
!>  contravariant (full mesh) "displacements" as part of the SIESTA project.
!>  \author S. P. Hirshman and R. Sanchez
!>  \date Aug 18, 2006
!>  added l_asym=T Dec, 2016 M. R. Cianciosa
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
      USE shared_data, ONLY: xc, l_push_s, l_push_u, l_push_v,           &
                             l_push_edge, l_natural, l_asym, col_scale
      USE timer_mod, ONLY: time_update_upperv
      USE descriptor_mod, ONLY: iam  

      CONTAINS
!>    \brief Parallel subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv_par
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE siesta_state, ONLY: clear_field_perts
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, m2=2
      INTEGER ::  nsmin, nsmax, ns_span, istat
      REAL(dp) :: ton, toff, r0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
                gvsupsmncf, gvsupumnsf, gvsupvmnsf,                      &
                gvsupsmnsf, gvsupumncf, gvsupvmncf
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow-2); nsmax=MIN(endglobrow+2,ns)
      ns_span = nsmax-nsmin+1

      CALL Clear_Field_Perts

!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components
       ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                 gvsupumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                 gvsupvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error #1 in UPDATE_UPPERV'

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

!Edge boundary conditions (B_tangential = 0 => vsups = 0)
       IF (nsmax .EQ. ns) THEN
          gvsupsmncf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumnsf(:,:,ns) = 0
             gvsupvmnsf(:,:,ns) = 0
          END IF
       END IF  !!nsmax == ns

       CALL toijsp_par(gvsupsmncf, jvsupsijf(:,:,nsmin:nsmax), f_none,   &  ! Calculate covariant velocities in real space (full mesh)
                       f_cos, f_full, 1, ns_span)                         
       CALL toijsp_par(gvsupumnsf, jvsupuijf(:,:,nsmin:nsmax), f_none,   &
                       f_sin, f_full, 1, ns_span)                        
       CALL toijsp_par(gvsupvmnsf, jvsupvijf(:,:,nsmin:nsmax), f_none,   &
                       f_sin, f_full, 1, ns_span)

       DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf,stat=istat)

       IF (.NOT.l_asym) GOTO 100

       ALLOCATE (gvsupsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                 gvsupumncf(0:mpol,-ntor:ntor,nsmin:nsmax),              &
                 gvsupvmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)

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
             gvsupsmnsf(m1,:,1) = r0*gvsupumncf(m1,:,1)
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

       CALL toijsp_par(gvsupsmnsf, jvsupsijf(:,:,nsmin:nsmax),           &
                       f_sum, f_sin, f_full, 1, ns_span)
       CALL toijsp_par(gvsupumncf, jvsupuijf(:,:,nsmin:nsmax),           &
                       f_sum, f_cos, f_full, 1, ns_span)
       CALL toijsp_par(gvsupvmncf, jvsupvijf(:,:,nsmin:nsmax),           &
                       f_sum, f_cos, f_full, 1, ns_span)

       DEALLOCATE (gvsupsmnsf, gvsupumncf, gvsupvmncf)

 100   CONTINUE
       CALL second0(toff)
       time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv_par

!>    \brief Serial subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv

      CALL update_upperv_par

      END SUBROUTINE update_upperv

      END MODULE siesta_displacement
