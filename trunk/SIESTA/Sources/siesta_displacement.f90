!*******************************************************************************
!>  \brief Module contained subroutines for computing the real-space
!>  contravariant (full mesh) "displacements" as part of the SIESTA project.
!>  \author S. P. Hirshman and R. Sanchez
!>  \date Aug 18, 2006
!>  added lasym=T Dec, 2016 M. R. Cianciosa
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
      USE nscalingtools, ONLY: startglobrow, endglobrow, ErrorAbort
      USE timer_mod, ONLY: time_update_upperv
      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
      INTEGER, PRIVATE ::  nsmin, nsmax, ns_span

      CONTAINS
!>    \brief Parallel subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv_par
      USE shared_data, ONLY: lasym
      USE siesta_state, ONLY: clear_field_perts
      USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,          &
                            jvsupsmnsf, jvsupumncf, jvsupvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  istat
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)
      ns_span = nsmax-nsmin+1

      CALL Clear_Field_Perts

!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components

      CALL InitDisplacement(jvsupsmncf, jvsupumnsf, jvsupvmnsf, isym)
      IF (lasym)                                                         &
      CALL InitDisplacement(jvsupsmnsf, jvsupumncf, jvsupvmncf, iasym)

      CALL second0(toff)
      time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv_par

!>    \brief Serial subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv

      CALL update_upperv_par

      END SUBROUTINE update_upperv

 
      SUBROUTINE InitDisplacement(jvsupsmnf, jvsupumnf, jvsupvmnf,       &
                                  iparity)
      USE quantities, ONLY: jvsupsijf, jvsupuijf, jvsupvijf
      USE shared_data, ONLY: l_push_s, l_push_u, l_push_v,               &
                             l_push_edge, l_natural, col_scale
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:,:), POINTER  ::                            &
                jvsupsmnf, jvsupumnf, jvsupvmnf
      INTEGER, INTENT(IN)  :: iparity
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, m2=2
      INTEGER  :: istat, m, n, c1, c2, c3,                               &
                  fours, fouru, fourv, fcomb 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                        &
                gvsupsmnf, gvsupumnf, gvsupvmnf
      REAL(dp) :: r0
!-----------------------------------------------
      IF (iparity .EQ. isym) THEN
         fours=f_cos; fouru=f_sin; fourv=f_sin; fcomb=f_none
         c1 = 1; r0 = hs_i/2
      ELSE
         fours=f_sin; fouru=f_cos; fourv=f_cos; fcomb=f_sum
         c1 = 4; r0 = -hs_i/2
      END IF
      c2 = c1+1;  c3 = c2+1

      IF (LBOUND(jvsupsmnf,3) .NE. 1)  CALL ErrorAbort('LBOUND1 WRONG IN InitDisp')
      IF (UBOUND(jvsupsmnf,3) .NE. ns) CALL ErrorAbort('UBOUND1 WRONG IN InitDisp')

      ALLOCATE (gvsupsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                gvsupumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                gvsupvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      IF (istat .NE. 0) CALL ErrorAbort('Allocation error #1 in UPDATE_UPPERV')

!Apply scaling in real space (keep xc components unscaled)
      gvsupsmnf = jvsupsmnf(:,:,nsmin:nsmax)*col_scale(:,:,c1,nsmin:nsmax)
      gvsupumnf = jvsupumnf(:,:,nsmin:nsmax)*col_scale(:,:,c2,nsmin:nsmax)
      gvsupvmnf = jvsupvmnf(:,:,nsmin:nsmax)*col_scale(:,:,c3,nsmin:nsmax)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
      NS1: IF (nsmin .EQ. 1) THEN
         gvsupsmnf(m0,:,1) = 0;  gvsupsmnf(m2:,:,1) = 0
         IF (.NOT.l_push_s) gvsupsmnf(m1,:,1) = 0

         gvsupumnf(m0,:,1) = 0;  gvsupumnf(m2:,:,1) = 0
         IF (.NOT.l_push_u) gvsupumnf(m1,:,1) = 0

         gvsupvmnf(m1:,:,1) = 0
         IF (.NOT.l_push_v) gvsupvmnf(m0,:,1) = 0

         IF (.NOT. l_natural) THEN
            gvsupsmnf(m1,:,1) = -r0*gvsupumnf(m1,:,1)
         END IF
       END IF NS1   

!m=0, n<0 (need to avoid singular Hessian)
       gvsupsmnf(m0,-ntor:-1,nsmin:nsmax) = 0       
       gvsupumnf(m0,-ntor:-1,nsmin:nsmax) = 0       
       gvsupvmnf(m0,-ntor:-1,nsmin:nsmax) = 0       

!Edge boundary conditions (B_tangential = 0 => vsups = 0)
       NSM: IF (nsmax .EQ. ns) THEN
          gvsupsmnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumnf(:,:,ns) = 0
             gvsupvmnf(:,:,ns) = 0
          END IF
       END IF NSM

       CALL toijsp_par(gvsupsmnf, jvsupsijf(:,:,nsmin:nsmax), fcomb,     &  ! Calculate covariant velocities in real space (full mesh)
                       fours, 1, ns_span)                         
       CALL toijsp_par(gvsupumnf, jvsupuijf(:,:,nsmin:nsmax), fcomb,     &
                       fouru, 1, ns_span)                        
       CALL toijsp_par(gvsupvmnf, jvsupvijf(:,:,nsmin:nsmax), fcomb,     &
                       fourv, 1, ns_span)

       DEALLOCATE (gvsupsmnf, gvsupumnf, gvsupvmnf,stat=istat)


      END SUBROUTINE InitDisplacement


      END MODULE siesta_displacement
