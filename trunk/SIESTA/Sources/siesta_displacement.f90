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
      MODULE siesta_displacement
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants        
      USE fourier 
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE timer_mod, ONLY: time_update_upperv
      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
      INTEGER, PRIVATE ::  nsmin, nsmax

      CONTAINS
!>    \brief Parallel/Serial subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv
      USE shared_data, ONLY: lasym
      USE siesta_state, ONLY: clear_field_perts
      USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,         &
                            jvsupsmnsf, jvsupumncf, jvsupvmncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  istat
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)

      CALL Clear_Field_Perts

!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components

      CALL InitDisplacement(jvsupsmncf, jvsupumnsf, jvsupvmnsf, isym)
      IF (lasym)                                                        &
      CALL InitDisplacement(jvsupsmnsf, jvsupumncf, jvsupvmncf, iasym)

      CALL second0(toff)
      time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv
      
 
      SUBROUTINE InitDisplacement(jvsupsmnf, jvsupumnf, jvsupvmnf, ipar)
      USE quantities, ONLY: jvsupsijf, jvsupuijf, jvsupvijf
      USE shared_data, ONLY: l_push_s, l_push_u, l_push_v,              &
                             l_push_edge, l_natural, col_scale
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:,:), POINTER  ::                           &
                jvsupsmnf, jvsupumnf, jvsupvmnf
      INTEGER, INTENT(IN)  :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, m2=2
      INTEGER  :: istat, m, n, mcol,                                    &
                  fours, fouru, fourv, fcomb 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
                gvsupsmnf, gvsupumnf, gvsupvmnf
      REAL(dp) :: r0
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fours=f_cos; fouru=f_sin; fourv=f_sin; fcomb=f_none
         mcol = 1; r0 = hs_i/2
      ELSE
         fours=f_sin; fouru=f_cos; fourv=f_cos; fcomb=f_sum
         mcol = 4; r0 = -hs_i/2
      END IF

      istat = LBOUND(jvsupsmnf,3)
      CALL ASSERT(istat.EQ.1,'LBOUND WRONG IN InitDisplacement')
      istat = UBOUND(jvsupsmnf,3)
      CALL ASSERT(istat.eq.ns,'UBOUND WRONG IN InitDisplacement')

      ALLOCATE (gvsupsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                gvsupumnf(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                gvsupvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in InitDisp')

!Apply scaling in real space (keep xc components unscaled)
      gvsupsmnf = jvsupsmnf(:,:,nsmin:nsmax)*col_scale(:,:,mcol,nsmin:nsmax)
      gvsupumnf = jvsupumnf(:,:,nsmin:nsmax)*col_scale(:,:,mcol+1,nsmin:nsmax)
      gvsupvmnf = jvsupvmnf(:,:,nsmin:nsmax)*col_scale(:,:,mcol+2,nsmin:nsmax)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
      NS1: IF (nsmin .EQ. 1) THEN
         gvsupsmnf(m0,:,1) = 0;  gvsupsmnf(m2:,:,1) = 0
         IF (.NOT.l_push_s) gvsupsmnf(m1,:,1) = 0

         gvsupumnf(m0,:,1) = 0;  gvsupumnf(m2:,:,1) = 0
         IF (.NOT.l_push_u) gvsupumnf(m1,:,1) = 0

         gvsupvmnf(m1:,:,1) = 0
         IF (.NOT.l_push_v) gvsupvmnf(m0,:,1) = 0

         IF (.NOT.l_natural) THEN
            gvsupsmnf(m1,:,1) = -r0*gvsupumnf(m1,:,1)
         END IF
       END IF NS1   

!m=0, n<0 set = 0 to avoid singular Hessian
       gvsupsmnf(m0,-ntor:-1,nsmin:nsmax) = 0       
       gvsupumnf(m0,-ntor:-1,nsmin:nsmax) = 0       
       gvsupvmnf(m0,-ntor:-1,nsmin:nsmax) = 0       

!Edge boundary conditions (B_tangential = 0 => vsups = 0)
       EDGE: IF (nsmax .EQ. ns) THEN
          gvsupsmnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumnf(:,:,ns) = 0
             gvsupvmnf(:,:,ns) = 0
          END IF
       END IF EDGE

!Calculate contravariant (SUP) velocities in real space on full mesh
       CALL toijsp_par(gvsupsmnf, jvsupsijf(:,:,nsmin:nsmax), fcomb,    &
                       fours)                         
       CALL toijsp_par(gvsupumnf, jvsupuijf(:,:,nsmin:nsmax), fcomb,    &
                       fouru)                        
       CALL toijsp_par(gvsupvmnf, jvsupvijf(:,:,nsmin:nsmax), fcomb,    &
                       fourv)

       DEALLOCATE (gvsupsmnf, gvsupumnf, gvsupvmnf,stat=istat)


      END SUBROUTINE InitDisplacement


      END MODULE siesta_displacement
