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
      USE island_params, ONLY: hs_i, ns=>ns_i
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants        
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE timer_mod, ONLY: time_update_upperv
      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
      INTEGER, PRIVATE ::  nsmin, nsmax

      CONTAINS
!>    \brief Parallel/Serial subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv
      USE shared_data, ONLY: lasym, gc_sup, xc, col_scale
      USE siesta_state, ONLY: clear_field_perts
      USE quantities, ONLY: gvsupsmncf => fsupsmncf, gvsupumnsf => fsupumnsf,  &
             gvsupvmnsf => fsupvmnsf, gvsupsmnsf => fsupsmnsf,          &
             gvsupumncf => fsupumncf, gvsupvmncf => fsupvmncf
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

!Use components of gc_sup as scratch array
      CALL ScaleDisplacement(gc_sup, xc, col_scale)

      CALL InitDisplacement(gvsupsmncf, gvsupumnsf, gvsupvmnsf, isym)
      IF (lasym)                                                        &
      CALL InitDisplacement(gvsupsmnsf, gvsupumncf, gvsupvmncf, iasym)

      CALL second0(toff)
      time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv
      
 
      SUBROUTINE InitDisplacement(gvsupsmnf, gvsupumnf, gvsupvmnf, ipar)
      USE fourier, ONLY: toijsp_par, f_cos, f_sin, f_none, f_sum
      USE quantities, ONLY: jvsupsijf, jvsupuijf, jvsupvijf, mpol, ntor
      USE shared_data, ONLY: l_push_s, l_push_u, l_push_v, l_pedge,     &
                             l_push_edge, col_scale
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns)  ::                     &
                gvsupsmnf, gvsupumnf, gvsupvmnf
      INTEGER, INTENT(IN)  :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1, m2=2
      INTEGER  :: istat, m, n, mcol,                                    &
                  fours, fouru, fourv, fcomb 
!-----------------------------------------------
      IF (ipar .EQ. isym) THEN
         fours=f_cos; fouru=f_sin; fourv=f_sin; fcomb=f_none
         gvsupumnf(m0,0,nsmin:nsmax) = 0; gvsupvmnf(m0,0,nsmin:nsmax) = 0
      ELSE
         fours=f_sin; fouru=f_cos; fourv=f_cos; fcomb=f_sum
         gvsupsmnf(m0,0,nsmin:nsmax) = 0
      END IF
      
      istat = LBOUND(gvsupsmnf,3)
      CALL ASSERT(istat.EQ.1,'LBOUND WRONG IN InitDisplacement')
      istat = UBOUND(gvsupsmnf,3)
      CALL ASSERT(istat.EQ.ns,'UBOUND WRONG IN InitDisplacement')

!m=0, n<0 set = 0 to avoid singular Hessian
      gvsupsmnf(m0,-ntor:-1,nsmin:nsmax) = 0       
      gvsupumnf(m0,-ntor:-1,nsmin:nsmax) = 0       
      gvsupvmnf(m0,-ntor:-1,nsmin:nsmax) = 0       

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
      NS1: IF (nsmin .EQ. 1) THEN
         gvsupsmnf(m0,:,1) = 0;  gvsupsmnf(m2:,:,1) = 0
         IF (.NOT.l_push_s) gvsupsmnf(m1,:,1) = 0

         gvsupumnf(m0,:,1) = 0;  gvsupumnf(m2:,:,1) = 0
         IF (.NOT.l_push_u) gvsupumnf(m1,:,1) = 0

         gvsupvmnf(m1:,:,1) = 0
         IF (.NOT.l_push_v) gvsupvmnf(m0,:,1) = 0
       END IF NS1   

!Edge boundary conditions (B_tangential = 0 => vsups = 0; J^s = 0 => F_v = 0)
!Note that F_v and F_u are both ~J^s at edge, so only one can be varied to
!prevent a dependent row in the Hessian matrix
       EDGE: IF (nsmax .EQ. ns) THEN
          gvsupsmnf(:,:,ns) = 0
          IF (l_pedge) gvsupumnf(:,:,ns) = 0
          IF (.NOT.l_push_edge) THEN
             gvsupumnf(:,:,ns) = 0
             gvsupvmnf(:,:,ns) = 0
          END IF
       END IF EDGE

!Calculate contravariant (SUP) velocities in real space on full mesh
       CALL toijsp_par(gvsupsmnf(:,:,nsmin:nsmax),                      &
                       jvsupsijf(:,:,nsmin:nsmax), fcomb, fours)                         
       CALL toijsp_par(gvsupumnf(:,:,nsmin:nsmax),                      &
                       jvsupuijf(:,:,nsmin:nsmax), fcomb, fouru)                        
       CALL toijsp_par(gvsupvmnf(:,:,nsmin:nsmax),                      &
                       jvsupvijf(:,:,nsmin:nsmax), fcomb, fourv)        

      END SUBROUTINE InitDisplacement
      

      SUBROUTINE ScaleDisplacement(xc_scratch, xc, colscale)
      USE hessian, ONLY: apply_colscale
      USE shared_data, ONLY: mblk_size
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN), DIMENSION(mblk_size,ns)  :: xc, colscale
      REAL(dp), INTENT(OUT)  :: xc_scratch(mblk_size, ns)
      
      xc_scratch(:, nsmin:nsmax) = xc(:, nsmin:nsmax)
      CALL Apply_ColScale(xc_scratch, colscale, nsmin, nsmax)

      END SUBROUTINE ScaleDisplacement

      END MODULE siesta_displacement
