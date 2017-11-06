!>  \brief Module contained subroutines and functions updating MHD forces and Wmhd
!!  \author S. P. Hirshman
!!  \date July, 2014
      MODULE shared_functions
      USE v3_utilities, ONLY: assert
      USE shared_data
      USE descriptor_mod, ONLY: INHESSIAN, nprocs, iam, SIESTA_COMM
      USE island_params, ns=>ns_i, hs=>hs_i
      USE hessian, ONLY: apply_precond, l_Compute_Hessian
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            fsupsmncf, fsupumnsf, fsupvmnsf,            &
                            fsubsmnsf, fsubumncf, fsubvmncf,            &
                            fsupsmnsf, fsupumncf, fsupvmncf, wb, wp
      USE timer_mod, ONLY: time_init_state, time_funci
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR,       &
                               PARSOLVER
      USE mpi_inc
      USE GMRES_LIB, ONLY: Truncate
      
      IMPLICIT NONE

      INTEGER, PRIVATE :: nsmin, nsmax
      
      CONTAINS

!>  \brief Parallel/Serial routine to compute forces for perturbed state.
!!  \author S. P. Hirshman, R. Sanchez, S. K. Seal
      SUBROUTINE funct_island
!-----------------------------------------------
!     PARALLEL (in N rows, PARSOLVER=T)
!     SERIAL   (PARSOLVER=F, startglobrow=1, endglobrow=ns)
!-----------------------------------------------
      USE siesta_bfield, ONLY: update_bfield
      USE siesta_pressure, ONLY: update_pres
      USE siesta_force, ONLY: update_force
      USE siesta_displacement, ONLY: update_upperv
      USE siesta_init, ONLY: init_state
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff, skston, skstoff
      LOGICAL  :: ltype
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      IF (INHESSIAN) THEN
         in_hess_nfunct=in_hess_nfunct+1
      ELSE
         out_hess_nfunct=out_hess_nfunct+1
      END IF

      l_update_state = .FALSE.
      IF (l_init_state) THEN
        CALL second0(skston)
        CALL init_state(.FALSE.)
        CALL second0(skstoff)
        time_init_state=time_init_state+(skstoff-skston)
        l_init_state = .FALSE.
      END IF

      CALL update_upperv
      CALL update_bfield (.FALSE.)
      CALL update_pres
      CALL update_force

      CALL ASSERT(gamma.ne.one,'SIESTA REQUIRES gamma != 1')
      wtotal = wb + wp/(gamma-1)

      IF (wtotal0 == -1) wtotal0 = wtotal
      gc = gnorm_i*gc

!     CALLED FROM GMRES: ADD BACK INITIAL FORCE gc0 or STORE UNPRECONDITIONED FORCES
      IF (ALLOCATED(gc0) .AND. (l_getfsq .OR. l_conjgrad) .AND.         &
         l_linearize .AND. .NOT.l_Compute_Hessian) gc = gc+gc0

      CALL ASSERT(.not.(l_getfsq.and.inhessian),                        &
                  'l_getfsq must be set to FALSE in Hessian')

      IF (l_getfsq) THEN
! COMPUTE PRECONDITIONED, VOLUME-AVERAGED FORCE fsq_total 
         CALL toupper_forces
         CALL GetFSQ(fsq_total)

!COMPUTE UN-PRECONDITIONED, VOLUME-AVERAGED FORCE fsq_total1
         IF (ANY(col_scale .NE. one)) THEN
            fsubsmncf(:,:,nsmin:nsmax)=fsubsmncf(:,:,nsmin:nsmax)/col_scale(:,:,1,nsmin:nsmax)
            fsubumnsf(:,:,nsmin:nsmax)=fsubumnsf(:,:,nsmin:nsmax)/col_scale(:,:,2,nsmin:nsmax)
            fsubvmnsf(:,:,nsmin:nsmax)=fsubvmnsf(:,:,nsmin:nsmax)/col_scale(:,:,3,nsmin:nsmax)
            IF (lasym) THEN
            fsubsmnsf(:,:,nsmin:nsmax)=fsubsmnsf(:,:,nsmin:nsmax)/col_scale(:,:,4,nsmin:nsmax)
            fsubumncf(:,:,nsmin:nsmax)=fsubumncf(:,:,nsmin:nsmax)/col_scale(:,:,5,nsmin:nsmax)
            fsubvmncf(:,:,nsmin:nsmax)=fsubvmncf(:,:,nsmin:nsmax)/col_scale(:,:,6,nsmin:nsmax)
            END IF

            CALL toupper_forces
            CALL GetFSQ(fsq_total1)
            
            fsubsmncf(:,:,nsmin:nsmax)=fsubsmncf(:,:,nsmin:nsmax)*col_scale(:,:,1,nsmin:nsmax)
            fsubumnsf(:,:,nsmin:nsmax)=fsubumnsf(:,:,nsmin:nsmax)*col_scale(:,:,2,nsmin:nsmax)
            fsubvmnsf(:,:,nsmin:nsmax)=fsubvmnsf(:,:,nsmin:nsmax)*col_scale(:,:,3,nsmin:nsmax)
            IF (lasym) THEN
            fsubsmnsf(:,:,nsmin:nsmax)=fsubsmnsf(:,:,nsmin:nsmax)*col_scale(:,:,4,nsmin:nsmax)
            fsubumncf(:,:,nsmin:nsmax)=fsubumncf(:,:,nsmin:nsmax)*col_scale(:,:,5,nsmin:nsmax)
            fsubvmncf(:,:,nsmin:nsmax)=fsubvmncf(:,:,nsmin:nsmax)*col_scale(:,:,6,nsmin:nsmax)
            END IF
         ELSE
            fsq_total1 = fsq_total
         END IF
      END IF

!      IF (l_Gmres .AND. l_ApplyPrecon) STOP 'l_ApplyPrecon must be FALSE in funct_island'
!
!     GMRES handles preconditioning itself: do NOT apply
!     Do not precondition when Hessian is being computed
      IF (l_ApplyPrecon) CALL apply_precond(gc)

      CALL second0(toff)
      time_funci = time_funci+(toff-ton)

      END SUBROUTINE funct_island


      SUBROUTINE GetFSQ(fsqout)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(OUT) :: fsqout
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: stmp(3), rtmp(3)
      INTEGER  :: js
!-----------------------------------------------
      stmp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)**2)
      stmp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)**2)
      stmp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)**2)
      IF (lasym) THEN
         stmp(1) = stmp(1) + SUM(fsubsmnsf(:,:,nsmin:nsmax)**2)
         stmp(2) = stmp(2) + SUM(fsubumncf(:,:,nsmin:nsmax)**2)
         stmp(3) = stmp(3) + SUM(fsubvmncf(:,:,nsmin:nsmax)**2)
      END IF
         
#if defined(MPI_OPT)         
      IF (PARSOLVER) THEN
         CALL MPI_ALLREDUCE(stmp,rtmp,3,MPI_REAL8,MPI_SUM,              &
                            SIESTA_COMM,MPI_ERR)
         stmp = rtmp
      END IF
#endif         
      fsqvs=hs*stmp(1); fsqvu=hs*stmp(2); fsqvv=hs*stmp(3)

!VOLUME-AVERAGE |F|**2
      stmp(1) = 0
      DO js = nsmin,nsmax
         stmp(1) = stmp(1) + vp_f(js)*SUM(                              &
                   fsupsmncf(:,:,js)*fsubsmncf(:,:,js)                  & 
                 + fsupumnsf(:,:,js)*fsubumnsf(:,:,js)                  &
                 + fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
         IF (lasym) THEN
         stmp(1) = stmp(1) + vp_f(js)*SUM(                              &
                   fsupsmnsf(:,:,js)*fsubsmnsf(:,:,js)                  &
                 + fsupumncf(:,:,js)*fsubumncf(:,:,js)                  &
                 + fsupvmncf(:,:,js)*fsubvmncf(:,:,js))
         END IF
      END DO
      stmp(2) = SUM(vp_f(nsmin:nsmax))
#if defined(MPI_OPT)         
      IF (PARSOLVER) THEN
         CALL MPI_ALLREDUCE(stmp,rtmp,2,MPI_REAL8,MPI_SUM,              &
                            SIESTA_COMM,MPI_ERR)
         stmp = rtmp
      END IF
#endif
      fsqout = stmp(1)/stmp(2)
      
      
      END SUBROUTINE GetFSQ


      SUBROUTINE toupper_forces
      USE island_params, ONLY: ntheta=>nu_i, nzeta=>nv_i,               &
                         mpol=>mpol_i, ntor=>ntor_i
      USE fourier, ONLY: toijsp_par, tomnsp_par, f_none, f_cos, f_sin,  &
                         f_sum
      USE metrics, ONLY: toupper
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                        &
                                    fsubsijf, fsubuijf, fsubvijf,       &
                                    fsupsijf, fsupuijf, fsupvijf
!-----------------------------------------------
      INTEGER :: nsmin, nsmax, istat
!-----------------------------------------------
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)

      ALLOCATE(fsubsijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsubuijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsubvijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupsijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupuijf(ntheta,nzeta,nsmin:nsmax),                      &
               fsupvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,' Allocation error in toupper_forces')

      CALL toijsp_par(fsubsmncf(:,:,nsmin:nsmax), fsubsijf, f_none, f_cos)
      CALL toijsp_par(fsubumnsf(:,:,nsmin:nsmax), fsubuijf, f_none, f_sin)
      CALL toijsp_par(fsubvmnsf(:,:,nsmin:nsmax), fsubvijf, f_none, f_sin)
 
      IF (lasym) THEN

      CALL toijsp_par(fsubsmnsf(:,:,nsmin:nsmax), fsubsijf, f_sum, f_sin)
      CALL toijsp_par(fsubumncf(:,:,nsmin:nsmax), fsubuijf, f_sum, f_cos)
      CALL toijsp_par(fsubvmncf(:,:,nsmin:nsmax), fsubvijf, f_sum, f_cos)
 
      END IF

!      
!     COMPUTE CONTRAVARIANT COMPONENTS FOR CHECKING OUTPUT
!
      CALL toupper(fsubsijf, fsubuijf, fsubvijf, fsupsijf, fsupuijf, fsupvijf) 
      DEALLOCATE(fsubsijf, fsubuijf, fsubvijf)

      CALL tomnsp_par(fsupsijf, fsupsmncf(:,:,nsmin:nsmax), f_cos)
      CALL tomnsp_par(fsupuijf, fsupumnsf(:,:,nsmin:nsmax), f_sin)
      CALL tomnsp_par(fsupvijf, fsupvmnsf(:,:,nsmin:nsmax), f_sin)

      IF (lasym) THEN

      CALL tomnsp_par(fsupsijf, fsupsmnsf(:,:,nsmin:nsmax), f_sin)
      CALL tomnsp_par(fsupuijf, fsupumncf(:,:,nsmin:nsmax), f_cos)
      CALL tomnsp_par(fsupvijf, fsupvmncf(:,:,nsmin:nsmax), f_cos)

      END IF

      DEALLOCATE(fsupsijf, fsupuijf, fsupvijf, stat=istat)

      END SUBROUTINE toupper_forces

!>   \brief Returns perturbed MHD+KINETIC energy
!!   \author S. P. Hirshman, R. Sanchez
      FUNCTION getwmhd (p)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(in) :: p(neqs)           !<displacement vector
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: getwmhd
!-----------------------------------------------
        
      xc = p
      l_init_state = .TRUE.
      l_linearize  = .FALSE.
      l_getwmhd    = .TRUE.
      l_getfsq     = .FALSE.

      CALL funct_island

      l_getwmhd    = .FALSE.

      getwmhd = wtotal

      RETURN

      END FUNCTION getwmhd 

!>    \brief Performs line search along xc vector for minimum fsq_total1
!!    \author S. P. Hirshman (2016)
      SUBROUTINE LineSearch(xcmin, fsq_min)
      USE perturbation, ONLY: ftol
      IMPLICIT NONE
      REAL(dp), INTENT(INOUT) :: xcmin(neqs)
      REAL(dp), INTENT(INOUT) :: fsq_min
      REAL(dp), PARAMETER :: facmin = 2
      INTEGER  :: iter, j

      IF (iam .EQ. 0 .and. lverbose) PRINT 90
 90   FORMAT(1x,'LINE SEARCH - SCAN ||X|| FOR MIN FSQ_NL',/,'-------------',   &
           /,1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|')

      xc = xcmin
      iter = 0

      DO j = 1, 100
         CALL funct_island
         IF (j .EQ. 1) THEN
            fsq_min = fsq_total1
            xc = facmin*xcmin
            CYCLE
         ELSE IF (j.NE.4 .AND. fsq_total1.LT.fsq_min) THEN
            xcmin = xc
            IF (fsq_total1 .GT. 0.98_dp*fsq_min) iter=iter+1
            IF (fsq_total1 .LT. 0.85_dp*fsq_min) iter = 0
            fsq_min = fsq_total1
         ELSE IF (j .GT. 4) THEN
            iter = iter+1
         END IF
         IF (iam .EQ. 0 .and. lverbose)                                        &
            PRINT 100, j-1, fsq_total1, SQRT(SUM(xc*xc)), MAXVAL(ABS(xc))
         IF (iter.GT.2 .OR. fsq_total1.LE.ftol) EXIT
         xc = xc/SQRT2
      END DO

 100  FORMAT(i5,3(3x,1pe12.3))
      END SUBROUTINE LineSearch


      END MODULE shared_functions
