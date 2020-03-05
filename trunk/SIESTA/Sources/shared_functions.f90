!>  \brief Module contained subroutines and functions updating MHD forces and Wmhd
!!  \author S. P. Hirshman
!!  \date July, 2014
      MODULE shared_functions
      USE v3_utilities, ONLY: assert
      USE shared_data
      USE utilities
      USE descriptor_mod, ONLY: INHESSIAN, nprocs, iam, SIESTA_COMM
      USE island_params, ns=>ns_i, hs=>hs_i, mpol=>mpol_i, ntor=>ntor_i
      USE hessian, ONLY: apply_precond, l_Compute_Hessian, apply_colscale
      USE timer_mod, ONLY: time_init_state, time_funci
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR,       &
                               PARSOLVER
      USE mpi_inc
      USE GMRES_LIB, ONLY: Truncate
      
      IMPLICIT NONE

      INTEGER, PRIVATE :: nsmin, nsmax
      
      CONTAINS

!>  \brief Parallel/Serial routine to compute forces for perturbed state.
!!  \author S. K. Seal and S. P. Hirshman
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
      USE quantities, ONLY: wb, wp
      USE siesta_state, ONLY: Clear_Field_Perts
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff, skston, skstoff
      LOGICAL  :: ltype
!-----------------------------------------------
      CALL second0(ton)
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(endglobrow, ns)

      IF (INHESSIAN) THEN
         in_hess_nfunct = in_hess_nfunct + 1
      ELSE
         out_hess_nfunct = out_hess_nfunct + 1
      END IF

      l_update_state = .FALSE.
      IF (l_init_state) THEN
         CALL second0(skston)
         CALL init_state(.FALSE.)
         CALL second0(skstoff)
         time_init_state = time_init_state + (skstoff - skston)
         l_init_state = .FALSE.
      END IF

      IF (ANY(xc .ne. 0) .or. ALLOCATED(buv_res)) THEN
         CALL update_upperv
         CALL update_bfield(.FALSE.)
         CALL update_pres
      ELSE
         CALL Clear_Field_Perts
      END IF
      CALL update_force

      CALL ASSERT(gamma.NE.one,'SIESTA REQUIRES gamma != 1')
      wtotal = wb + wp/(gamma - 1)
      IF (wtotal0 .eq. -1) THEN
         wtotal0 = wtotal
      END IF

      gc = gnorm_i*gc

!     CALLED FROM GMRES: ADD BACK INITIAL FORCE gc0 or STORE UNPRECONDITIONED FORCES
      IF (ALLOCATED(gc0) .AND. (l_getfsq .OR. l_conjgrad) .AND.         &
         l_linearize .AND. .NOT.l_Compute_Hessian) gc = gc+gc0

      CALL ASSERT(.NOT.(l_getfsq.AND.inhessian),                        &
                  'l_getfsq must be set to FALSE in Hessian')

      IF (l_getfsq) THEN
! COMPUTE PRECONDITIONED, VOLUME-AVERAGED FORCE fsq_total 
         CALL GetFSQ(fsq_total)

!COMPUTE UN-PRECONDITIONED, VOLUME-AVERAGED FORCE fsq_total1
         IF (ANY(col_scale .NE. one)) THEN
            col_scale(:,:,:,nsmin:nsmax) = one/col_scale(:,:,:,nsmin:nsmax)
            CALL Apply_ColScale(gc, col_scale, nsmin, nsmax)

            CALL GetFSQ(fsq_total1)

            col_scale(:,:,:,nsmin:nsmax) = one/col_scale(:,:,:,nsmin:nsmax)            
            CALL Apply_ColScale(gc, col_scale, nsmin, nsmax)
         ELSE
            fsq_total1 = fsq_total
         END IF
      END IF

!
!     GMRES handles preconditioning itself: do NOT apply
!     Do not precondition when Hessian is being computed
      IF (l_ApplyPrecon) CALL apply_precond(gc)

      CALL second0(toff)
      time_funci = time_funci+(toff-ton)

      END SUBROUTINE funct_island


      SUBROUTINE GetFSQ(fsqout)
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            fsupsmncf, fsupumnsf, fsupvmnsf,            &
                            fsubsmnsf, fsubumncf, fsubvmncf,            &
                            fsupsmnsf, fsupumncf, fsupvmncf,            &
                            toupper_forces
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(OUT) :: fsqout
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: temp(3)
      INTEGER  :: js
!-----------------------------------------------
      temp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)**2)
      temp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)**2)
      temp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)**2)
      IF (lasym) THEN
         temp(1) = temp(1) + SUM(fsubsmnsf(:,:,nsmin:nsmax)**2)
         temp(2) = temp(2) + SUM(fsubumncf(:,:,nsmin:nsmax)**2)
         temp(3) = temp(3) + SUM(fsubvmncf(:,:,nsmin:nsmax)**2)
      END IF
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,3,MPI_REAL8,MPI_SUM,             &
                            SIESTA_COMM,MPI_ERR)
      END IF
#endif
      fsqvs = hs*temp(1)
      fsqvu = hs*temp(2)
      fsqvv = hs*temp(3)

      CALL toupper_forces

!VOLUME-AVERAGE |F|**2
      temp(1) = 0
      DO js = nsmin, nsmax
         temp(1) = temp(1)                                                     &
                 + vp_f(js)*SUM(fsupsmncf(:,:,js)*fsubsmncf(:,:,js) +          &
                                fsupumnsf(:,:,js)*fsubumnsf(:,:,js) +          &
                                fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
         IF (lasym) THEN
            temp(1) = temp(1)                                                  &
                    + vp_f(js)*SUM(fsupsmnsf(:,:,js)*fsubsmnsf(:,:,js) +       &
                                   fsupumncf(:,:,js)*fsubumncf(:,:,js) +       &
                                   fsupvmncf(:,:,js)*fsubvmncf(:,:,js))
         END IF
      END DO
#if defined(MPI_OPT)         
      IF (PARSOLVER) THEN
!  FIXME: All reduce is not deterministic. This causes a divergent run sequence.
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,1,MPI_REAL8,MPI_SUM,             &
                            SIESTA_COMM,MPI_ERR)
      END IF
#endif

      fsqout = temp(1)/SUM(vp_f)

      END SUBROUTINE GetFSQ


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

      END FUNCTION getwmhd 
      

!>    \brief Performs line search along xc vector for minimum fsq_total1
!!    \author S. P. Hirshman (2016)
      SUBROUTINE LineSearch(xcmin, fsq_min)
      USE siesta_namelist, ONLY: ftol
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(INOUT) :: xcmin(neqs)
      REAL(dp), INTENT(INOUT) :: fsq_min
!-----------------------------------------------
      REAL(dp) :: facmin
      INTEGER  :: iter, j
      LOGICAL  :: lwrite
!-----------------------------------------------
      facmin = 3
      l_init_state = .FALSE.
      l_getfsq = .TRUE.
      
      lwrite = (iam .EQ. 0 .and. lverbose) 
      IF (lwrite) PRINT 90
 90   FORMAT(/,1x,'LINE SEARCH - SCAN ||X|| FOR MIN FSQ_NL',/,1x,15('-'), &
           /,1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|',10x,'FAC')

      xc = facmin*xcmin
      iter = 0
      fsq_min = HUGE(fsq_min)

      DO j = 1, 100
         CALL funct_island
         IF (fsq_total1 .lt. fsq_min) THEN
            xcmin = xc
!  FIXME: Remove one.
            IF (fsq_total1 .GT. 0.98_dp*fsq_min) THEN
               iter = iter + 1
            END IF
            IF (fsq_total1 .LT. 0.85_dp*fsq_min) THEN
               iter = 0
            END IF
            fsq_min = fsq_total1
         ELSE IF (j .gt. 4) THEN
            iter = iter + 1
         END IF
         IF (iam .eq. 0 .and. lverbose) THEN
            WRITE (*,1000) j, fsq_total1, SQRT(SUM(xc*xc)), MAXVAL(ABS(xc)),   &
                           facmin
            facmin = facmin/SQRT2
         END IF
         IF (iter .gt. 2 .OR. fsq_total1 .le. ftol) EXIT
         xc = xc/SQRT2
      END DO

1000  FORMAT(i5,4(3x,1pe12.3))
 
      END SUBROUTINE LineSearch

!>    \brief Initializes pointers for xc (displacement variables) and gc (force) arrays
      SUBROUTINE init_ptrs (xtarget, ptr1, ptr2, ptr3, ptr4, ptr5, ptr6)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), TARGET, INTENT(IN)        :: xtarget(0:mpol,-ntor:ntor,ndims,ns) !<input 3D, ndims-component arrays
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr1, ptr4  !<1st (s) array subsection (lasym=F,T)
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr2, ptr5  !<2nd (u) array subsection
      REAL(dp), POINTER, DIMENSION(:,:,:) :: ptr3, ptr6  !<3rd (v) array subsection

!-----------------------------------------------
      ptr1 => xtarget(:,:,1,:)
      ptr2 => xtarget(:,:,2,:)
      ptr3 => xtarget(:,:,3,:)
      IF (ndims .EQ. 6) THEN
         ptr4 => xtarget(:,:,4,:)
         ptr5 => xtarget(:,:,5,:)
         ptr6 => xtarget(:,:,6,:)
      END IF

      END SUBROUTINE init_ptrs

      END MODULE shared_functions
