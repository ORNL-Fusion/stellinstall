
!>  \brief Module contained subroutines and functions computing MHD forces and energy
!!  \author S. P. Hirshman
!!  \date July, 2014
      MODULE shared_functions
      USE shared_data
      USE descriptor_mod, ONLY: INHESSIAN, nprocs, iam, SIESTA_COMM,     &
                                in_hess_nfunct, out_hess_nfunct,         &
                                SIESTA_COMM
      USE island_params, ns=>ns_i, hs=>hs_i
      USE hessian, ONLY: apply_precond, l_Compute_Hessian
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,             &
                            fsupsmncf, fsupumnsf, fsupvmnsf,             &
                            fsubsmnsf, fsubumncf, fsubvmncf,             &
                            fsupsmnsf, fsupumncf, fsupvmncf,             &
                            wb, wp
      USE timer_mod, ONLY: time_init_state, time_funci
      USE siesta_init, ONLY: init_state
      USE nscalingtools, ONLY: rcountsNS, dispNS, PARFUNCTISL,           &
                               startglobrow, endglobrow, MPI_ERR,        &
                               PARSOLVER
      CONTAINS

!>  \brief Parallel routine to compute MHD forces (both linear AND nonlinear)
!!  depending on l_linearize = T/F
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
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff, skston, skstoff, fs, fa
      INTEGER  :: nsmin, nsmax, js
      REAL(dp) :: stmp(3), rtmp(3), sendbuf(ns), recbuf(ns)
      LOGICAL  :: ltype
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      IF (INHESSIAN) THEN
         in_hess_nfunct=in_hess_nfunct+1
      ELSE
         out_hess_nfunct=out_hess_nfunct+1
      END IF

      IF (l_init_state) THEN
        CALL second0(skston)
        CALL init_state(.FALSE.)
        CALL second0(skstoff)
        time_init_state=time_init_state+(skstoff-skston)
        l_init_state = .FALSE.
      END IF

      CALL update_upperv
      CALL update_bfield(.FALSE.)
      CALL update_pres
      CALL update_force

      IF (gamma .eq. 1._dp) CALL ErrorAbort('SIESTA REQUIRES gamma != 1')
      wtotal = wb + wp/(gamma-1)

      IF (wtotal0 == -1) wtotal0 = wtotal
      gc = gnorm_i*gc

!     CALLED FROM GMRES: ADD BACK INITIAL FORCE gc0 or STORE UNPRECONDITIONED FORCES
      IF (ALLOCATED(gc0) .AND. l_linearize .AND.                         &
         .NOT.l_Compute_Hessian .AND. (l_getfsq .OR. l_conjgrad)) THEN
         gc = gc+gc0
      END IF

      IF (l_getfsq .AND. INHESSIAN) CALL ErrorAbort('l_getfsq must be set to FALSE in Hessian')

!     COMPUTE UN-PRECONDITIONED FORCE RESIDUES
      IF (l_getfsq) THEN
         CALL toupper_forces_par
!     VOLUME-AVERAGE |F|**2
         stmp(1) = 0
         DO js = nsmin,nsmax
            stmp(1) = stmp(1) + vp_f(js)*SUM(                            &
                      fsupsmncf(:,:,js)*fsubsmncf(:,:,js)                & 
                    + fsupumnsf(:,:,js)*fsubumnsf(:,:,js)                &
                    + fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
            IF (lasym) THEN
            stmp(1) = stmp(1) + vp_f(js)*SUM(                            &
                      fsupsmnsf(:,:,js)*fsubsmnsf(:,:,js)                &
                    + fsupumncf(:,:,js)*fsubumncf(:,:,js)                &
                    + fsupvmncf(:,:,js)*fsubvmncf(:,:,js))
            END IF

         END DO
         stmp(2) = SUM(vp_f(nsmin:nsmax))
         IF (PARSOLVER) THEN
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(stmp,rtmp,2,MPI_REAL8,MPI_SUM,SIESTA_COMM,MPI_ERR)
         stmp = rtmp
#endif
         END IF
         fsq_total = stmp(1)/stmp(2)
         IF (ANY(col_scale .NE. one)) THEN
            fsubsmncf(:,:,nsmin:nsmax)=fsubsmncf(:,:,nsmin:nsmax)/col_scale(:,:,1,nsmin:nsmax)
            fsubumnsf(:,:,nsmin:nsmax)=fsubumnsf(:,:,nsmin:nsmax)/col_scale(:,:,2,nsmin:nsmax)
            fsubvmnsf(:,:,nsmin:nsmax)=fsubvmnsf(:,:,nsmin:nsmax)/col_scale(:,:,3,nsmin:nsmax)
            IF (lasym) THEN
            fsubsmnsf(:,:,nsmin:nsmax)=fsubsmnsf(:,:,nsmin:nsmax)/col_scale(:,:,4,nsmin:nsmax)
            fsubumncf(:,:,nsmin:nsmax)=fsubumncf(:,:,nsmin:nsmax)/col_scale(:,:,5,nsmin:nsmax)
            fsubvmncf(:,:,nsmin:nsmax)=fsubvmncf(:,:,nsmin:nsmax)/col_scale(:,:,6,nsmin:nsmax)
            END IF

            CALL toupper_forces_par
!    
!     THIS WAY GIVES GREATER PROCESSOR # INDEPENDENCE OF CONVERGENGE SEQUENCE   
!
            DO js = nsmin,nsmax
               sendbuf(js) = vp_f(js)*SUM(                               &
                         fsupsmncf(:,:,js)*fsubsmncf(:,:,js)             & 
                       + fsupumnsf(:,:,js)*fsubumnsf(:,:,js)             &
                       + fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
            END DO
            IF (lasym) THEN
            DO js = nsmin,nsmax
               sendbuf(js) =  sendbuf(js) + vp_f(js)*SUM(                &
                         fsupsmnsf(:,:,js)*fsubsmnsf(:,:,js)             & 
                       + fsupumncf(:,:,js)*fsubumncf(:,:,js)             &
                       + fsupvmncf(:,:,js)*fsubvmncf(:,:,js))
            END DO
            END IF
          
            IF (PARSOLVER) THEN
#if defined(MPI_OPT)
            CALL MPI_ALLGATHERV(sendbuf(nsmin), nsmax-nsmin+1, MPI_REAL8,&
                                recbuf, rcountsNS, dispNS, MPI_REAL8,    &
                                SIESTA_COMM, MPI_ERR)
            fsq_total1 = SUM(recbuf)/rtmp(2)
#endif
            ELSE
            fsq_total1 = SUM(sendbuf)/stmp(2)
            END IF

            stmp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)**2)
            stmp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)**2)
            stmp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)**2)
            IF (lasym) THEN
            stmp(1) = stmp(1) + SUM(fsubsmnsf(:,:,nsmin:nsmax)**2)
            stmp(2) = stmp(2) + SUM(fsubumncf(:,:,nsmin:nsmax)**2)
            stmp(3) = stmp(3) + SUM(fsubvmncf(:,:,nsmin:nsmax)**2)
            END IF

            IF (PARSOLVER) THEN
#if defined(MPI_OPT)
            CALL MPI_ALLREDUCE(stmp,rtmp,3,MPI_REAL8,MPI_SUM,SIESTA_COMM,MPI_ERR)
            stmp = rtmp
#endif
            END IF
            fsqvs=hs*stmp(1); fsqvu=hs*stmp(2); fsqvv=hs*stmp(3)

            fsubsmncf(:,:,nsmin:nsmax)=fsubsmncf(:,:,nsmin:nsmax)*col_scale(:,:,1,nsmin:nsmax)
            fsubumnsf(:,:,nsmin:nsmax)=fsubumnsf(:,:,nsmin:nsmax)*col_scale(:,:,2,nsmin:nsmax)
            fsubvmnsf(:,:,nsmin:nsmax)=fsubvmnsf(:,:,nsmin:nsmax)*col_scale(:,:,3,nsmin:nsmax)
            IF (lasym) THEN
            fsubsmnsf(:,:,nsmin:nsmax)=fsubsmnsf(:,:,nsmin:nsmax)*col_scale(:,:,4,nsmin:nsmax)
            fsubumncf(:,:,nsmin:nsmax)=fsubumncf(:,:,nsmin:nsmax)*col_scale(:,:,5,nsmin:nsmax)
            fsubvmncf(:,:,nsmin:nsmax)=fsubvmncf(:,:,nsmin:nsmax)*col_scale(:,:,6,nsmin:nsmax)
            END IF
         ELSE         
            stmp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)**2)
            stmp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)**2)
            stmp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)**2)
            IF (lasym) THEN
            stmp(1) = stmp(1) + SUM(fsubsmnsf(:,:,nsmin:nsmax)**2)
            stmp(2) = stmp(2) + SUM(fsubumncf(:,:,nsmin:nsmax)**2)
            stmp(3) = stmp(3) + SUM(fsubvmncf(:,:,nsmin:nsmax)**2)
            END IF

            IF (PARSOLVER) THEN
#if defined(MPI_OPT)
            CALL MPI_ALLREDUCE(stmp,rtmp,3,MPI_REAL8,MPI_SUM,SIESTA_COMM,MPI_ERR)
            stmp = rtmp
#endif
            END IF
            fsqvs=hs*stmp(1); fsqvu=hs*stmp(2); fsqvv=hs*stmp(3)
            fsq_total1 = fsq_total
         END IF
      END IF

!      IF (l_Gmres .AND. l_ApplyPrecon) CALL ErrorAbort('l_ApplyPrecon must be FALSE in funct_island')
!
!     GMRES handles preconditioning itself: do NOT apply
!     Do not precondition when Hessian is being computed
      IF (l_ApplyPrecon) CALL apply_precond(gc)

      CALL second0(toff)
      time_funci = time_funci+(toff-ton)

      END SUBROUTINE funct_island

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


!>   \brief Returns perturbed MHD+KINETIC energy
!!   \author S. P. Hirshman, R. Sanchez
      FUNCTION getwmhd (p)
      USE GMRES_LIB, ONLY: Truncate
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN) :: p(neqs)           !<displacement vector
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

      CALL Truncate(wtotal, 12)
      getwmhd = wtotal

      END FUNCTION getwmhd 

      END MODULE shared_functions
