      SUBROUTINE add_resistive_E
      USE stel_kinds
      USE stel_constants, ONLY: zero
      USE shared_data, ONLY: fsq_total1, xc, l_getfsq, l_init_state, lverbose
      USE shared_functions, ONLY: funct_island
      USE perturbation, ONLY: ftol, eta_factor
      USE descriptor_mod, ONLY: iam
      USE siesta_bfield, ONLY: update_bfield
      USE siesta_init, ONLY: init_state
      USE nscalingtools, ONLY: startglobrow, endglobrow, PARFUNCTISL
      USE siesta_state, ONLY: update_state_par, update_state
      USE siesta_displacement, ONLY: update_upperv
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER   :: nits = 10
      REAL(dp)             :: save_eta, fmhd
      INTEGER              :: itime, j
      LOGICAL, PARAMETER   :: lcurrent_only=.TRUE.
!-----------------------------------------------
      IF (fsq_total1.LE.ftol) RETURN

      fmhd = fsq_total1
!
!     Update resistive contribution to B with new currents
!       
!     No perturbation, vsupX = 0
      xc = 0

      CALL init_state(.FALSE.)
      CALL update_upperv

      save_eta = eta_factor
      eta_factor = eta_factor/MAX(1,nits)

!     Diffuse B-field (eta*j) from last update_state call
      DO itime=1,nits

         CALL init_state(lcurrent_only)                                    !Use new B-field to get KsubXij currents
         CALL update_bfield(.TRUE.)                                        !Resistive update to B-field
         CALL update_state_par(.FALSE., zero, zero)                        !Updates B-field with resistive piece

      END DO

      eta_factor = save_eta
      
      l_getfsq = .TRUE.
      l_init_state = .TRUE.
      CALL funct_island

      IF (iam .EQ. 0) THEN
         DO j = 6, 33, 27
         IF (.not.lverbose .and. j .eq. 6) CYCLE
            WRITE (j,'(/,a,i3)') ' UPDATING RESISTIVE E-FIELD: ITERATIONS=',nits
            WRITE (j, '(a,1p2e12.3)')                                    &
          ' JUMP IN FSQ DUE TO RESISTIVE DIFFUSION: ', fmhd, fsq_total1
         END DO
      END IF

!SPH101116: limit jump next time
      IF (fsq_total1 .GT. 4*fmhd) THEN
         eta_factor = eta_factor*(4*fmhd/fsq_total1)
      ELSE IF (fsq_total1 .LT. 1.1_dp*fmhd) THEN
         eta_factor = eta_factor*(1.1_dp*fmhd/fsq_total1)
      END IF
      

!SPH101116 NEED THIS FOR GMRES TO INITIALIZE CORRECTLY
      fsq_total1 = fmhd

      END SUBROUTINE add_resistive_E
