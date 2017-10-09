!>  \brief Module contained subroutines for updating (from t to t+delta_t) the magnetic field and pressure as part of the SIESTA project.
!!  Stores updated values of JPMN?H, JBSUPXMN?H in Quantities
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Sep 15, 2006
      MODULE siesta_state
#define DUMP_STATE
      USE stel_kinds
      USE quantities
      USE diagnostics_mod
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio_s, bsbu_ratio_a
      USE timer_mod, ONLY: time_update_state
      USE siesta_init, ONLY: init_state
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE mpi_inc

      IMPLICIT NONE
      LOGICAL, PUBLIC :: lfirst=.TRUE.      
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1

      CONTAINS

!>  \brief subroutine (parallel and serial) for updating the SIESTA state
      SUBROUTINE update_state (lprint, fsq_total, ftol)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: fsq_total, ftol
      LOGICAL, INTENT(IN)   :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m1=1
      INTEGER     :: istat, js, n1, n2
      REAL(dp) :: ton, toff, r0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: work4
!-----------------------------------------------
      CALL second0(ton)
       
      CALL UpdateFields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,     &
                        djbsupsmnsh,djbsupumnch,djbsupvmnch,djpmnch) 
      IF (lasym)                                                        &
      CALL UpdateFields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh,     &
                        djbsupsmnch,djbsupumnsh,djbsupvmnsh,djpmnsh)

!     Gather new total fields onto all processors
      IF (PARSOLVER) CALL GatherFields
       
!     Reset perturbations
      CALL Clear_Field_Perts

      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

!     fsq_total == 0 when updating from add_resistivity
      IF (fsq_total.NE.zero .AND. fsq_total.LT.10*ftol)                 &
         CALL write_profiles(fsq_total)                                 ! SPH: write pmn, jbsupsmn profiles

      IF (.NOT.lprint) RETURN

      CALL second0(ton)
      CALL Update_Diagnostics(jbsupsmnsh, jbsupumnch, jbsupumnch,       &
                              jpmnch, bs0(1:6), bu0(1:6),               &
                              bsbu_ratio_s, pwr_spec_s, isym)
      IF (lasym)                                                        &
      CALL Update_Diagnostics(jbsupsmnch, jbsupumnsh, jbsupvmnsh,       &
                              jpmnsh, bs0(7:12), bu0(7:12),             &
                              bsbu_ratio_a, pwr_spec_a, iasym)
      lfirst = .FALSE.

      lcurr_init=.TRUE.                             
      CALL init_state(.TRUE.)                                            ! Compute co/contra variant currents for divJ, bdotj

      n1 = startglobrow;  n2 = endglobrow

      CALL divb(n1, n2)
      CALL divj(n1, n2)                                                   !(RS) Added, 09/10/07  = compute rms of divergence of J
      CALL bgradp(n1, n2)
      CALL tflux
      CALL bdotj_par      
      CALL beta
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         IF (.NOT.lverbose .AND. ISTAT.EQ.6) CYCLE
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),            &
                             divb_rms, toroidal_flux, wp/wb,            &
                             bgradp_rms, max_bgradp, min_bgradp,        &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110   FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,   &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                       &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,        &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ', 1pe11.3,         &
             ' B.GRAD-P (max): ', 1pe11.3,' B.GRAD-P (min): ',1pe11.3,  & 
           /,' (J*B)/|JxB| (rms): ', 1pe11.3,                           &
             ' (J_par)/|J_tot| (rms): ', 1pe11.3,                       &
             '   DIV-J (rms): ', 1pe11.3)
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)
       
      END SUBROUTINE update_state

!>  \brief Subroutine for updating the fields (B,p) of the SIESTA state
!!   based on dJBSUPXMNH perturbations computed in update_bfield and update_pressure
      SUBROUTINE UpdateFields(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,   &
                              djbsupsmnh,djbsupumnh,djbsupvmnh,djpmnh)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   ::                      &
          jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,                       &
         djbsupsmnh,djbsupumnh,djbsupvmnh,djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n1, n2
!-----------------------------------------------
      n1 = startglobrow !LBOUND(djbsupsmnh,3)
      n2 = endglobrow   !UBOUND(djbsupsmnh,3)

      jbsupsmnh(:,:,n1:n2) = jbsupsmnh(:,:,n1:n2) + djbsupsmnh(:,:,n1:n2)
      jbsupumnh(:,:,n1:n2) = jbsupumnh(:,:,n1:n2) + djbsupumnh(:,:,n1:n2)
      jbsupvmnh(:,:,n1:n2) = jbsupvmnh(:,:,n1:n2) + djbsupvmnh(:,:,n1:n2)
      jpmnh(:,:,n1:n2)     = jpmnh(:,:,n1:n2)     + djpmnh(:,:,n1:n2)

      END SUBROUTINE UpdateFields

  
!>  \brief Subroutine for writing SIESTA diagnostic (screen) output per iteration
      SUBROUTINE Update_Diagnostics(jbsupsmnh, jbsupumnh, jbsupvmnh,    &
                                    jpmnh, bs, bu, bsbu_ratio,          &
                                    pwr_spec, ipar)
      USE island_params, ns=>ns_i, hs=>hs_i
      USE descriptor_mod, ONLY: iam, SIESTA_COMM
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN) ::           &
                jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns,4), INTENT(INOUT) ::  pwr_spec
      REAL(dp), INTENT(OUT) :: bs(6), bu(6), bsbu_ratio
      INTEGER  :: ipar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*(4)      :: def(2) = (/'SYM ', 'ASYM'/)
      INTEGER, PARAMETER :: m1=1
      REAL(dp) :: r0, diag_b, diag_p, d1(12), d2(12)
      INTEGER  :: js, itype, n1, n2
!-----------------------------------------------
      n1 = startglobrow;  n2 = endglobrow

      IF (ipar .EQ. isym) THEN
         itype=1; r0 = hs/2
      ELSE
         itype=2; r0 =-hs/2
      END IF

      bs(1) = SQRT(SUM((jbsupsmnh(m1,:,2) - r0*jbsupumnh(m1,:,2))**2))/r0
      bu(1) = SQRT(SUM((jbsupsmnh(m1,:,2) + r0*jbsupumnh(m1,:,2))**2))/r0
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
      d1(1) = bs(1); d1(2) = bu(1)
      CALL MPI_BCAST(d1,2,MPI_REAL8,0,SIESTA_COMM,MPI_ERR)
      bs(1) = d1(1); bu(1) = d1(2)
      END IF
#endif
      IF (ABS(bu(1)) .GT. 1.E-10_dp) THEN
         bsbu_ratio = bs(1)/bu(1)
         d1 = 0
         DO js=n1, MIN(6, n2)
            d1(js)   = SQRT(SUM(jbsupsmnh(m1,:,js)**2))/ABS(r0)
            d1(js+6) = SQRT(SUM(jbsupumnh(m1,:,js)**2))
         END DO
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
         d2 = d1
         CALL MPI_ALLREDUCE(d2, d1, 12, MPI_REAL8, MPI_SUM,              &
                            SIESTA_COMM, MPI_ERR)
         END IF
#endif
         bs(1:6) = d1(1:6);  bu(1:6) = d1(7:12)
      ELSE 
         bs(1:6)=0; bu(1:6)=0; bsbu_ratio_s = 1
      END IF

#if defined(DUMP_STATE)
      IF (lfirst) THEN
         pwr_spec(:,:,n1:n2,1) = jbsupsmnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,2) = jbsupumnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,3) = jbsupvmnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,4) = jpmnh(:,:,n1:n2)
      ELSE
         d1(1) = SUM((jbsupsmnh(:,:,n1:n2)-pwr_spec(:,:,n1:n2,1))**2     &
                +    (jbsupumnh(:,:,n1:n2)-pwr_spec(:,:,n1:n2,2))**2     &
                +    (jbsupvmnh(:,:,n1:n2)-pwr_spec(:,:,n1:n2,3))**2)              
         d1(2) = SUM((jpmnh(:,:,n1:n2)-pwr_spec(:,:,n1:n2,4))**2)
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
         d2 = d1
         CALL MPI_ALLREDUCE(d2, d1, 2, MPI_REAL8, MPI_SUM,               &
                            SIESTA_COMM, MPI_ERR)
		 END IF
#endif		 
         diag_b = d1(1);  diag_p = d1(2)
         IF (iam .EQ. 0) THEN
         DO js = 6, 33, 27
         IF (.NOT.lverbose .AND. js.EQ.6) CYCLE
         WRITE(js, 100) 'POWER SPECTRA(', def(itype),                    &
                        ') -- dB: ', diag_b, ' dP: ', diag_p
         END DO 
         END IF
      END IF
 100  FORMAT(1x, 2a,2(a,1p,e10.3))
#endif
      END SUBROUTINE Update_Diagnostics


!>  \brief Subroutine for resetting the perturbations of the SIESTA state
      SUBROUTINE Clear_Field_Perts
!-----------------------------------------------
      djbsupsmnsh = 0
      djbsupumnch = 0
      djbsupvmnch = 0
      djpmnch     = 0

      IF (lasym) THEN

      djbsupsmnch = 0
      djbsupumnsh = 0
      djbsupvmnsh = 0
      djpmnsh     = 0

      END IF

      END SUBROUTINE Clear_Field_Perts

      END MODULE siesta_state
