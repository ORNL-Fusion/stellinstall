!>  \brief Module contained subroutines for updating (from t to t+delta_t) the magnetic field and pressure as part of the SIESTA project.
!!  Stores updated values of PMNCF, BSUBMNF, BSUPMNF, and BSUPIJF in Quantities Module
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Sep 15, 2006
      MODULE siesta_state
#define DUMP_STATE
      USE stel_kinds
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio_s, bsbu_ratio_a, lverbose
      USE quantities
      USE timer_mod, ONLY: time_update_state
      USE siesta_init, ONLY: init_state_par, init_allocate_arrays
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
!      LOGICAL, PRIVATE :: lfirst=.TRUE. ! MRC need to reset this value in run_siesta
      LOGICAL :: lfirst=.TRUE.

      CONTAINS

!>  \brief Parallel subroutine for updating the SIESTA state
      SUBROUTINE update_state_par (lprint, fsq_total, ftol)
!-----------------------------------------------
!   D u m m y   A r g u m e n t sf
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: fsq_total, ftol
      LOGICAL, INTENT(IN)   :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)

      CALL PadFields
      CALL update_state(lprint, fsq_total, ftol)
 
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

      IF (.NOT.lprint) THEN
         CALL Init_Allocate_Arrays (.TRUE.)                               !Realloc perturbation arrays
         RETURN
      END IF

      END SUBROUTINE update_state_par

!>  \brief Serial subroutine for updating the SIESTA state
      SUBROUTINE update_state (lprint, fsq_total, ftol)
      USE siesta_init, ONLY: init_state
      USE diagnostics_mod, ONLY: toroidal_flux, toroidal_flux0,          &
          divb_rms, bgradp_rms, max_bgradp, min_bgradp, bdotj_rms,       &
          bdotj2_rms, divj_rms, lcurr_init, divb, divj, bgradp, tflux,   &
          bdotj, beta, write_profiles, divb_par, divj_par, bgradp_par,   &
          bdotj_par
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(in)  :: fsq_total, ftol
      LOGICAL, INTENT(in)   :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat, js, m, n
      REAL(dp)    :: ton, toff, r0, testa, tests
!-----------------------------------------------
      CALL second0(ton)

!Update BSUP?MN?H*jacobh, pressure-h*jacobh harmonics based on vsup?mn?f perturbations
      CALL Update_Fields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,     &
                        djbsupsmnsh,djbsupumnch,djbsupvmnch,djpmnch)
      IF (lasym)                                                         &
      CALL Update_Fields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh,     &
                        djbsupsmnch,djbsupumnsh,djbsupvmnsh,djpmnsh)

      CALL Clear_Field_Perts

      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

!fsq_total == 0 when updating from add_resistivity
      IF (fsq_total.NE.zero .AND. fsq_total.LT.10*ftol)                  &
      CALL write_profiles(fsq_total)                                       ! write pmn, jbsupsmn profiles

      IF (.NOT.lprint) RETURN

      CALL Update_Diagnostics(jbsupsmnsh, jbsupumnch, jbsupumnch,        &
                              jpmnch, bs0(1:6), bu0(1:6), bsbu_ratio_s,  &
                              pwr_spec_s, isym)
      IF (lasym)                                                         &
      CALL Update_Diagnostics(jbsupsmnch, jbsupumnsh, jbsupvmnsh,        &
                              jpmnsh, bs0(7:12), bu0(7:12), bsbu_ratio_a,&
                              pwr_spec_a, iasym)
      lfirst = .FALSE.

      lcurr_init=.TRUE.                             
      CALL init_state_par (.TRUE.)                                        !Update to new state in real space

      CALL divb_par
      CALL divj_par                                                       !(RS) Added, 09/10/07  = compute rms of divergence of J
      CALL bgradp_par
      CALL tflux
      CALL bdotj_par                                                      !(RS) Added, 08/30/07  = compute rms of parallel current
!      CALL beta
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         IF (.not.lverbose .and. ISTAT .eq. 6) CYCLE
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),             &
                             divb_rms, toroidal_flux, wp/wb,             &
                             bgradp_rms, max_bgradp, min_bgradp,         &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110  FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,     &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                        &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,         &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ', 1pe11.3,          &
             ' B.GRAD-P (max): ', 1pe11.3,' B.GRAD-P (min): ',1pe11.3,   & 
           /,' (J*B)/|JxB| (rms): ', 1pe11.3,                            &
             ' (J_par)/|J_tot| (rms): ', 1pe11.3,                        &
             '   DIV-J (rms): ', 1pe11.3)
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)
       
      END SUBROUTINE update_state


!>  \brief Subroutine for updating the fields (B,p) of the SIESTA state
      SUBROUTINE Update_Fields(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,   &
                              djbsupsmnh,djbsupumnh,djbsupvmnh,djpmnh)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   ::                       &
          jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,                        &
         djbsupsmnh,djbsupumnh,djbsupvmnh,djpmnh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n1, n2
!-----------------------------------------------
      n1=MAX(1,startglobrow-1); n2=MIN(ns,endglobrow+2)
!      IF (n1 .NE. LBOUND(djbsupsmnh,3)) CALL ErrorAbort('Update_Fields djbsupX wrong LB')
!      IF (n2 .NE. UBOUND(djbsupsmnh,3)) CALL ErrorAbort('Update_Fields djbsupX wrong UB')
      IF (LBOUND(jbsupsmnh,3) .NE. 1)   CALL ErrorAbort('Update_Fields jbsupX wrong LB')
      IF (UBOUND(jbsupsmnh,3) .NE. ns)  CALL ErrorAbort('Update_Fields jbsupX wrong UB')

      jbsupsmnh(:,:,n1:n2) = jbsupsmnh(:,:,n1:n2) + djbsupsmnh(:,:,n1:n2)  ! update B^s harmonics  (half mesh)     
      jbsupumnh(:,:,n1:n2) = jbsupumnh(:,:,n1:n2) + djbsupumnh(:,:,n1:n2)  ! update B^u harmonics  (half mesh)  
      jbsupvmnh(:,:,n1:n2) = jbsupvmnh(:,:,n1:n2) + djbsupvmnh(:,:,n1:n2)  ! update B^v harmonics  (half mesh) 
      jpmnh(:,:,n1:n2)     = jpmnh(:,:,n1:n2)     + djpmnh(:,:,n1:n2)      ! update pressure harmonics  (half mesh) 

      END SUBROUTINE Update_Fields


      SUBROUTINE Update_Diagnostics(jbsupsmnh, jbsupumnh, jbsupvmnh,     &
                                    jpmnh, bs, bu, bsbu_ratio, pwr_spec, &
                                    iparity)
      USE island_params, ns=>ns_i, hs=>hs_i
      USE mpi_inc
      USE descriptor_mod, ONLY: SIESTA_COMM
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN) ::           &
                jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns,4), INTENT(INOUT) ::  pwr_spec
      REAL(dp), INTENT(OUT) :: bs(6), bu(6), bsbu_ratio
      INTEGER  :: iparity
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*(4)      :: def(2) = (/'SYM ', 'ASYM'/)
      INTEGER, PARAMETER :: m1=1
      REAL(dp) :: r0, diag_b, diag_p, d1(12), d2(12)
      INTEGER  :: js, itype, n1, n2
!-----------------------------------------------
      n1 = startglobrow;  n2 = endglobrow

      IF (iparity .EQ. isym) THEN
         itype=1; r0 = hs/2
      ELSE
         itype=2; r0 =-hs/2
      END IF

      bs(1) = SQRT(SUM((jbsupsmnh(m1,:,2) - r0*jbsupumnh(m1,:,2))**2))/r0
      bu(1) = SQRT(SUM((jbsupsmnh(m1,:,2) + r0*jbsupumnh(m1,:,2))**2))/r0
#if defined(MPI_OPT)
      d1(1) = bs(1); d1(2) = bu(1)
      CALL MPI_BCAST(d1,2,MPI_REAL8,0,MPI_COMM_WORLD,MPI_ERR)
      bs(1) = d1(1); bu(1) = d1(2)
#endif
      IF (ABS(bu(1)) .GT. 1.E-10_dp) THEN
         bsbu_ratio = bs(1)/bu(1)
         d1 = 0
         DO js=n1, MIN(6, n2)
            d1(js)   = SQRT(SUM(jbsupsmnh(m1,:,js)**2))/ABS(r0)
            d1(js+6) = SQRT(SUM(jbsupumnh(m1,:,js)**2))
         END DO
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(d1, d2, 12, MPI_REAL8, MPI_SUM,              &
                            SIESTA_COMM, MPI_ERR)
#else
         d2 = d1
#endif
         bs(1:6) = d2(1:6);  bu(1:6) = d2(7:12)
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
         CALL MPI_ALLREDUCE(d1, d2, 2, MPI_REAL8, MPI_SUM,               &
                            MPI_COMM_WORLD, MPI_ERR)
#else
         d2 = d1
#endif
         diag_b = d2(1);  diag_p = d2(2)
         IF (iam .EQ. 0) THEN
         DO js = 6, 33, 27
         IF (.not.lverbose .and. js .eq. 6) CYCLE
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


