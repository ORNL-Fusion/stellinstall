!>  \brief Module contained subroutines for updating (from t to t+delta_t) the magnetic field and pressure as part of the SIESTA project.
!!  Stores updated values of PMNCF, BSUBMNF, BSUPMNF, and BSUPIJF in Quantities Module
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Sep 15, 2006
      MODULE siesta_state
#define DUMP_STATE
      USE stel_kinds
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio_s, bsbu_ratio_a, lverbose
      USE timer_mod, ONLY: time_update_state
      USE siesta_init, ONLY: init_state_par, init_allocate_arrays
      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE :: isym=0, iasym=1
!      LOGICAL, PRIVATE :: lfirst=.TRUE.
      LOGICAL :: lfirst=.TRUE.

      CONTAINS

!>  \brief Parallel subroutine for updating the SIESTA state
      SUBROUTINE update_state_par (lprint, fsq_total, ftol)
      USE quantities
      USE shared_data, ONLY: delta_t !!!!!!!MRC
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
#if defined(SKS)
      CALL GATHER_FIELDS(djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch)
      IF (lasym)                                                         &
      CALL GATHER_FIELDS(djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh)
#endif
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

      CALL update_state(lprint, fsq_total, ftol)

      IF (.NOT.lprint) THEN
         CALL Init_Allocate_Arrays (.TRUE.)                               !Realloc perturbation arrays
         RETURN
      END IF

      END SUBROUTINE update_state_par

!>  \brief Serial subroutine for updating the SIESTA state
      SUBROUTINE update_state (lprint, fsq_total, ftol)
      USE quantities
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
      jbsupsmnh = jbsupsmnh + djbsupsmnh                                   ! update B^s harmonics  (half mesh)     
      jbsupumnh = jbsupumnh + djbsupumnh                                   ! update B^u harmonics  (half mesh)  
      jbsupvmnh = jbsupvmnh + djbsupvmnh                                   ! update B^v harmonics  (half mesh) 
      jpmnh     = jpmnh     + djpmnh                                       ! update pressure harmonics  (half mesh) 

!Reset perturbations
      djbsupsmnh = 0
      djbsupumnh = 0
      djbsupvmnh = 0
      djpmnh     = 0

      END SUBROUTINE Update_Fields


      SUBROUTINE Update_Diagnostics(jbsupsmnh, jbsupumnh, jbsupvmnh,     &
                                    jpmnh, bs, bu, bsbu_ratio, pwr_spec, &
                                    iparity)
      USE quantities, ONLY: mpol, ntor, iam
      USE island_params, ns=>ns_i, hs=>hs_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN) ::           &
                jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh
      REAL(dp), DIMENSION(0:mpol,-ntor:ntor,ns,4) ::  pwr_spec
      REAL(dp), INTENT(OUT) :: bs(6), bu(6), bsbu_ratio
      INTEGER  :: iparity
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*(4)      :: def(2) = (/'SYM ', 'ASYM'/)
      INTEGER, PARAMETER :: m1=1
      REAL(dp) :: r0, diag_b, diag_p
      INTEGER  :: js, itype
!-----------------------------------------------
      IF (iparity .EQ. isym) THEN
         itype=1; r0 = hs/2
      ELSE
         itype=2; r0 =-hs/2
      END IF

      bs(1) = SQRT(SUM((jbsupsmnh(m1,:,2) - r0*jbsupumnh(m1,:,2))**2))/r0
      bu(1) = SQRT(SUM((jbsupsmnh(m1,:,2) + r0*jbsupumnh(m1,:,2))**2))/r0

      IF (ABS(bu(1)) .GT. 1.E-10_dp) THEN
         bsbu_ratio = bs(1)/bu(1)
         DO js = 2, MIN(ns,6)
            bs(js) = SQRT(SUM(jbsupsmnh(m1,:,js)**2))/ABS(r0)
            bu(js) = SQRT(SUM(jbsupumnh(m1,:,js)**2))
         END DO
      ELSE 
         bs(1:6)=0; bu(1:6)=0; bsbu_ratio_s = 1
      END IF

#if defined(DUMP_STATE)
      IF (lfirst) THEN
         pwr_spec(:,:,:,1) = jbsupsmnh; pwr_spec(:,:,:,2) = jbsupumnh
         pwr_spec(:,:,:,3) = jbsupvmnh; pwr_spec(:,:,:,4) = jpmnh
      ELSE IF (iam .EQ. 0) THEN
         diag_b = SUM((jbsupsmnh-pwr_spec(:,:,:,1))**2                   &
                +     (jbsupumnh-pwr_spec(:,:,:,2))**2                   &
                +     (jbsupvmnh-pwr_spec(:,:,:,3))**2)              
         diag_p = SUM((jpmnh-pwr_spec(:,:,:,4))**2)
         DO js = 6, 33, 27
         IF (.not.lverbose .and. js .eq. 6) CYCLE
         WRITE(js, 100) 'POWER SPECTRA(', def(itype),                    &
                        ') -- dB: ', diag_b, ' dP: ', diag_p
         END DO
      END IF
 100  FORMAT(1x, 2a,2(a,1p,e10.3))
#endif


      END SUBROUTINE Update_Diagnostics


!>  \brief Subroutine for resetting the perturbations of the SIESTA state
      SUBROUTINE Clear_Field_Perts
      USE quantities
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


