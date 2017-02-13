!>  \brief Module contained subroutines for updating (from t to t+delta_t) the magnetic field and pressure as part of the SIESTA project.
!!  Stores updated values of PMNCF, BSUBMNF, BSUPMNF, and BSUPIJF in Quantities Module
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Sep 15, 2006
      MODULE siesta_state
#define DUMP_STATE
      USE stel_kinds
      USE quantities
      USE diagnostics_mod
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio_s, bsbu_ratio_a
      USE timer_mod, ONLY: time_update_state
      USE hessian, ONLY: gather_array
#if defined(DUMP_STATE)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                         &
         deljbss, deljbuc, deljbvc, deljbsc, deljbus, deljbvs,           &
         delpc, delps
#endif
      CONTAINS

!>  \brief Parallel subroutine for updating the SIESTA state
       SUBROUTINE update_state_par (lprint, fsq_total, ftol)
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE mpi_inc
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), INTENT(IN)  :: fsq_total, ftol
       LOGICAL, INTENT(IN)      :: lprint
#if defined(SKS)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(dp) :: ton, toff, tests, testa
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpss, tmpuc, tmpvc, tmppc
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpsc, tmpus, tmpvs, tmpps

       INTEGER :: nsmin, nsmax, istat
       LOGICAL :: lfirst
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      ALLOCATE (tmpss(0:mpol,-ntor:ntor,ns), tmpuc(0:mpol,-ntor:ntor,ns),&
                tmpvc(0:mpol,-ntor:ntor,ns), tmppc(0:mpol,-ntor:ntor,ns),&
                stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error #1 in update_state_par'

      
      CALL gather_array(djbsupsmnsh(:,:,nsmin:nsmax), tmpss)
      CALL gather_array(djbsupumnch(:,:,nsmin:nsmax), tmpuc)
      CALL gather_array(djbsupvmnch(:,:,nsmin:nsmax), tmpvc)
      CALL gather_array(djpmnch(:,:,nsmin:nsmax), tmppc)

      IF (l_asym) THEN

      ALLOCATE (tmpsc(0:mpol,-ntor:ntor,ns), tmpus(0:mpol,-ntor:ntor,ns),&
                tmpvs(0:mpol,-ntor:ntor,ns), tmpps(0:mpol,-ntor:ntor,ns),&
                stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error #2 in update_state_par'

      CALL gather_array(djbsupsmnch(:,:,nsmin:nsmax), tmpsc)
      CALL gather_array(djbsupumnsh(:,:,nsmin:nsmax), tmpus)
      CALL gather_array(djbsupvmnsh(:,:,nsmin:nsmax), tmpvs)
      CALL gather_array(djpmnsh(:,:,nsmin:nsmax), tmpps)

      END IF

      CALL Init_Allocate_Arrays (.FALSE.)

      djbsupsmnsh = tmpss
      djbsupumnch = tmpuc
      djbsupvmnch = tmpvc
      djpmnch     = tmppc
      DEALLOCATE (tmpss, tmpuc, tmpvc, tmppc)

      IF (l_asym) THEN
      djbsupsmnch = tmpsc
      djbsupumnsh = tmpus
      djbsupvmnsh = tmpvs
      djpmnsh     = tmpps

      DEALLOCATE (tmpsc, tmpus, tmpvs, tmpps)
      END IF

      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

      CALL update_state(.FALSE., fsq_total, ftol)

#if defined(DUMP_STATE)
      lfirst = .NOT.ALLOCATED(deljbss)
      IF (lfirst) THEN
         ALLOCATE(deljbss(0:mpol,-ntor:ntor,ns),                         &
                  deljbuc(0:mpol,-ntor:ntor,ns),                         &
                  deljbvc(0:mpol,-ntor:ntor,ns),                         &
                  delpc(0:mpol,-ntor:ntor,ns))
         deljbss = jbsupsmnsh; deljbuc = jbsupumnch
         deljbvc = jbsupvmnch; delpc = jpmnch
      END IF

      IF (.NOT.lprint) THEN
         CALL Init_Allocate_Arrays (.TRUE.)                               !Clears parallel perts
         RETURN
      END IF

      IF (l_asym) THEN
      IF (lfirst) THEN
         ALLOCATE(deljbsc(0:mpol,-ntor:ntor,ns),                         &
                  deljbus(0:mpol,-ntor:ntor,ns),                         &
                  deljbvs(0:mpol,-ntor:ntor,ns),                         &
                  delps(0:mpol,-ntor:ntor,ns))
         deljbsc = jbsupsmnch; deljbus = jbsupumnsh
         deljbvs = jbsupvmnsh; delps = jpmnsh
      END IF
      END IF
      IF (iam.EQ.0 .AND. .NOT.lfirst) THEN

         tests = SUM((jbsupsmnsh-deljbss)**2 + (jbsupumnch-deljbuc)**2 + &
                     (jbsupvmnch-deljbvc)**2)
         IF (l_asym) THEN     
         testa = SUM((jbsupsmnch-deljbsc)**2 + (jbsupumnsh-deljbus)**2 + &
                     (jbsupvmnsh-deljbvs)**2)
         ELSE
         testa = 0
         ENDIF
         PRINT 100,'dB-POWER(SYM): ', tests, ' dB-POWER(ASYM): ', testa

         tests = SUM((jpmnch-delpc)**2)
         IF (l_asym) THEN  
         testa = SUM((jpmnsh-delps)**2)
         ELSE
         testa = 0
         END IF
         PRINT 100,'dP-POWER(SYM): ', tests, ' dP-POWER(ASYM): ', testa
      END IF
 100  FORMAT(2(1x,a,1p,e10.3))
#endif

      CALL second0(ton)

      CALL init_state_par (.TRUE.)                          
!   Compute co-variant and contravariant components of current (times jacobian) - need for divJ, BdotJ
      lcurr_init=.TRUE.
      CALL divb_par
      CALL divj_par
      CALL bgradp_par
      CALL tflux
      CALL bdotj_par      
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),             &
                             divb_rms, toroidal_flux, wp/wb,             &
                             bgradp_rms, max_bgradp, min_bgradp,         &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110   FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,    &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                        &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,         &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ', 1pe11.3,          &
             ' B.GRAD-P (max): ', 1pe11.3,' B.GRAD-P (min): ',1pe11.3,   & 
           /,' (J*B)/|JxB| (rms): ', 1pe11.3,                            &
             ' (J_par)/|J_tot| (rms): ', 1pe11.3,                        &
             '   DIV-J (rms): ', 1pe11.3)

      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)
#endif       
      END SUBROUTINE update_state_par

!>  \brief Serial subroutine for updating the SIESTA state
      SUBROUTINE update_state (lprint, fsq_total, ftol)
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: fsq_total, ftol
       LOGICAL, INTENT(in)      :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m1=1
       INTEGER     :: istat, js, m, n
       REAL(dp) :: ton, toff, r0, xFSQ
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: work4
!-----------------------------------------------
       CALL second0(ton)

!      Update BSUP?MN?H*jacobh, pressure-h*jacobh harmonics based on vsup?mn?f perturbations
!  
       jbsupsmnsh = jbsupsmnsh + djbsupsmnsh                              ! update B^s harmonics  (half mesh)     
       jbsupumnch = jbsupumnch + djbsupumnch                              ! update B^u harmonics  (half mesh)  
       jbsupvmnch = jbsupvmnch + djbsupvmnch                              ! update B^v harmonics  (half mesh) 
       jpmnch     = jpmnch     + djpmnch                                  ! update pressure harmonics  (half mesh) 

       IF (l_asym) THEN

       jbsupsmnch = jbsupsmnch + djbsupsmnch                              ! update B^s harmonics  (half mesh)
       jbsupumnsh = jbsupumnsh + djbsupumnsh                              ! update B^u harmonics  (half mesh)
       jbsupvmnsh = jbsupvmnsh + djbsupvmnsh                              ! update B^v harmonics  (half mesh)
       jpmnsh     = jpmnsh     + djpmnsh                                  ! update pressure harmonics  (half mesh) 

       END IF

!      Reset perturbations
       CALL Clear_Field_Perts

       CALL second0(toff)
       time_update_state = time_update_state + (toff-ton)

!      fsq_total == 0 when updating from add_resistivity
       IF (fsq_total.NE.zero .AND. fsq_total.LT.10*ftol) THEN
          CALL write_profiles(fsq_total)                                 ! SPH: write pmn, jbsupsmn profiles
       END IF

      r0 = hs_i/2

      bs0(1) = SQRT(SUM((jbsupsmnsh(m1,:,2) - r0*jbsupumnch(m1,:,2))**2))/r0
      bu0(1) = SQRT(SUM((jbsupsmnsh(m1,:,2) + r0*jbsupumnch(m1,:,2))**2))/r0

      IF (l_asym) THEN
      bs0(7) = SQRT(SUM((jbsupsmnch(m1,:,2) + r0*jbsupumnsh(m1,:,2))**2))/r0
      bu0(7) = SQRT(SUM((jbsupsmnch(m1,:,2) - r0*jbsupumnsh(m1,:,2))**2))/r0
      END IF

      IF (bu0(1) .GT. 1.E-10_dp) THEN
         bsbu_ratio_s = bs0(1)/bu0(1)
         DO js = 2, MIN(ns,6)
            bs0(js) = SQRT(SUM(jbsupsmnsh(m1,:,js)**2))/r0
            bu0(js) = SQRT(SUM(jbsupumnch(m1,:,js)**2))
         END DO
      ELSE 
         bs0(1:6)=0; bu0(1:6)=0; bsbu_ratio_s = 1
      END IF

      IF (l_asym .AND. bu0(7) .GT. 1.E-10_dp) THEN
         bsbu_ratio_a = bs0(7)/bu0(7)
         DO js = 2, MIN(ns,6)
            bs0(js+6) = SQRT(SUM(jbsupsmnch(m1,:,js)**2))/r0
            bu0(js+6) = SQRT(SUM(jbsupumnsh(m1,:,js)**2))
         END DO
      ELSE 
         bs0(7:12)=0; bu0(7:12)=0; bsbu_ratio_a = 1
      END IF

      IF (.NOT.lprint) RETURN

      CALL second0(ton)

      lcurr_init=.TRUE.                             
      CALL init_state(.TRUE.)                                            ! Compute co/contra variant currents for divJ, bdotj

      CALL divb
      CALL divj                                                          ! (RS) Added, 09/10/07  = compute rms of divergence of J
      CALL bgradp
      CALL tflux
      CALL bdotj(bdotjijh)                                               ! (RS) Added, 08/30/07  = compute rms of parallel current
      DEALLOCATE(bdotjijh)
      CALL beta
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),             &
                             divb_rms, toroidal_flux, wp/wb,             &
                             bgradp_rms, max_bgradp, min_bgradp,         &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110   FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,    &
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


!>  \brief Subroutine for resetting the perturbations of the SIESTA state
      SUBROUTINE Clear_Field_Perts
      IMPLICIT NONE
!-----------------------------------------------
      djbsupsmnsh = 0
      djbsupumnch = 0
      djbsupvmnch = 0
      djpmnch     = 0

      IF (l_asym) THEN

      djbsupsmnch = 0
      djbsupumnsh = 0
      djbsupvmnsh = 0
      djpmnsh     = 0

      END IF

      END SUBROUTINE Clear_Field_Perts

      END MODULE siesta_state
