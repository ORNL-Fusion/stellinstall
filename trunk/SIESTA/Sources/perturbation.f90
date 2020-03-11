      MODULE PERTURBATION
      USE v3_utilities, ONLY: assert, assert_eq
      USE quantities
      USE descriptor_mod, ONLY: iam
      USE siesta_state, ONLY: update_state, update_state
      USE shared_data, ONLY: ngmres_type, iortho, lcolscale, lasym,            &
                             hesspass_test, buv_res, unit_out
      USE nscalingtools, ONLY: startglobrow, endglobrow

      IMPLICIT NONE
      
      PRIVATE

      INTEGER                                 :: nsmin
      INTEGER                                 :: nsmax
      INTEGER                                 :: irad_scale = 11
      INTEGER, PUBLIC                         :: niter_max                                              !width of perturbation (relative to 1)
      REAL(dp), ALLOCATABLE                   :: jdotb(:,:,:) !parallel current
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bmnc_res
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bmns_res     !resonant components of JDOTB

      PUBLIC :: init_data, add_perturb

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Load the namelist
!>
!>  READS siesta.jcf (JOB CONTROL FILE) for data (generalized to
!>  siesta_<name>.jcf, where <name> is the command line argument (ignored if not
!>  present!)
!>
!>  @returns The preconditioner control value.
!-------------------------------------------------------------------------------
      FUNCTION init_data()
      USE date_and_computer
      USE siesta_namelist

      IMPLICIT NONE

!  Declare Arguments
      INTEGER            :: init_data

!  Local variables
      INTEGER            :: istat
      INTEGER            :: index1
      INTEGER            :: imon
      INTEGER            :: numargs
      CHARACTER(LEN=10)  :: date0
      CHARACTER(LEN=10)  :: time0
      CHARACTER(LEN=10)  :: zone0
      CHARACTER(LEN=256) :: temp
      CHARACTER(LEN=100) :: command_arg
      CHARACTER(LEN=256) :: filename

!  Start of executable code

!  Read Command line arguments if present.
      CALL getcarg(1, command_arg, numargs)

      IF (numargs .GT. 0) THEN
         IF (LEN_TRIM(temp) .NE. 0) THEN
            IF (command_arg(1:7) .NE. 'siesta_' .AND.                          &
                command_arg(1:7) .NE. 'SIESTA_') THEN
               temp = 'siesta_' // TRIM(command_arg)
            END IF
            istat = LEN_TRIM(temp)
            IF (temp(istat - 3:istat) .NE. '.jcf' .AND.                        &
                temp(istat - 3:istat) .NE. '.JCF') THEN
               temp = TRIM(temp) // '.jcf'
            END IF
         END IF
      ELSE
         temp = 'siesta.jcf'
      END IF

      CALL siesta_namelist_read(temp)
      init_data = nprecon

      irad_scale = MIN(irad_scale, nsin)
      ngmres_type = MIN(2, MAX(1, ngmres_type))
      iortho = MIN(3, MAX(0, iortho))
      niter_max      = niter

!  Separate out the name form the wout file name. Wout files sould have the form
!  of wout_name.nc. Remove the leading wout_ and trailing .nc.
      IF ((wout_file(1:5) .ne. 'WOUT_' .and.                                   &
           wout_file(1:5) .ne. 'wout_') .or.  &
          wout_file(LEN_TRIM(wout_file) - 2:LEN_TRIM(wout_file)) .ne.          &
          '.nc') THEN
         CALL assert(.false., 'ERROR: ' // TRIM(wout_file) //                  &
     &                        ' is not a valid wout file name.')
      ELSE
         temp = wout_file(6:LEN_TRIM(wout_file) - 3)
      END IF

      istat = 0
      IF (iam .EQ. 0) THEN
	     WRITE (filename,1000) TRIM(temp), nsin, mpolin, ntorin
	     OPEN (unit=unit_out, file=filename, iostat=istat)
      END IF
      CALL assert_eq(istat, 0, 'ERROR WRITING SIESTA OUTPUT FILE')
      IF (iam .NE. 0) RETURN

      IF (lverbose) THEN
         WRITE (*,1001) '4.0', 100917
      END IF
      WRITE (unit_out,1001) '4.0', 100917

      WRITE (unit_out, 1002) TRIM(temp)
      IF (lverbose) THEN
         WRITE (*,1003) TRIM(temp)
      END IF

!  Format date and time.
      CALL DATE_AND_TIME(date0,time0,zone0)
      READ (date0(5:6),2000) imon
      WRITE (unit_out,1004) months(imon), date0(7:8), date0(1:4),              &
                            time0(1:2), time0(3:4), time0(5:6)

      DO istat = 1, SIZE(HelPert)
         IF (HelPert(istat).NE.zero .OR. HelPertA(istat).NE.zero) THEN
            WRITE (unit_out, 1005) istat, mres(istat), HelPert(istat),         &
                                   HelPertA(istat)
         END IF
      END DO

      WRITE (unit_out, 1006) ngmres_type, iortho, lColScale
      CALL FLUSH(unit_out)

1000  FORMAT('output_',a,'_',i4.4,2('X',i3.3),'.txt')
1001  FORMAT(62('-'),/,                                                        &
             1x,'SIESTA MHD EQUILIBRIUM CODE v',a,' (',i6,')', /,1x,           &
             'Scalable Island Equilibrium Solver for Toroidal Applications',/, &
             62('-'),/)
1002  FORMAT('CASE: ',a)
1003  FORMAT('OUTPUT FILE (SCREEN DUMP): output_',a,'.txt')
1004  FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2,/)
1005  FORMAT(i2,' mres: ',i4,' HelPert: ',1pe9.2,' HelPertA: ',1pe9.2)
1006  FORMAT(/,' ngmres_type: ', i4,' iOrtho: ', i4, ' lColScale: ', l2)

2000  FORMAT(i2)

      END FUNCTION

      SUBROUTINE add_perturb (xc, getwmhd)
      USE siesta_error
      USE fourier, ONLY: f_cos, f_sin, f_none
      USE siesta_namelist, ONLY: HelPert, HelPertA, eta_factor, lresistive, mres
      USE island_params, ONLY: fourier_context

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp)            :: xc(:)
      REAL(dp), EXTERNAL  :: getwmhd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER,PARAMETER  :: ipert_sign=2
      INTEGER            :: js, l, iprint, mres0, nres0, isign1, icount, &
                            irscale, imin(2), jstart, istat, jsave
      REAL(dp)           :: w0, w1, wmhd, eta_factor_save, p_width_min
      REAL(dp)           :: normal, HelPert0, HelPert0A, HP, rad,       &
                            chip0, phip0
      LOGICAL            :: lresist_save
      REAL(dp)           :: p_width

!-----------------------------------------------
!
!     Add helical RESISTIVE flux perturbation on full mesh
!     Approximate nonideal E_sub ~ f(mu+nv) B_sub in update_bfield (E ~ B)
!
      IF (ntor.EQ.0 .OR. (ALL(HelPert.EQ.zero)                          &
                    .AND. ALL(HelPertA.EQ.zero))) RETURN

!     INITIAL STATE: allocate buv_res so bsubXijf will be computed in init_state
      nsmin = MAX(1,startglobrow - 1)
      nsmax = MIN(ns,endglobrow + 1)

      ALLOCATE (buv_res(ntheta,nzeta,nsmin:nsmax),                      &
                jdotb(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,'ALLOCATION ERROR IN add_perturbation')
 
      lresist_save = lresistive
      lresistive=.FALSE.
      eta_factor_save = eta_factor
      eta_factor=1
      xc = 0
      w0 = getwmhd(xc)                                                   !Get latest current and bsubXijf
      
      IF (iam .EQ. 0) THEN
         DO iprint = 6, unit_out, unit_out-6
            IF (.NOT.lverbose .AND. iprint.EQ.6) CYCLE
            WRITE (iprint, '(/,a,/,2a)')                                   &
            ' Adding helical magnetic field perturbations',                &
            ' 10^6 X Del-W    mres    nres     HelPert',                   &
            '     rad  |m*chip+n*phip|     iota   radial width'
         END DO
      ENDIF

#if defined(JDOTB_PERT)
!     Compute jdotb to extract resonant components of sqrt(g)*J dot B
      jdotb = ksupsijf0(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)      &
            + ksupuijf0(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)      &
            + ksupvijf0(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

      ALLOCATE (bmnc_res(0:mpol,-ntor:ntor, nsmin:nsmax),               &
                bmns_res(0:mpol,-ntor:ntor, nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,'ISTAT != 0 IN ADD_PERT')
      
      bmns_res = 0
      CALL fourier_context%tomnsp(jdotb, bmnc_res, f_cos)
      IF (lasym) THEN
         CALL fourier_context%tomnsp(jdotb, bmns_res, f_sin)
      END IF
#endif      

!     Compute responses to resonant resistitive perturbations
      lresistive = .TRUE.
      normal = hs_i**2

      RES_TEST: DO icount = 1, SIZE(HelPert)
         HelPert0 = ABS(HelPert(icount)/normal)
         HelPert0A= ABS(HelPertA(icount)/normal)
         mres0 = ABS(mres(icount))
         IF ((HelPert0  .EQ. zero .AND.                                       &
              HelPert0A .EQ. zero) .OR. mres0.GT.mpol) CYCLE

!Scan in radius to determine primary resonance

         NRES_TEST: DO nres0 = -ntor, ntor
            jstart = 3   !Avoid 0/0 chip/phip at origin
            JSURF_TEST: DO WHILE (jstart .LT. ns-1)
         
               IF (.NOT.FindResonance(mres0, nres0, rad, jstart)) EXIT
               w1 = 2*w0
               imin(1) = 1
               imin(2) = 1

               PERT_SIGN: DO isign1 = 1, 2

!     NOTE: in update_bfield, where perturbed b-field is computed, this is multiplied
!           by eta_prof = rho*(1-rho) so it vanishes at both ends. Multiply again
!           here by that factor to assure radial derivatives entering in dB also vanish
!           at endpoints s=0,1
 
                  SCALE_TEST: DO irscale = 1, irad_scale
                     CALL GetResPert(irscale, isign1, mres0, nres0, rad,       &
                                     HelPert0, HelPert0A, chip0, phip0, p_width)
                     xc = 0
                     wmhd = getwmhd(xc)
                     IF (wmhd .LT. w1) THEN
                        imin(1) = irscale
                        imin(2) = isign1
                        w1 = wmhd
                        p_width_min = p_width
                     ENDIF
                  END DO SCALE_TEST
               END DO PERT_SIGN

!     Make sure energy decreases
               wmhd = w1
!     recompute perturbation buv_res here
               irscale = imin(1)
               isign1  = imin(2)
         
               IF (iam .EQ. 0) THEN
                  HP = HelPert0
                  IF (isign1 .EQ. 2) HP = -HP

                  DO iprint = 6, unit_out, unit_out-6
                  IF (.NOT.lverbose .AND. iprint.EQ.6) CYCLE
                  WRITE (iprint, 100) 10**6*(wmhd-w0), mres0, nres0,          &
                          HP*normal, rad, ABS(mres0*chip0+nres0*nfp*phip0),   &
                          chip0/phip0, p_width_min
                  END DO
 100  FORMAT(1p,e12.3,2(i8),3x,1pe10.2,0p,f8.2,f11.2,4x,2(f11.2))            
                  CALL FLUSH(unit_out)
               END IF
         
               CALL GetResPert(irscale, isign1, mres0, nres0, rad,            &
                               HelPert0, HelPert0A, chip0, phip0, p_width)

               wmhd = getwmhd(xc)
               IF (ABS(wmhd-w1) .GT. 1.E-12_dp) THEN
                  CALL siesta_error_set_error(siesta_error_general,           &
                                              'Error1 in Perturb')
               END IF
               xc = 0
               CALL update_state(.FALSE., zero, zero)
!Compute dW relative to perturbed state
               w0 = wmhd
   
            END DO JSURF_TEST                                                  ! Look for multiple resonances
         END DO NRES_TEST
      END DO RES_TEST

      IF (iam .EQ. 0) THEN
         DO iprint = 6, unit_out, unit_out-6
            IF (.NOT.lverbose .AND. iprint.EQ.6) CYCLE
            WRITE (iprint, *)
         END DO
      END IF

      DEALLOCATE (buv_res, jdotb, bmnc_res, bmns_res, stat=istat)
      lresistive = lresist_save
      eta_factor = eta_factor_save

      END SUBROUTINE add_perturb


      LOGICAL FUNCTION FindResonance(mres0, nres0, resrad, jstart)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)    :: mres0
      INTEGER, INTENT(in)    :: nres0
      INTEGER, INTENT(inout) :: jstart
      REAL(dp), INTENT(out)  :: resrad
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                :: js, jmin
      REAL(dp)               :: del1, del2, rho1, rho2, delmin
!-----------------------------------------------
!
!     Starting at jstart, find place where mres0*chip+n*nfp*phip is minimum (changes sign)
!
      del2 = 0
      resrad = 0
      delmin = -1
      IF (mres0*nres0 .EQ. 0) THEN
          FindResonance = .FALSE.
          RETURN
      END IF
!     Check for zero crossing
      jmin = jstart
      DO js = jstart, ns - 1
         del1 = mres0*chipf_i(js) + nres0*nfp*phipf_i(js)
         del1 = del1/MAX(ABS(mres0*chipf_i(js)),ABS(nres0*nfp*phipf_i(js)),1.E-20_dp)
         IF (delmin .EQ. -1 .OR. ABS(del1) .LT. delmin) THEN
            IF (delmin .EQ. -1) del2 = del1
            delmin = ABS(del1)
            jmin = js
         END IF
         
         IF (del1*del2 < zero) THEN
            jstart = js+1
            rho2 = hs_i*(js-2)
            rho1 = hs_i*(js-1)
            resrad = (rho2*ABS(del1) + rho1*ABS(del2))/(ABS(del1) + ABS(del2))
            delmin=0
            EXIT
         END IF
         del2 = del1
      END DO

!IF NO ZERO-CROSSING, RESONANCE MIGHT BE VERY CLOSE
      IF (delmin .LT. 1.E-3_dp) THEN
         jstart = jmin+1
         resrad = hs_i*(jmin-1)
         FindResonance=.TRUE.
      ELSE
         jstart=ns-1
         resrad = 0
         FindResonance=.FALSE.
      END IF

      END FUNCTION FindResonance


      SUBROUTINE GetResPert(irscale, isign1, mres0, nres0, rad,                &
                            HelPert0, HelPert0A, chip0, phip0, p_width)
      USE fourier, ONLY: f_cos, f_sin, f_sum, f_none
      USE siesta_namelist, ONLY: mpolin, ntorin
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)    :: irscale, isign1
      INTEGER, INTENT(in)    :: mres0, nres0
      REAL(dp), INTENT(in)   :: rad
      REAL(dp), INTENT(out)  :: chip0, phip0
      REAL(dp), INTENT(in)   :: HelPert0, HelPert0A
      REAL(dp), INTENT(out)  :: p_width
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: js, istat
      REAL(dp)    :: rho, rhores, pert_prof(ns), HelP_local,                   &
                     HelPA_local, locrad
      REAL(dp), ALLOCATABLE   :: bmn_res(:,:)
      REAL(dp), ALLOCATABLE   :: resc(:,:,:), ress(:,:,:)
!-----------------------------------------------
      IF (mres0.GT.mpolin .OR. ABS(nres0).GT.ntorin) RETURN

      HelP_local = HelPert0
      HelPA_local = HelPert0A
      IF (isign1 .EQ. 2) THEN
         HelP_local = -HelPert0
         HelPA_local = -HelPert0A
      END IF

!	Compute radial form factor (window-function)
      js = INT(rad/hs_i) + 1
      IF (js.LT.1 .OR. js.GT.ns) RETURN
      chip0 = chipf_i(js)
      phip0 = phipf_i(js)
      locrad = REAL(irscale - 1,dp)/(irad_scale - 1)
      p_width = 0.1_dp + 0.9_dp*locrad                !decay length, ~delta-fcn for irscale=1 to broad (full radius)
      locrad = one/p_width
      DO js = 1,ns
         rho = hs_i*(js-1)
         rhores = rho-rad
         rhores = locrad*rhores
         pert_prof(js) = one/(one + rhores*rhores)      !tearing parity; mult by alpha for equal area
      END DO

!THIS DOESN'T SEEM TO WORK AS WELL AS THE SIMPLE FORM FACTOR
!SO KEEP JDOTB_PERT UNDEFINED UNTIL WE CAN IMPROVE IT
#if defined(JDOTB_PERT)
      ALLOCATE(resc(0:mpol,-ntor:ntor, nsmin:nsmax),                    &
               ress(0:mpol,-ntor:ntor, nsmin:nsmax), stat=istat)

      rho = MAX(MAXVAL(ABS(bmnc_res(mres0,nres0,:))),                   &
                MAXVAL(ABS(bmns_res(mres0,nres0,:))))

      resc = 0
      ress = 0
      resc(mres0,nres0,:) = HelP_local*bmnc_res(mres0,nres0,:)
      
      CALL fourier_context%toijsp(resc, buv_res, f_none, f_cos)
      
      IF (lasym) THEN
         ress(mres0,nres0,:) = HelPA_local*bmns_res(mres0,nres0,:)
         CALL fourier_context%toijsp(ress, buv_res, f_sum, f_sin)
      END IF
      
      IF (rho .NE. zero) THEN
         buv_res = buv_res/rho
      END IF

!	Apply radial form-factor
      DO js = nsmin, nsmax
         buv_res(:,:,js) = buv_res(:,:,js)*pert_prof(js)
      END DO

      DEALLOCATE (resc, ress, stat=istat)
#else      
      ALLOCATE (bmn_res(0:mpol,-ntor:ntor), stat=istat)
      CALL ASSERT(istat.eq.0,'ISTAT != 0 IN GETRESPERT')
      bmn_res = 0
         
      bmn_res(mres0,nres0) = HelP_local
      CALL fourier_context%toijsp(bmn_res, buv_res(:,:,nsmin), f_none, f_cos)
      
      IF (lasym) THEN
         bmn_res(mres0,nres0) = HelPA_local
         CALL fourier_context%toijsp(bmn_res, buv_res(:,:,nsmin), f_sum, f_sin)
      END IF
      DO js = nsmin + 1, nsmax
         buv_res(:,:,js) = buv_res(:,:,nsmin)*pert_prof(js)
      END DO

      buv_res(:,:,nsmin) = buv_res(:,:,nsmin)*pert_prof(nsmin)
      
      DEALLOCATE (bmn_res, stat=istat)
#endif

      END SUBROUTINE GetResPert
  
      SUBROUTINE Compute_Resonant_Epar(xc, epar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(out) :: epar(ns,ntheta,nzeta)
      REAL(dp), INTENT(in)  :: xc(ns,0:mpol,-ntor:ntor,3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE :: eparmn(:,:,:)
      REAL(dp)              :: rad
      INTEGER               :: mres, nres, istat, jstart
!-----------------------------------------------
!
!     xc(1,...,1) ARE THE RESONANT AMPLITUDES OF E||
!     eparmn      RESONANT RESPONSE TO xc
!     epar        IS THE RESULTING E*B/B^2 IN REAL SPACE
!
      ALLOCATE (eparmn(ns,0:mpol,-ntor:ntor), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in Compute_Resonant_Epar')
      eparmn = 0
      epar = 0

      DO mres=0,mpol
         DO nres=-ntor,ntor
            jstart = 2
            DO WHILE (jstart .lt. ns-1)
               IF (.NOT.FindResonance(mres, nres, rad, jstart)) EXIT
			   CALL ASSERT(.FALSE.,'AddResonantE not implemented')
!NOT IMPLEMENTED     CALL AddResonantE (eparmn(:,mres,nres), rad, jstart)
            END DO
         END DO
      END DO

      DEALLOCATE(eparmn)

      END SUBROUTINE Compute_Resonant_Epar


      END MODULE PERTURBATION
