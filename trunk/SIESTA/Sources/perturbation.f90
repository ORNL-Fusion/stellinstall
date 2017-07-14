      MODULE PERTURBATION
      USE quantities
      USE descriptor_mod, ONLY: iam
      USE siesta_state, ONLY: update_state, update_state_par
      USE shared_data, ONLY: ngmres_type, iortho, lcolscale,             &
                             lasym, hesspass_test, mupar_test,           &
                             buv_res, lrecon
!#if defined(SKS)
      USE nscalingtools, ONLY: PARFUNCTISL
!#endif
      
      INTEGER  :: irad_scale=11
      INTEGER  :: nsin=101, mpolin=12, ntorin=3
      INTEGER  :: mres1, nres1, mres2, nres2                                !Backwards compat
      INTEGER  :: mres(20)
      INTEGER  :: niter_max
      REAL(dp) :: HelPert1, HelPert2, rad1, rad2                            !Backwards compat
      REAL(dp) :: HelPert(20)
      REAL(dp) :: ftol = 1.E-20_dp
      REAL(dp) :: eta_factor = 1.E-2_dp
      REAL(dp) :: p_width                                                   !width of perturbation (relative to 1)
      LOGICAL  :: ladd_pert = .TRUE.
      LOGICAL  :: lresistive = .TRUE.                                       !RS: If true, add resistive electric field.
      LOGICAL  :: lrestart = .FALSE.
      CHARACTER(LEN=256) :: wout_file, restart_ext

!NETCDF PARAMETERS
      CHARACTER(LEN=*), PARAMETER :: vn_nsin='nrad', vn_mpolin='mpol',   &
         vn_ntorin='ntor', vn_wout='wout file', vn_lasym='lasym',        &
         vn_bsupss='JBsupss(m,n,r)',vn_bsupuc='JBsupuc(m,n,r)',          &
         vn_bsupvc='JBsupvc(m,n,r)', vn_bsupsc='JBsupsc(m,n,r)',         &
         vn_bsupus='JBsupus(m,n,r)', vn_bsupvs='JBsupvs(m,n,r)',         &
         vn_presc='presc(m,n,r)', vn_press='press(m,n,r)'
      CHARACTER(LEN=*), DIMENSION(3), PARAMETER ::                       &
                  r3dim = (/'m-mode','n-mode','radius'/)

!SCRATCH ARRAYS
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupss, bsupuc, bsupvc, &
                presc, bsupsc, bsupus, bsupvs, press

!OUTPUT FILE PARAMETERS
      INTEGER:: nrs=200, nzs=200, nphis=2                                   ! BField resolution (CYL):  R and Z
      INTEGER:: nvs=100, nus=150, nss=150                                   ! BField resolution (VMEC): V (and PHI), U, S      
      LOGICAL :: l_tracing = .FALSE.                                        ! Do not produce file with CYL field for Poincare plots.
      LOGICAL :: l_silo_output = .FALSE.
      LOGICAL :: l_silo3D = .FALSE.
      LOGICAL :: l_output_alliter = .FALSE.                                 ! Default: do not dump output at every iteration on file (RS, Jul'10)
      LOGICAL :: l_VMEC_Uniform

      CHARACTER(LEN=256) :: filename

      CONTAINS

      SUBROUTINE init_data (nprecon)
      USE date_and_computer
      USE Hessian, ONLY: levmarq_param, mupar, levmarq_param0, mupar0
      IMPLICIT NONE
      INTEGER, INTENT(out) :: nprecon
      INTEGER              :: istat, index1, imon, niter
      CHARACTER(LEN=10)    :: date0, time0, zone0
      CHARACTER(LEN=256)   :: temp, short_name
      NAMELIST /SIESTA_INFO/ ladd_pert, lresistive, lrestart, niter,     &
               nsin, mpolin, ntorin, ftol, wout_file, restart_ext,       &
               mres1, nres1, rad1, helpert1, mres2, nres2, rad2,         &
               helpert2, mres, HelPert, levmarq_param, mupar,            &
               eta_factor, nphis, nrs, nzs, nvs, nus, nss, l_tracing,    &
               ngmres_type, iortho, lcolscale,                           &
               l_silo_output, l_silo3D, nprecon, l_output_alliter,       & 
               l_VMEC_Uniform, lasym, HessPass_Test, mupar_test, lrecon

!-----------------------------------------------
!     
!     READS siesta.jcf (JOB CONTROL FILE) for data
!
      levmarq_param = 1.E-3_dp
      mupar         = 0
      niter         = 10                !Number iterations after diagonal prec
      HelPert1      = 0
      HelPert2      = 0
      mres          = 0
      HelPert       = 0
      lresistive    = .TRUE.
      lcolscale     = .TRUE.
      mupar_test    = 0
      nprecon       = 0                 !set > 0 to skip diagonal preconditioner

      OPEN (unit=30, file='siesta.jcf', delim='apostrophe', iostat=istat)
      IF (istat .NE. 0) STOP 'SIESTA.JCF file not found'

      READ (30, NML=SIESTA_INFO, iostat=istat) 
      IF (istat .ne. 0) THEN
         IF (iam .EQ. 0) THEN
            STOP 'Error reading siesta.jcf file'
         ELSE
            STOP
         END IF
      END IF
      CLOSE (UNIT=30)

      irad_scale = MIN(irad_scale, nsin)
      ngmres_type = MIN(2, MAX(1,ngmres_type))
      iortho = MIN(3, MAX(0,iortho))
      nss = nsin
      levmarq_param0 = levmarq_param
      mupar0         = mupar
      niter_max      = niter

      temp = wout_file
      index1 = INDEX(temp,'WOUT')
      IF (index1 .eq. 0) index1 = INDEX(temp,'wout')
      IF (index1 .eq. 0) STOP 'ERROR: WRONG WOUT FILE NAME'
      temp = wout_file(index1:)
!
!     STRIP FILE EXTENSION: NAME COULD INCLUDE '.' IN IT!
!
      index1 = SCAN(temp, '.',BACK=.TRUE.)
      IF (index1 == 0) THEN
         temp = TRIM(temp(6:))
      ELSE
         temp = TRIM(temp(6:index1-1))
      END IF

      IF (HelPert1.ne.zero .or. HelPert2.ne.zero) THEN
         HelPert(1) = HelPert1; mres(1) = mres1
         HelPert(2) = HelPert2; mres(2) = mres2
      END IF

      IF (iam .ne. 0) RETURN

!     WOUT FILE MUST BE OF FORM "wout_NAME.nc (or txt)" or "wout_NAME"
      WRITE (filename,'(3a,i4.4,a,2(i3.3,a))') "output_", TRIM(temp),    &
             '_',nsin,'X',mpolin,'X',ntorin,".txt"
      OPEN (unit=33, file=filename,iostat=istat)
      IF (istat .ne. 0) STOP 'ERROR WRITING SIESTA OUTPUT FILE'

      DO istat = 6, 33, 27
	     WRITE (istat, 10)                                           & 
               'SIESTA MHD EQUILIBRIUM CODE v3.0 (091216)',              &
               'Scalable Island Equilibrium ',                           &
               'Solver for Toroidal Applications' 
         IF (HelPert1.ne.zero .or. HelPert2.ne.zero) THEN
            WRITE (istat, *) 'Please update jcf file to new format'
            WRITE (istat, *) 'Replace mres1,2 with mres array, HelPert1,2 with HelPert array'
         END IF
      END DO

10    FORMAT (62('-'),/,1x,a,/,2(1x,a),/,62('-'),/)

      WRITE (33, *)'CASE: ', TRIM(temp)
      PRINT *,'OUTPUT FILE (SCREEN DUMP): ', "output_" // TRIM(temp) // ".txt"

      CALL DATE_AND_TIME(date0,time0,zone0)
      READ (date0(5:6),'(i2)')imon
      WRITE (33,20) months(imon),date0(7:8),date0(1:4),                  &
     &              time0(1:2),time0(3:4),time0(5:6)

20    FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2,/)

      DO istat = 1, SIZE(HelPert)
         IF (HelPert(istat) .ne. zero) THEN
            WRITE (33, 30) istat, mres(istat), HelPert(istat)
         END IF
      END DO

      WRITE (33, 32) ngmres_type, iortho, lColScale

30    FORMAT(i2,' mres: ',i4,' HelPert: ',1pe9.2)
32    FORMAT(/,' ngmres_type: ', i4,' iOrtho: ', i4, ' lColScale: ', l2)

      END SUBROUTINE init_data

      SUBROUTINE add_perturb (xc, getwmhd)
      USE nscalingtools, ONLY: startglobrow, endglobrow
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp)            :: xc(:)
      REAL(dp), EXTERNAL  :: getwmhd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER,PARAMETER  :: ipert_sign=2
      INTEGER            :: js, l, iprint, mres0, nres0, isign, icount,  &
                            irscale, imin(2), jstart, istat, jsave
      INTEGER            :: nsmin, nsmax
      REAL(dp)           :: w0, w1, wmhd, eta_factor_save, p_width_min
      REAL(dp)           :: normal, HelPert0, HP, rad, chip0, phip0
      LOGICAL            :: lresist_save
!-----------------------------------------------
!
!     Add helical RESISTIVE flux perturbation on full mesh
!     Approximate nonideal E|| ~ f(mu+nv) B with B ~ Bsubv
!
      IF (ntor.EQ.0 .OR. ALL(HelPert.EQ.zero)) RETURN

!     INITIAL STATE
 
      lresist_save = lresistive
      lresistive=.FALSE.
      eta_factor_save = eta_factor
      eta_factor=1
      xc = 0
      IF (PARFUNCTISL) THEN
         nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)
      END IF
      w0 = getwmhd(xc, PARFUNCTISL)

      ALLOCATE (buv_res(ntheta,nzeta,ns), stat=istat)
      IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN add_perturbation'

      IF (iam .eq. 0) THEN
      DO iprint = 6, 33, 27
         WRITE (iprint, '(/,a,/,2a)')                                    &
         ' Adding helical magnetic field perturbations',                 &
         ' 10^6 X Del-W    mres    nres     HelPert',                    &
         '     rad  |m*chip+n*phip|     iota   radial width'
      END DO
      ENDIF

!     Resonant resistitive perturbations first
      lresistive=.TRUE.
      normal = hs_i**2

      RES_TEST: DO icount = 1, SIZE(HelPert)
         HelPert0 = ABS(HelPert(icount)/normal)
         mres0 = ABS(mres(icount))
         IF (HelPert0.EQ.zero .OR. mres0.GT.mpol) CYCLE

!Scan in radius to determine primary resonance

         NRES_TEST: DO nres0=-ntor,ntor
         jstart = 3   !Avoid 0/0 chip/phip at origin
         JSURF_TEST: DO WHILE (jstart .LT. ns-1)
         
         IF (.NOT.FindResonance(mres0, nres0, rad, jstart)) EXIT
         w1 = 2*w0         

      PERT_SIGN: DO isign = 1, 2

!     NOTE: in update_bfield, where perturbed b-field is computed, this is multiplied
!           by eta_prof = rho*(1-rho) so it vanishes at both ends. Multiply again
!           here by that factor to assure radial derivatives entering in dB also vanish
!           at endpoints s=0,1
         SCALE_TEST: DO irscale = 1, irad_scale
            CALL GetResPert(irscale, isign, mres0, nres0, rad, HelPert0, chip0, phip0)
            xc = 0
            wmhd = getwmhd(xc, PARFUNCTISL)

            IF (wmhd .LT. w1) THEN
               imin(1) = irscale
               imin(2) = isign
               w1 = wmhd
               p_width_min = p_width
            ENDIF
            END DO SCALE_TEST
         END DO PERT_SIGN

!     Make sure energy decreases
         wmhd = w1
!     recompute perturbation buv_rex here
         irscale = imin(1)
         isign   = imin(2)
         
         IF (iam .EQ. 0) THEN
            HP = HelPert0
            IF (isign .EQ. 2) HP = -HP

            DO iprint = 6, 33, 27
            WRITE (iprint, '(1p,e12.3,2(i8),3x,1pe10.2,0p,f8.2,f11.2,4x,2(f11.2))') &
                   10**6*(wmhd-w0), mres0, nres0, HP*normal, rad,         &
                   ABS(mres0*chip0+nres0*nfp*phip0), chip0/phip0, p_width_min
            END DO
         END IF
         
         CALL GetResPert(irscale, isign, mres0, nres0, rad, HelPert0, chip0, phip0)

         wmhd = getwmhd(xc, PARFUNCTISL)
         IF (ABS(wmhd-w1) .GT. 1.E-12_dp) THEN
            IF (iam .EQ. 0) PRINT *,'Error1 in Perturb'
            STOP 
         END IF
         xc = 0
         CALL update_state_par(.FALSE., zero, zero)
!Compute dW relative to perturbed state
         w0 = wmhd
   
         END DO JSURF_TEST                                                  ! Look for multiple resonances
         END DO NRES_TEST
      END DO RES_TEST

      IF (iam .EQ. 0) THEN
         DO iprint = 6, 33, 27
            WRITE (iprint, *)
         END DO
      END IF

      DEALLOCATE (buv_res, stat=istat)
      lresistive = lresist_save
      eta_factor = eta_factor_save

      END SUBROUTINE add_perturb


      LOGICAL FUNCTION FindResonance(mres0, nres0, resrad, jstart)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: mres0
      INTEGER, INTENT(in)     :: nres0
      INTEGER, INTENT(inout)   :: jstart
      REAL(rprec), INTENT(out) :: resrad
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                   :: js, jmin
      REAL(rprec)               :: del1, del2, rho1, rho2, delmin
!-----------------------------------------------
!
!     Starting at jstart, find place where mres0*chip+n*nfp*phip is minimum (changes sign)
!
      del2 = 0
      resrad = 0
      delmin = -1
      IF (mres0*nres0 .EQ. 0) THEN
          FindResonance=.FALSE.; RETURN
      END IF

!     Check for zero crossing
      DO js = jstart, ns-1
         del1 = mres0*chipf_i(js)+nres0*nfp*phipf_i(js)
         del1 = del1/MAX(ABS(mres0*chipf_i(js)),ABS(nres0*nfp*phipf_i(js)),1.E-20_dp)
         IF (delmin.EQ.-1 .OR. ABS(del1).LT.delmin) THEN
            IF (delmin .EQ. -1) del2 = del1
            delmin = ABS(del1)
            jmin = js
         END IF
         
         IF (del1*del2 < zero) THEN
            jstart = js+1
            rho2 = hs_i*(js-2)
            rho1 = hs_i*(js-1)
            resrad = (rho2*ABS(del1)+rho1*ABS(del2))/(ABS(del1)+ABS(del2))
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


      SUBROUTINE GetResPert(irscale, isign, mres0, nres0, rad,           &
                            HelPert0, chip0, phip0)
      USE fourier, ONLY: toijsp2, f_cos, f_sin
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)       :: irscale, isign
      INTEGER, INTENT(in)       :: mres0, nres0
      REAL(rprec), INTENT(in)   :: rad
      REAL(rprec), INTENT(out)  :: chip0, phip0
      REAL(rprec), INTENT(in)   :: HelPert0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: js, istat
      REAL(rprec) :: rho, rhores, pert_prof(ns), HelP_local, locrad
      REAL(rprec), ALLOCATABLE   :: bmn_res(:,:)
      REAL(rprec), ALLOCATABLE   :: temp_uv(:,:)
!-----------------------------------------------
      IF (mres0.GT.mpolin .OR. ABS(nres0).GT.ntorin) RETURN

      HelP_local = HelPert0
      IF (isign .EQ. 2) HelP_local = -HelPert0

      ALLOCATE (bmn_res(0:mpol,-ntor:ntor), stat=istat)
      IF (istat .NE. 0) STOP 'ISTAT != 0 IN GETRESPERT'
      bmn_res = 0

      js = INT(rad/hs_i) + 1
      IF (js.LT.1 .OR. js.GT.ns) RETURN
      chip0 = chipf_i(js)
      phip0 = phipf_i(js)
      locrad = REAL(irscale-1,dp)/(irad_scale-1)
      p_width = 0.1_dp + 0.9_dp*locrad                !decay length, ~delta-fcn for irscale=1 to broad (full radius)
      locrad = one/p_width
      DO js = 1,ns
         rho = hs_i*(js-1)
         rhores = rho-rad
         rhores = locrad*rhores
         pert_prof(js) = one/(one+rhores*rhores)      !tearing parity; mult by alpha for equal area
!         pert_prof(js) = rhores*pert_prof(js)        !odd-parity (rippling mode)
      END DO

      pert_prof = HelP_local*pert_prof
         
      bmn_res(mres0,nres0) = 1
      CALL toijsp2(bmn_res, buv_res, f_cos, 1)
      IF (lasym) THEN
         ALLOCATE (temp_uv(ntheta,nzeta), stat=istat)
         CALL toijsp2(bmn_res, temp_uv, f_sin, 1)
         buv_res(:,:,1) = buv_res(:,:,1) + temp_uv;
         DEALLOCATE(temp_uv)
      END IF
      DO js = 2,ns
         buv_res(:,:,js) = buv_res(:,:,1)*pert_prof(js)
      END DO

      buv_res(:,:,1) = pert_prof(1)*buv_res(:,:,1)

      DEALLOCATE (bmn_res, stat=istat)

      END SUBROUTINE GetResPert
        

      SUBROUTINE write_restart_file
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,  &
                            ns, mpol, ntor
      USE descriptor_mod, ONLY: iam
#if defined(MRC)
      USE restart_mod
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: istat, js, m, n
      CHARACTER*(256)   :: filename
!-----------------------------------------------
    
      IF (iam .NE. 0) RETURN
#if defined(MRC)
      istat = restart_write("siesta_" // TRIM(restart_ext) // ".nc", wout_file)
      RETURN
#elif defined(NETCDF)
      CALL write_restart_ncfile (istat)
      IF (istat .EQ. 0) RETURN
#endif

      filename = "siesta_" // TRIM(restart_ext) // ".txt"
      OPEN (unit=15, file=TRIM(filename), status='replace',              &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) STOP 'ERROR WRITING SIESTA RESTART FILE'

!     PRINT *,'WRITING RESTART FILE: ', TRIM(filename)
!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      WRITE (15,'(L1)') lasym
      WRITE (15,100) nsin, mpolin, ntorin
      WRITE (15,200) wout_file

      DO js=1,nsin
      DO m =0,mpolin
      DO n = -ntorin,ntorin
      WRITE (15,300) jbsupsmnsh(m,n,js), jbsupumnch(m,n,js),             &
     &               jbsupvmnch(m,n,js), jpmnch(m,n,js)
      ENDDO
      ENDDO
      ENDDO

      IF (lasym) THEN

      DO js=1,nsin
         DO m =0,mpolin
            DO n = -ntorin,ntorin
               WRITE (15,300) jbsupsmnch(m,n,js), jbsupumnsh(m,n,js),    &
                              jbsupvmnsh(m,n,js), jpmnsh(m,n,js)
            ENDDO
         ENDDO
      ENDDO

      ENDIF

      CLOSE (unit=15)

 100  FORMAT(3i8)
 200  FORMAT(a256)
 300  FORMAT(1p4e24.16)

      END SUBROUTINE write_restart_file
     
      
#if defined(NETCDF)
      SUBROUTINE write_restart_ncfile (istat)
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,  &
                            ns, mpol, ntor
      USE ezcdf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: ncdf, js, m, n
      CHARACTER*(256)   :: filename
!-----------------------------------------------
      filename = "siesta_" // TRIM(restart_ext) // ".nc"

      CALL cdf_open(ncdf,filename,'w',istat)
      IF (istat .NE. 0) RETURN

!     PRINT *,'WRITING RESTART FILE: ', TRIM(filename)

!================================
! Define Variables
!================================
!  Scalars 
      CALL cdf_define(ncdf, vn_lasym, lasym)
      CALL cdf_define(ncdf, vn_nsin, nsin)
      CALL cdf_define(ncdf, vn_mpolin, mpolin)
      CALL cdf_define(ncdf, vn_ntorin, ntorin)
      CALL cdf_define(ncdf, vn_wout, wout_file)

!  3D arrays
      CALL cdf_define(ncdf, vn_bsupss, jbsupsmnsh, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupuc, jbsupumnch, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupvc, jbsupvmnch, dimname=r3dim)
      CALL cdf_define(ncdf, vn_presc,  jpmnch,     dimname=r3dim)

      IF (lasym) THEN

      CALL cdf_define(ncdf, vn_bsupsc, jbsupsmnch, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupus, jbsupumnsh, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupvs, jbsupvmnsh, dimname=r3dim)
      CALL cdf_define(ncdf, vn_press,  jpmnsh,     dimname=r3dim)

      END IF
!================================
! Write Variables
!================================
      CALL cdf_write(ncdf, vn_lasym, lasym)
      CALL cdf_write(ncdf, vn_nsin, nsin)
      CALL cdf_write(ncdf, vn_mpolin, mpolin)
      CALL cdf_write(ncdf, vn_ntorin, ntorin)
      CALL cdf_write(ncdf, vn_wout, wout_file)

      CALL cdf_write(ncdf, vn_bsupss, jbsupsmnsh)
      CALL cdf_write(ncdf, vn_bsupuc, jbsupumnch)
      CALL cdf_write(ncdf, vn_bsupvc, jbsupvmnch)
      CALL cdf_write(ncdf, vn_presc,  jpmnch)

      IF (lasym) THEN
      CALL cdf_write(ncdf, vn_bsupsc, jbsupsmnch)
      CALL cdf_write(ncdf, vn_bsupus, jbsupumnsh)
      CALL cdf_write(ncdf, vn_bsupvs, jbsupvmnsh)
      CALL cdf_write(ncdf, vn_press,  jpmnsh)
      END IF

      CALL cdf_close(ncdf)

      END SUBROUTINE write_restart_ncfile
#endif

      SUBROUTINE read_restart_file (nprecon)
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch
#if defined(MRC)
      USE restart_mod
#endif
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nprecon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1
      INTEGER            :: istat, m, n, js, ns, ntor, mpol, nmin
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp
      REAL(rprec)        :: r0
      CHARACTER*(256)    :: filename
!-----------------------------------------------
#if defined(MRC)
      istat = restart_read("siesta_" // TRIM(restart_ext) // ".nc",      &
                           mpolin, ntorin, nsin, wout_file)
      RETURN
#elif defined(NETCDF)
      CALL read_restart_ncfile (istat, ns, ntor, mpol)
      nmin = MIN(ntor, ntorin)
      IF (istat .EQ. 0) GOTO 111
#endif

      filename = "siesta_" // TRIM(restart_ext) // ".txt"
      OPEN (unit=15, file=TRIM(filename), status='old',                  &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) THEN 
         PRINT *,'SIESTA RESTART FILE WAS NOT READ: STARTING AT BEGINNING'
         RETURN
      END IF

      IF (iam .EQ. 0) THEN
         PRINT *,'READING RESTART FILE: ', TRIM(filename)
      ENDIF
!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      READ (15,'(L1)') lasym
      READ (15,100) ns, mpol, ntor
      nmin = MIN(ntor, ntorin)

      ALLOCATE (bsupss(0:mpol,-ntor:ntor,ns),                            &
                bsupuc(0:mpol,-ntor:ntor,ns),                            &
                bsupvc(0:mpol,-ntor:ntor,ns),                            & 
                presc(0:mpol,-ntor:ntor,ns), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in read_restart_file'

      READ (15,200) wout_file

      DO js=1,ns
      DO m=0,mpol
      DO n=-ntor,ntor
      READ (15,300) bsupss(m,n,js), bsupuc(m,n,js),                      &
                    bsupvc(m,n,js), presc(m,n,js)
      ENDDO
      ENDDO
      ENDDO

      IF (lasym) THEN
      ALLOCATE (bsupsc(0:mpol,-ntor:ntor,ns),                            &
                bsupus(0:mpol,-ntor:ntor,ns),                            &
                bsupvs(0:mpol,-ntor:ntor,ns),                            & 
                press(0:mpol,-ntor:ntor,ns), stat=istat)
      DO js=1,ns
      DO m=0,mpol
      DO n=-ntor,ntor
      READ (15,300) bsupsc(m,n,js), bsupus(m,n,js),                      &
                    bsupvs(m,n,js), press(m,n,js)
      ENDDO
      ENDDO
      ENDDO
      END IF

      CLOSE (unit=15)

 111  CONTINUE

      jbsupsmnsh(:,:,1) = 0
      jbsupvmnch(:,:,1) = 0
      jpmnch(:,:,1) = 0

      CALL interpit(bsupss, jbsupsmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupsmnsh(m1,-nmin:nmin,1) = jbsupsmnsh(m1,-nmin:nmin,2)
      CALL interpit(bsupuc, jbsupumnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupumnch(m1,-nmin:nmin,1) = jbsupumnch(m1,-nmin:nmin,2)
      CALL interpit(bsupvc, jbsupvmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupvmnch(m0,-nmin:nmin,1) = jbsupvmnch(m0,-nmin:nmin,2)
      CALL interpit(presc, jpmnch,     ns, nsin, mpol, mpolin, ntor, ntorin)
      jpmnch(m0,-nmin:nmin,1) = jpmnch(m0,-nmin:nmin,2)
      DEALLOCATE (bsupss, bsupuc, bsupvc, presc, stat=istat)

      IF (lasym) THEN
      jbsupsmnch(:,:,1) = 0
      jbsupvmnsh(:,:,1) = 0
      jpmnsh(:,:,1) = 0

      CALL interpit(bsupsc, jbsupsmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupsmnch(m1,-nmin:nmin,1) = jbsupsmnch(m1,-nmin:nmin,2)
      CALL interpit(bsupus, jbsupumnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupumnsh(m1,-nmin:nmin,1) = jbsupumnsh(m1,-nmin:nmin,2)
      CALL interpit(bsupvs, jbsupvmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupvmnsh(m0,-nmin:nmin,1) = jbsupvmnsh(m0,-nmin:nmin,2)
      CALL interpit(press, jpmnsh,     ns, nsin, mpol, mpolin, ntor, ntorin)
      jpmnsh(m0,-nmin:nmin,1) = jpmnsh(m0,-nmin:nmin,2)
      DEALLOCATE (bsupsc, bsupus, bsupvs, press, stat=istat)
      END IF

      WRITE (33, 400) mpol, ntor, ns
      nprecon = 2   !USED TO SKIP DIAGONAL PRECONDITIONER IN EVOLVE

 100  FORMAT(3i8)
 200  FORMAT(a256)
 300  FORMAT(1p4e24.16)
 400  FORMAT(/,' RESTARTED FROM RUN PARAMETERS M: ',i3,' N: ',i3,' NS: ', i3,/)

      END SUBROUTINE read_restart_file

#if defined(NETCDF)
      SUBROUTINE read_restart_ncfile (istat, ns, ntor, mpol)
      USE ezcdf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: istat, ns, mpol, ntor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: ncdf
      CHARACTER*(256)   :: filename
!-----------------------------------------------
      filename = "siesta_" // TRIM(restart_ext) // ".nc"

      CALL cdf_open(ncdf,filename,'r',istat)
      IF (istat .NE. 0) RETURN

      CALL cdf_read(ncdf, vn_lasym, lasym)
      CALL cdf_read(ncdf, vn_nsin, ns)
      CALL cdf_read(ncdf, vn_mpolin, mpol)
      CALL cdf_read(ncdf, vn_ntorin, ntor)
      CALL cdf_read(ncdf, vn_wout, wout_file)

      ALLOCATE (bsupss(0:mpol,-ntor:ntor,ns),                            &
                bsupuc(0:mpol,-ntor:ntor,ns),                            &
                bsupvc(0:mpol,-ntor:ntor,ns),                            & 
                presc(0:mpol,-ntor:ntor,ns), stat=istat)
      IF (lasym) THEN
      ALLOCATE (bsupsc(0:mpol,-ntor:ntor,ns),                            &
                bsupus(0:mpol,-ntor:ntor,ns),                            &
                bsupvs(0:mpol,-ntor:ntor,ns),                            & 
                press(0:mpol,-ntor:ntor,ns), stat=istat)
      END IF
      IF (istat .NE. 0) STOP 'Memory allocation error in read_restart_ncfile'

!================================
! Read Variables
!================================

      CALL cdf_read(ncdf, vn_bsupss, bsupss)
      CALL cdf_read(ncdf, vn_bsupuc, bsupuc)
      CALL cdf_read(ncdf, vn_bsupvc, bsupvc)
      CALL cdf_read(ncdf, vn_presc,  presc)

      IF (lasym) THEN
      CALL cdf_read(ncdf, vn_bsupsc, bsupsc)
      CALL cdf_read(ncdf, vn_bsupus, bsupus)
      CALL cdf_read(ncdf, vn_bsupvs, bsupvs)
      CALL cdf_read(ncdf, vn_press,  press)
      END IF

      CALL cdf_close(ncdf)

      END SUBROUTINE read_restart_ncfile
#endif
  
#if !defined(MRC)
      SUBROUTINE interpit(aold, anew, ns_old, ns_new, m1, m2, n1, n2)
      USE stel_constants, ONLY: one
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: n1, n2, m1, m2, ns_old, ns_new
      REAL(rprec), INTENT(in)  :: aold(0:m1,-n1:n1,ns_old)
      REAL(rprec), INTENT(out) :: anew(0:m2,-n2:n2,ns_new)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER   :: p5 = 0.5_dp
      REAL(rprec), ALLOCATABLE :: temp1(:,:,:)
      INTEGER                  :: m12, n12, js_old, js_new, istat
      REAL(rprec)              :: hs_old, hs_new, snew, sold, x1
!-----------------------------------------------

      m12 = MIN(m1,m2);  n12 = MIN(n1, n2)
      
      anew = 0

      IF (ns_old .EQ. ns_new) THEN
         anew(0:m12,-n12:n12,:) = aold(0:m12,-n12:n12,:)
         RETURN
      END IF

      ALLOCATE (temp1(0:m12,-n12:n12,ns_old), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error #1 in INTERPIT'
!
!     FIRST INTERP TO FULL MESH FROM HALF-MESH
!     ASSUMES js=1 FULL-MESH VALUE VALUE IS STORED IN aold(1)
!
      temp1(:,:,1) = 0
      temp1(:,:,2:ns_old-1) = p5*(aold(:,:,2:ns_old-1) + aold(:,:,3:ns_old))
      temp1(:,:,ns_old) = 2*aold(:,:,ns_old) - temp1(:,:,ns_old-1)
      temp1(0:1,:,1)  = 2*aold(0:1,:,2) - temp1(0:1,:,2)
      hs_old = one/(ns_old-1)
      hs_new = one/(ns_new-1)

      DO js_new=1, ns_new
         snew = (js_new-1)*hs_new
         js_old = 1 + INT(snew/hs_old)
         x1 = snew - (js_old-1)*hs_old
         x1 = x1/hs_old
         IF (x1 .GT. one) x1 = 1
         IF (x1 .LT. zero) x1 = 0
         IF (js_old .gt. ns_old) STOP 'js_old >= ns_old in INTERPIT'
         IF (js_old .eq. ns_old) THEN
             anew(:,:,js_new) = temp1(:,:,js_old)
         ELSE
             anew(:,:,js_new) = (one-x1)*temp1(:,:,js_old)               &
                              + x1*temp1(:,:,js_old+1)
         END IF
      END DO
!      
!     AVERAGE BACK TO HALF-MESH
!
      DEALLOCATE (temp1)
      ALLOCATE (temp1(0:m12,-n12:n12,ns_new), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error #2 in INTERPIT'
      temp1 = anew
      anew(:,:,1) = 0
      anew(:,:,2:ns_new) = p5*(temp1(:,:,2:ns_new) + temp1(:,:,1:ns_new-1))
      DEALLOCATE(temp1)

      END SUBROUTINE interpit
#endif
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
      IF (istat .ne. 0) STOP 'Allocation error in Compute_Resonant_Epar'
      eparmn = 0
      epar = 0

      DO mres=0,mpol
         DO nres=-ntor,ntor
            jstart = 2
            DO WHILE (jstart .lt. ns-1)
               IF (.NOT.FindResonance(mres, nres, rad, jstart)) EXIT
               CALL AddResonantE (eparmn(:,mres,nres), rad, jstart)
            END DO
         END DO
      END DO

      DEALLOCATE(eparmn)

      END SUBROUTINE Compute_Resonant_Epar

      SUBROUTINE AddResonantE(eres, rad, jstart)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(inout) :: eres(ns)
      REAL(dp), INTENT(in)    :: rad
      INTEGER,  INTENT(in)    :: jstart
!-----------------------------------------------

      END SUBROUTINE AddResonantE


      END MODULE PERTURBATION
