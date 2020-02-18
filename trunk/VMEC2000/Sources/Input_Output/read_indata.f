!*******************************************************************************
!>  @file read_indata.f
!>  @brief Initializes VMEC after reading in the namelist input file.
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Read in the namelist input file and initialize VMEC.
!>
!>  @param[in]  in_file  Path to the namelist input file.
!>  @param[out] iunit    IO unit of the opened namelist input.
!>  @param[out] ier_flag Error status of the namelist read.
!-------------------------------------------------------------------------------
      SUBROUTINE read_indata(in_file, iunit, ier_flag)
      USE vmec_main
      USE vmec_params
      USE vacmod
      USE safe_open_mod
      USE parallel_include_module
      USE rho_bndy

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)    :: in_file
      INTEGER, INTENT(out)             :: ier_flag
      INTEGER, INTENT(out)             :: iunit

!  Local Variables
      INTEGER                          :: ireadseq
      INTEGER                          :: iosnml
      CLASS (rho_bndy_class), POINTER  :: bndy => null()

!  Start of executable code
      iunit = indata0 + rank
      CALL safe_open(iunit, ireadseq, in_file, 'old', 'formatted')
      IF (ireadseq .ne. 0) THEN
         WRITE (6,1000) TRIM(in_file), ireadseq
         ier_flag = input_error_flag
         RETURN
      ENDIF

      CALL read_namelist(iunit, iosnml, 'indata')
      IF (iosnml .ne. 0) THEN
         WRITE (6,1001) iosnml
         ier_flag = input_error_flag
         RETURN
      ENDIF

      CALL read_namelist(iunit, iosnml, 'mseprofile')

      IF (lrecon .and. itse.le.0 .and. imse.le.0) THEN
         lrecon = .false.
      END IF
      IF (lfreeb .and. mgrid_file.eq.'NONE') THEN
         lfreeb = .false.
      END IF

      IF (bloat .eq. zero) THEN
         bloat = one
      END IF
      IF (bloat .ne. one .and. ncurr .ne. 1) THEN
         ier_flag = 3
         RETURN
      ENDIF

!  Compute real space grid sizes.
      mpol = ABS(mpol)
      ntor = ABS(ntor)

!  FIXME: These should be asserts.
      IF (mpol .gt. mpold) THEN
         STOP 'mpol>mpold: lower mpol'
      END IF
      IF (ntor .gt. ntord) THEN
         STOP 'ntor>ntord: lower ntor'
      END IF

      mpol1 = mpol - 1
      ntor1 = ntor + 1

!  Number of theta grid points (>=2*mpol+6)
      IF (ntheta .le. 0) THEN
         ntheta = 2*mpol + 6
      END IF
      ntheta1 = 2*(ntheta/2)
      ntheta2 = 1 + ntheta1/2
      IF (ntor .eq. 0) THEN
         lthreed = .false.
      END IF
      IF (ntor .gt. 0) THEN
         lthreed = .true.
      END IF

!  Number of zeta grid points (=1 IF ntor=0)
      IF (ntor.eq.0 .and. nzeta .eq. 0) THEN
         nzeta = 1
      END IF
      IF (nzeta .le. 0) THEN
         nzeta = 2*ntor + 4
      END IF

!  Size of rmnc, zmns, ...
      mnmax = ntor1 + mpol1*(1 + 2*ntor)
!  Size of rmncc, rmnss, ...
      mnsize = mpol*ntor1

      mf = mpol + 1
      nf = ntor
      nu = ntheta1
      nv = nzeta
      mf1 = 1 + mf
      nf1 = 2*nf + 1
      mnpd = mf1*nf1
      nfper = nfp

!  Set rb* zb* values from rho boundary.
      IF (ANY(rhoc .ne. 0)) THEN
         bndy => rho_bndy_class(mpol)
         CALL bndy%to_rz(mpol, ntor, rhoc(-ntor:ntor,0:mpol-1),                &
     &                   raxis_cc(0:ntor), zaxis_cs(0:ntor),                   &
     &                   rbc(-ntor:ntor,0:mpol-1),                             &
     &                   zbs(-ntor:ntor,0:mpol-1), rho_bndy_sym)
         IF (lasym) THEN
            CALL bndy%to_rz(mpol, ntor, rhos(-ntor:ntor,0:mpol-1),             &
     &                      raxis_cs(0:ntor), zaxis_cc(0:ntor),                &
     &                      rbs(-ntor:ntor,0:mpol-1),                          &
     &                      zbc(-ntor:ntor,0:mpol-1), rho_bndy_asym)
         END IF
         DEALLOCATE(bndy)
      END IF

!  Setup indexing for packed-array structure of xc and gc arrays.
      rcc = 1
      zsc = 1
      rss = 0
      rsc = 0
      rcs = 0
      zcc = 0
      zss = 0
      zcs = 0
      IF (.not.lasym) THEN
         ntheta3 = ntheta2
         mnpd2 = mnpd
         IF (lthreed) THEN
            ntmax = 2
            rss = 2
            zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         ntheta3 = ntheta1
         mnpd2 = 2*mnpd
         IF (lthreed) THEN
             ntmax = 4
             rss = 2
             rsc = 3
             rcs = 4
             zcs = 2
             zcc = 3
             zss = 4
         ELSE
             ntmax = 2
             rsc = 2
             zcc = 2
         END IF
      END IF

      nuv = nu*nv
      nu2 = nu/2 + 1
      nu3 = ntheta3
      nznt = nzeta*ntheta3
      nuv3 = nznt

!  Obsolete old format. May not be reading in ac.
!!Old FORMAT: may not be reading in ac
      IF (ncurr .eq. 1 .and. ALL(ac .eq. cbig)) THEN
         ac = ai
      END IF
      WHERE (ac .eq. cbig)
          ac = zero
      END WHERE

1000  FORMAT(' In VMEC, error opening input file: ',a,'. Iostat = ',i4)
1001  FORMAT(' In VMEC, indata NAMELIST error: iostat = ',i4)

      END SUBROUTINE
