!*******************************************************************************
!>  @file restart_mod.f90
!>  @brief Contains module @ref restart_mod.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Contains routines for writting the restart file.
!*******************************************************************************
      MODULE restart_mod
      USE ezcdf
      USE stel_kinds
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE metrics, ONLY: tolowerh
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: lasym, lrecon, unit_out
      USE v3_utilities, ONLY: assert_eq
      USE stel_constants, ONLY: mu0

      IMPLICIT NONE

!*******************************************************************************
!  restart module parameters
!*******************************************************************************
!>  Version number.
      INTEGER, PARAMETER           :: restart_version = 0

!  Flag will be placed in the last bits. This way the version number and flags 
!  can fit in the same memory.
!>  Bit position for the lasym flag.
      INTEGER, PARAMETER           :: restart_lasym = 31
!>  Bit position for the lrecon flag.
      INTEGER, PARAMETER           :: restart_lrecon = 30

!>  Radial Dimension names.
      CHARACTER (LEN=*), DIMENSION(1), PARAMETER ::                            &
     &   radial_dim = (/ 'radius' /)
!>  Fourier Dimension names.
      CHARACTER (LEN=*), DIMENSION(3), PARAMETER ::                            &
     &   restart_dims = (/ 'm-mode', 'n-mode', radial_dim(1) /)

!>  Name for the restart file number of radial points.
      CHARACTER (len=*), PARAMETER :: vn_nsin = 'nrad'
!>  Name for the restart file number of poloidal modes.
      CHARACTER (len=*), PARAMETER :: vn_mpolin = 'mpol'
!>  Name for the restart file number of toroidal modes.
      CHARACTER (len=*), PARAMETER :: vn_ntorin = 'ntor'
!>  Name for the restart file number of field periods.
      CHARACTER (len=*), PARAMETER :: vn_nfpin = 'nfp'
!>  Name for the restart file number of wout file modes.
      CHARACTER (len=*), PARAMETER :: vn_wout = 'wout_file'
!>  Name for the restart file number of state flags modes.
      CHARACTER (len=*), PARAMETER :: vn_flags = 'state_flags'

!>  Name for the restart file jbsupss.
      CHARACTER (len=*), PARAMETER :: vn_jbsupss = 'JBsupssh(m,n,r)'
!>  Name for the restart file jbsupuc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupuc = 'JBsupuch(m,n,r)'
!>  Name for the restart file jbsupvc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupvc = 'JBsupvch(m,n,r)'
!>  Name for the restart file jbsupsc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupsc = 'JBsupsch(m,n,r)'
!>  Name for the restart file jbsupus.
      CHARACTER (len=*), PARAMETER :: vn_jbsupus = 'JBsupush(m,n,r)'
!>  Name for the restart file jbsupus.
      CHARACTER (len=*), PARAMETER :: vn_jbsupvs = 'JBsupvsh(m,n,r)'
!>  Name for the restart file jpresc.
      CHARACTER (len=*), PARAMETER :: vn_jpresc = 'jpresch(m,n,r)'
!>  Name for the restart file jpress.
      CHARACTER (len=*), PARAMETER :: vn_jpress = 'jpressh(m,n,r)'

!>  Name for the restart file rmnc.
      CHARACTER (len=*), PARAMETER :: vn_rmnc = 'rmnc(m,n,r)'
!>  Name for the restart file rmns.
      CHARACTER (len=*), PARAMETER :: vn_rmns = 'rmns(m,n,r)'
!>  Name for the restart file zmnc.
      CHARACTER (len=*), PARAMETER :: vn_zmnc = 'zmnc(m,n,r)'
!>  Name for the restart file zmns.
      CHARACTER (len=*), PARAMETER :: vn_zmns = 'zmns(m,n,r)'

!>  Name for the restart file chipf.
      CHARACTER (len=*), PARAMETER :: vn_chipf = 'chipf(r)'
!>  Name for the restart file phipf.
      CHARACTER (len=*), PARAMETER :: vn_phipf = 'phipf(r)'

!>  Name for the restart file bsupsmns.
      CHARACTER (len=*), PARAMETER :: vn_bsupsmns = 'bsupsmnsh(m,n,r)'
!>  Name for the restart file bsupsmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupsmnc = 'bsupsmnch(m,n,r)'
!>  Name for the restart file bsupumns.
      CHARACTER (len=*), PARAMETER :: vn_bsupumns = 'bsupumnsh(m,n,r)'
!>  Name for the restart file bsupumnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupumnc = 'bsupumnch(m,n,r)'
!>  Name for the restart file bsupvmns.
      CHARACTER (len=*), PARAMETER :: vn_bsupvmns = 'bsupvmnsh(m,n,r)'
!>  Name for the restart file bsupvmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupvmnc = 'bsupvmnch(m,n,r)'
!>  Name for the restart file bsubsmns.
      CHARACTER (len=*), PARAMETER :: vn_bsubsmns = 'bsubsmnsh(m,n,r)'
!>  Name for the restart file bsubsmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubsmnc = 'bsubsmnch(m,n,r)'
!>  Name for the restart file bsubumns.
      CHARACTER (len=*), PARAMETER :: vn_bsubumns = 'bsubumnsh(m,n,r)'
!>  Name for the restart file bsubumnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubumnc = 'bsubumnch(m,n,r)'
!>  Name for the restart file bsubvmns.
      CHARACTER (len=*), PARAMETER :: vn_bsubvmns = 'bsubvmnsh(m,n,r)'
!>  Name for the restart file bsubvmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubvmnc = 'bsubvmnch(m,n,r)'
!>  Name for the restart file pmns.
      CHARACTER (len=*), PARAMETER :: vn_pmns = 'pmnsh(m,n,r)'
!>  Name for the restart file pmnc.
      CHARACTER (len=*), PARAMETER :: vn_pmnc = 'pmnch(m,n,r)'
!>  Name for the restart file jksupsmns.
      CHARACTER (len=*), PARAMETER :: vn_jksupsmns = 'jksupsmnsf(m,n,r)'
!>  Name for the restart file jksupsmnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupsmnc = 'jksupsmncf(m,n,r)'
!>  Name for the restart file jksupumns.
      CHARACTER (len=*), PARAMETER :: vn_jksupumns = 'jksupumnsf(m,n,r)'
!>  Name for the restart file jksupumnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupumnc = 'jksupumncf(m,n,r)'
!>  Name for the restart file jksupvmns.
      CHARACTER (len=*), PARAMETER :: vn_jksupvmns = 'jksupvmnsf(m,n,r)'
!>  Name for the restart file jksupvmnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupvmnc = 'jksupvmncf(m,n,r)'

!>  Name for the restart file p_factor.
      CHARACTER (len=*), PARAMETER :: vn_p_factor = 'p_factor'
!>  Name for the restart file b_factor.
      CHARACTER (len=*), PARAMETER :: vn_b_factor = 'b_factor'
!>  Name for the restart file wb.
      CHARACTER (len=*), PARAMETER :: vn_wb = 'wb'
!>  Name for the restart file wb.
      CHARACTER (len=*), PARAMETER :: vn_wp = 'wp'
!>  Name for the restart file rmajor.
      CHARACTER (len=*), PARAMETER :: vn_rmajor = 'rmajor'

!>  Name for the restart file p_max.
      CHARACTER (len=*), PARAMETER :: vn_p_max = 'p_max'
!>  Name for the restart file p_min.
      CHARACTER (len=*), PARAMETER :: vn_p_min = 'p_min'

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the restart file.
!>
!>  Reads the restart information and initalizes SIESTA quantities.
!>
!>  @param[in]    restart_ext Restart extention.
!>  @param[inout] wout_file   Name of the wout file.
!>  @param[in]    mpolin      Namelist number of polodal modes.
!>  @param[in]    ntorin      Namelist number of toroidal modes.
!>  @param[in]    nsin        Namelist number of radial grid points.
!>  @param[in]    nfpin       Namelist number of field periods.
!>  @returns Preconditioner control flag.
!-------------------------------------------------------------------------------
      FUNCTION restart_read(restart_ext, wout_file, mpolin, ntorin, nfpin,     &
                            nsin)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch,                                &
                            b_factor, p_factor, alloc_quantities
      USE island_params, ONLY: chipf => chipf_i, phipf => phipf_i,             &
     &                         wb => wb_i, wp => wp_i, nfp_i, gamma,           &
     &                         gnorm => gnorm_i, rmajor => rmajor_i
      USE vmec_info, ONLY: rmnc => rmnc_i, zmns => zmns_i,                     &
     &                     rmns => rmns_i, zmnc => zmnc_i,                     &
     &                     vmec_info_construct_island
      USE metrics, ONLY: LoadGrid
      USE fourier, ONLY: orthonorm, init_fourier, m0, m1

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: restart_read
      CHARACTER (len=*), INTENT(in)            :: restart_ext
      CHARACTER (len=*), INTENT(inout)         :: wout_file
      INTEGER, INTENT(in)                      :: mpolin
      INTEGER, INTENT(in)                      :: ntorin
      INTEGER, INTENT(in)                      :: nsin
      INTEGER, INTENT(in)                      :: nfpin

!  local variables
      INTEGER                                  :: flags
      INTEGER                                  :: ncid
      INTEGER                                  :: ns
      INTEGER                                  :: mpol
      INTEGER                                  :: ntor
      INTEGER                                  :: nfp
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: tempmn_r
      REAL (dp), DIMENSION(:), ALLOCATABLE     :: temp_r
      INTEGER                                  :: nmin
      INTEGER                                  :: status
      CHARACTER (LEN=256)                      :: filename
      INTEGER                                  :: s

!  Start of executable code
      restart_read = 1

      CALL alloc_quantities

      filename = 'siesta_' // TRIM(restart_ext) // '.nc'
      CALL cdf_open(ncid, TRIM(filename), 'r', status)
      CALL assert_eq(0, status, 'Failed to read restart file: ' //             &
     &                          TRIM(filename))

      CALL cdf_read(ncid, vn_flags, flags)

!  The namelist inut file can change the size of the radial grid and the 
!  poloidal and toroidal modes.
      CALL cdf_read(ncid, vn_nsin, ns)
      CALL cdf_read(ncid, vn_mpolin, mpol)
      CALL cdf_read(ncid, vn_ntorin, ntor)
      CALL cdf_read(ncid, vn_nfpin, nfp)

      IF (nfpin .lt. 1) THEN
         nfp_i = nfp
      ELSE
         nfp_i = nfpin
      END IF

      nmin = MIN(ntor, ntorin)

      CALL cdf_read(ncid, vn_wout, wout_file)

!  The namelist input file may turn the asymmetric terms on and off.
      ALLOCATE(tempmn_r(0:mpol,-ntor:ntor,ns))
      ALLOCATE(temp_r(ns))
      CALL vmec_info_construct_island(mpolin, ntorin, nsin, lasym)

      ALLOCATE(chipf(nsin))
      CALL cdf_read(ncid, vn_chipf, temp_r)
      CALL interpit_1d(temp_r, chipf, ns, nsin, .false., 1)

      ALLOCATE(phipf(nsin))
      CALL cdf_read(ncid, vn_phipf, temp_r)
      CALL interpit_1d(temp_r, phipf, ns, nsin, .false., 1)

      CALL cdf_read(ncid, vn_jbsupss, tempmn_r)
      jbsupsmnsh(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupsmnsh, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfpin, .true.)
      jbsupsmnsh(m1,-nmin:nmin,1) = jbsupsmnsh(m1,-nmin:nmin,2)

      CALL cdf_read(ncid, vn_jbsupuc, tempmn_r)
      jbsupumnch(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupumnch, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfpin, .true.)
      jbsupumnch(m1,-nmin:nmin,1) = jbsupumnch(m1,-nmin:nmin,2)

      CALL cdf_read(ncid, vn_jbsupvc, tempmn_r)
      jbsupvmnch(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupvmnch, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfpin, .true.)
      jbsupvmnch(m0,-nmin:nmin,1) = jbsupvmnch(m0,-nmin:nmin,2)

      CALL cdf_read(ncid, vn_jpresc,  tempmn_r)
      jpmnch(:,:,1) = 0
      CALL interpit(tempmn_r, jpmnch,     ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfpin, .true.)
      jpmnch(m0,-nmin:nmin,1) = jpmnch(m0,-nmin:nmin,2)

      CALL cdf_read(ncid, vn_rmnc, tempmn_r)
      CALL interpit(tempmn_r, rmnc, ns, nsin, mpol, mpolin,                    &
                    ntor, ntorin, nfp, nfpin, .false.)

      CALL cdf_read(ncid, vn_zmns, tempmn_r)
      CALL interpit(tempmn_r, zmns, ns, nsin, mpol, mpolin,                    &
                    ntor, ntorin, nfp, nfpin, .false.)

      IF (BTEST(flags, restart_lasym) .and. lasym) THEN
         jbsupsmnch(:,:,1) = 0
         CALL cdf_read(ncid, vn_jbsupsc, tempmn_r)
         CALL interpit(tempmn_r, jbsupsmnch, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfpin, .true.)
         jbsupsmnch(m1,-nmin:nmin,1) = jbsupsmnch(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupus, tempmn_r)
         CALL interpit(tempmn_r, jbsupumnsh, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfpin, .true.)
         jbsupumnsh(m1,-nmin:nmin,1) = jbsupumnsh(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupvs, tempmn_r)
         CALL interpit(tempmn_r, jbsupvmnsh, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfpin, .true.)
         jbsupvmnsh(m0,-nmin:nmin,1) = jbsupvmnsh(m0,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jpress,  tempmn_r)
         CALL interpit(tempmn_r, jpmnsh,     ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfpin, .true.)
         jpmnsh(m0,-nmin:nmin,1) = jpmnsh(m0,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_rmns, tempmn_r)
         CALL interpit(tempmn_r, rmns, ns, nsin, mpol, mpolin,                 &
                       ntor, ntorin, nfp, nfpin, .false.)

         CALL cdf_read(ncid, vn_zmns, tempmn_r)
         CALL interpit(tempmn_r, zmnc, ns, nsin, mpol, mpolin,                 &
                       ntor, ntorin, nfp, nfpin, .false.)

      ELSE IF (lasym) THEN
         jbsupsmnch = 0
         jbsupumnsh = 0
         jbsupvmnsh = 0
         jpmnsh = 0
         rmns = 0
         rmnc = 0
      END IF

!  Read normalization factors.
      CALL cdf_read(ncid, vn_wb, wb)
      CALL cdf_read(ncid, vn_wp, wp)
      CALL cdf_read(ncid, vn_rmajor, rmajor)

      CALL cdf_close(ncid)

      CALL init_fourier

      DO s = 1, nsin
         rmnc(:,:,s) = rmnc(:,:,s)/orthonorm
         zmns(:,:,s) = zmns(:,:,s)/orthonorm
         IF (lasym) THEN
            rmns(:,:,s) = rmns(:,:,s)/orthonorm
            zmnc(:,:,s) = zmnc(:,:,s)/orthonorm
         END IF
      END DO

      CALL LoadGrid(status)

      DEALLOCATE(tempmn_r)
      DEALLOCATE(temp_r)

!  Init quantities.
      p_factor = 1.0/ABS(wb)
      b_factor = SQRT(p_factor)
      gnorm = ABS(wb)/(wb + wp/(gamma - 1.0))

      WRITE (unit_out, 1000) mpol, ntor, ns

      restart_read = 2

1000  FORMAT(/,' RESTARTED FROM RUN PARAMETERS M: ',i3,' N: ',i3,' NS: ', i3,/)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write the restart file.
!>
!>  Writes the restart information.
!>
!>  @param[in] restart_ext Restart file extension.
!>  @param[in] wout_file   Name of the wout file.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_write(restart_ext, wout_file)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch,                                &
                            b_factor,   p_factor, jacobh
      USE fourier
      USE island_params, ONLY: nfp => nfp_i, chipf => chipf_i,                 &
                               phipf => phipf_i, wb => wb_i, wp => wp_i,       &
                               rmajor => rmajor_i
      USE shared_data, ONLY: r1_i, z1_i
      USE utilities, ONLY: set_bndy_half_to_full, set_bndy_half_to_full_ds

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)            :: restart_ext
      CHARACTER (len=*), INTENT(in)            :: wout_file

!  local variables
      INTEGER                                  :: flags
      INTEGER                                  :: ncid
      INTEGER                                  :: nblock
      INTEGER                                  :: s
      INTEGER                                  :: m
      INTEGER                                  :: n
      INTEGER                                  :: status

      REAL (dp)                                :: r0
      REAL (dp)                                :: tempscalar
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: dbuds
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: dbvds
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: tempmn_w
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pijh

!  Start of executable code
      CALL cdf_open(ncid, 'siesta_' // TRIM(restart_ext) // '.nc', 'w', status)
      CALL assert_eq(0, status, 'Failed to write restart file.')

      flags = restart_version
      IF (lasym) THEN
         flags = IBSET(flags, restart_lasym)
      END IF
      IF (lrecon) THEN
         flags = IBSET(flags, restart_lrecon)
      END IF

      ALLOCATE(tempmn_w(0:mpol,-ntor:ntor,ns))

      CALL cdf_define(ncid, vn_flags, flags)
      CALL cdf_define(ncid, vn_nsin, ns)
      CALL cdf_define(ncid, vn_mpolin, mpol)
      CALL cdf_define(ncid, vn_ntorin, ntor)
      CALL cdf_define(ncid, vn_nfpin, nfp)
      CALL cdf_define(ncid, vn_wout, wout_file)

      CALL cdf_define(ncid, vn_jbsupss, jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupuc, jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupvc, jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jpresc,  jpmnch,     dimname=restart_dims)

      CALL cdf_define(ncid, vn_rmnc, tempmn_w, dimname=restart_dims)
      CALL cdf_define(ncid, vn_zmns, tempmn_w, dimname=restart_dims)

      CALL cdf_define(ncid, vn_chipf, chipf, dimname=radial_dim)
      CALL cdf_define(ncid, vn_phipf, phipf, dimname=radial_dim)

      CALL cdf_define(ncid, vn_wb, wb)
      CALL cdf_define(ncid, vn_wp, wp)

      IF (lasym) THEN
         CALL cdf_define(ncid, vn_jbsupsc, jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupus, jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupvs, jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jpress,  jpmnsh,     dimname=restart_dims)

         CALL cdf_define(ncid, vn_rmns, tempmn_w, dimname=restart_dims)
         CALL cdf_define(ncid, vn_zmnc, tempmn_w, dimname=restart_dims)
      END IF

      RECON0: IF (lrecon) THEN
!  Using variable with jacobian only for the dimension sizes and types. Before
!  writting, the jacobian normalization will be removed.
         CALL cdf_define(ncid, vn_bsupsmns,  jbsupsmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsupumnc,  jbsupumnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsupvmnc,  jbsupvmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubsmns,  jbsupsmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubumnc,  jbsupumnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubvmnc,  jbsupvmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupsmns, jbsupsmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupumnc, jbsupumnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupvmnc, jbsupvmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_pmnc,      jpmnch,     dimname=restart_dims)
         CALL cdf_define(ncid, vn_b_factor,  b_factor)
         CALL cdf_define(ncid, vn_p_factor,  p_factor)
         CALL cdf_define(ncid, vn_rmajor,    rmajor)
         CALL cdf_define(ncid, vn_p_max,     b_factor)
         CALL cdf_define(ncid, vn_p_min,     p_factor)
         IF (lasym) THEN
            CALL cdf_define(ncid, vn_bsupsmnc, jbsupsmnch, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsupumns, jbsupumnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsupvmns, jbsupvmnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubsmnc, jbsupsmnch, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubumns, jbsupumnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubvmns, jbsupvmnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupsmnc, jbsupsmnch,                    &
                            dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupumns, jbsupumnsh,                    &
                            dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupvmns, jbsupvmnsh,                    &
                            dimname=restart_dims)
            CALL cdf_define(ncid, vn_pmns, jpmnsh, dimname=restart_dims)
         END IF
      END IF RECON0

      CALL cdf_write(ncid, vn_flags, flags)
      CALL cdf_write(ncid, vn_nsin, ns)
      CALL cdf_write(ncid, vn_mpolin, mpol)
      CALL cdf_write(ncid, vn_ntorin, ntor)
      CALL cdf_write(ncid, vn_nfpin, nfp)
      CALL cdf_write(ncid, vn_wout, wout_file)

      CALL cdf_write(ncid, vn_jbsupss, jbsupsmnsh)
      CALL cdf_write(ncid, vn_jbsupuc, jbsupumnch)
      CALL cdf_write(ncid, vn_jbsupvc, jbsupvmnch)
      CALL cdf_write(ncid, vn_jpresc,  jpmnch)

      CALL cdf_write(ncid, vn_chipf, chipf)
      CALL cdf_write(ncid, vn_phipf, phipf)

      CALL tomnsp(r1_i, tempmn_w, f_cos)
      DO s = 1, ns
         tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
      END DO
      CALL cdf_write(ncid, vn_rmnc, tempmn_w)
      CALL tomnsp(z1_i, tempmn_w, f_sin)
      DO s = 1, ns
         tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
      END DO
      CALL cdf_write(ncid, vn_zmns, tempmn_w)

      IF (lasym) THEN
         CALL cdf_write(ncid, vn_jbsupsc, jbsupsmnch)
         CALL cdf_write(ncid, vn_jbsupus, jbsupumnsh)
         CALL cdf_write(ncid, vn_jbsupvs, jbsupvmnsh)
         CALL cdf_write(ncid, vn_jpress,  jpmnsh)

         CALL tomnsp(r1_i, tempmn_w, f_sin)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_rmns, tempmn_w)
         CALL tomnsp(z1_i, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_zmnc, tempmn_w)
      END IF

      CALL cdf_write(ncid, vn_wb, wb)
      CALL cdf_write(ncid, vn_wp, wp)
      CALL cdf_write(ncid, vn_rmajor, rmajor)

      RECON1: IF (lrecon) THEN
         ALLOCATE(bsupsijh(ntheta,nzeta,ns))
         ALLOCATE(bsupuijh(ntheta,nzeta,ns))
         ALLOCATE(bsupvijh(ntheta,nzeta,ns))
         ALLOCATE(pijh(ntheta,nzeta,ns))

         bsupsijh = 0
         bsupuijh = 0
         bsupvijh = 0
         pijh = 0

         CALL toijsp(jbsupsmnsh, bsupsijh, f_none, f_sin)
         CALL toijsp(jbsupumnch, bsupuijh, f_none, f_cos)
         CALL toijsp(jbsupvmnch, bsupvijh, f_none, f_cos)
         CALL toijsp(jpmnch, pijh, f_none, f_cos)

         CALL cdf_write(ncid, vn_p_factor, p_factor)
         CALL cdf_write(ncid, vn_b_factor, b_factor)

         IF (lasym) THEN
            CALL toijsp(jbsupsmnch, bsupsijh, f_sum, f_cos)
            CALL toijsp(jbsupumnsh, bsupuijh, f_sum, f_sin)
            CALL toijsp(jbsupvmnsh, bsupvijh, f_sum, f_sin)
            CALL toijsp(jpmnsh, pijh, f_sum, f_sin)
         END IF

         bsupsijh(:,:,1) = 0
         bsupuijh(:,:,1) = 0
         bsupvijh(:,:,1) = 0
         pijh(:,:,1) = 0

!  Remove the jacobian term.
         bsupsijh = bsupsijh/jacobh
         bsupuijh = bsupuijh/jacobh
         bsupvijh = bsupvijh/jacobh
         pijh = pijh/jacobh

         tempscalar = MAXVAL(pijh(:,:,2:))
         CALL cdf_write(ncid, vn_p_max, tempscalar)
         tempscalar = MINVAL(pijh(:,:,2:))
         CALL cdf_write(ncid, vn_p_min, tempscalar)

         tempmn_w = 0

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
         CALL tomnsp(bsupsijh, tempmn_w, f_sin)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsupsmns, tempmn_w)

         CALL tomnsp(bsupuijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsupumnc, tempmn_w)

         CALL tomnsp(bsupvijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsupvmnc, tempmn_w)

!  Remove the orthonorm and p_factor so this quantity can be directly summed.
         CALL tomnsp(pijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)     &
                            / (p_factor*mu0)
         END DO
         CALL cdf_write(ncid, vn_pmnc, tempmn_w)

         ALLOCATE(bsubsijh(ntheta,nzeta,ns))
         ALLOCATE(bsubuijh(ntheta,nzeta,ns))
         ALLOCATE(bsubvijh(ntheta,nzeta,ns))

         CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                           &
                       bsubsijh, bsubuijh, bsubvijh, 1, ns)

!  Need these to compute the currents later.
         ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnh(0:mpol,-ntor:ntor,ns))

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
         CALL tomnsp(bsubsijh, bsubsmnh, f_sin)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubsmnh(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsubsmns, tempmn_w)

         CALL tomnsp(bsubuijh, bsubumnh, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubumnh(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsubumnc, tempmn_w)

         CALL tomnsp(bsubvijh, bsubvmnh, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubvmnh(:,:,s)     &
                            / b_factor
         END DO
         CALL cdf_write(ncid, vn_bsubvmnc, tempmn_w)

!  Need full mesh values to compute currents on the full mesh.
         ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnf(0:mpol,-ntor:ntor,ns))

         CALL to_full_mesh(bsubsmnh, bsubsmnf)
         CALL to_full_mesh(bsubumnh, bsubumnf)
         CALL to_full_mesh(bsubvmnh, bsubvmnf)

         CALL set_bndy_half_to_full(bsubsmnh, bsubsmnh, 1, f_sin)
         CALL set_bndy_half_to_full(bsubumnh, bsubumnh, 1, f_cos)
         CALL set_bndy_half_to_full(bsubvmnh, bsubvmnh, 1, f_cos)

         ALLOCATE(jksupsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupvmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbuds(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbvds(0:mpol,-ntor:ntor,ns))

!  Half grid radial derivatives become full grid.
         CALL GradientFull(dbuds, bsubumnh)
         CALL GradientFull(dbvds, bsubvmnh)

         CALL set_bndy_half_to_full_ds(bsubumnh, dbuds, 1, f_cos)
         CALL set_bndy_half_to_full_ds(bsubvmnh, dbuds, 1, f_cos)

         DO n = -ntor, ntor
            DO m = 0, mpol
               jksupsmnf(m,n,:) = -m*    bsubvmnf(m,n,:) + n*nfp*bsubumnf(m,n,:)
               jksupumnf(m,n,:) =  n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
               jksupvmnf(m,n,:) =  dbuds(m,n,:)          - m*    bsubsmnf(m,n,:)
               jksupsmnf(m,n,:) =  orthonorm(m,n)*jksupsmnf(m,n,:)/b_factor
               jksupumnf(m,n,:) =  orthonorm(m,n)*jksupumnf(m,n,:)/b_factor
               jksupvmnf(m,n,:) =  orthonorm(m,n)*jksupvmnf(m,n,:)/b_factor
            END DO
         END DO

         CALL cdf_write(ncid, vn_jksupsmns, jksupsmnf)
         CALL cdf_write(ncid, vn_jksupumnc, jksupumnf)
         CALL cdf_write(ncid, vn_jksupvmnc, jksupvmnf)

         IF (lasym) THEN
!  Remove the orthonorm and b_factor so these quantities can be directly summed.
            CALL tomnsp(bsupsijh, tempmn_w, f_cos)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsupsmnc, tempmn_w)

            CALL tomnsp(bsupuijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsupumns, tempmn_w)

            CALL tomnsp(bsupvijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsupvmns, tempmn_w)

!  Remove the orthonorm and p_factor so this quantity can be directly summed.
            CALL tomnsp(pijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)  &
                               / (p_factor*mu0)
            END DO
            CALL cdf_write(ncid, vn_pmns, tempmn_w)

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
            CALL tomnsp(bsubsijh, bsubsmnh, f_cos)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubsmnh(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsubsmnc, tempmn_w)

            CALL tomnsp(bsubuijh, bsubumnh, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubumnh(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsubumns, tempmn_w)

            CALL tomnsp(bsubvijh, bsubvmnh, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubvmnh(:,:,s)  &
                               / b_factor
            END DO
            CALL cdf_write(ncid, vn_bsubvmns, tempmn_w)

            CALL to_full_mesh(bsubsmnh, bsubsmnf)
            CALL to_full_mesh(bsubumnh, bsubumnf)
            CALL to_full_mesh(bsubvmnh, bsubvmnf)

            CALL set_bndy_half_to_full(bsubsmnh, bsubsmnh, 1, f_cos)
            CALL set_bndy_half_to_full(bsubumnh, bsubumnh, 1, f_sin)
            CALL set_bndy_half_to_full(bsubvmnh, bsubvmnh, 1, f_sin)

            CALL GradientFull(dbuds, bsubumnh)
            CALL GradientFull(dbvds, bsubvmnh)

            CALL set_bndy_half_to_full_ds(bsubumnh, dbuds, 1, f_sin)
            CALL set_bndy_half_to_full_ds(bsubvmnh, dbuds, 1, f_sin)

            DO n = -ntor, ntor
               DO m = 0, mpol
                  jksupsmnf(m,n,:) = m*bsubvmnf(m,n,:) - n*nfp*bsubumnf(m,n,:)
                  jksupumnf(m,n,:) = -n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
                  jksupvmnf(m,n,:) = dbuds(m,n,:) + m*bsubsmnf(m,n,:)
                  jksupsmnf(m,n,:) = orthonorm(m,n)*jksupsmnf(m,n,:)/b_factor
                  jksupumnf(m,n,:) = orthonorm(m,n)*jksupumnf(m,n,:)/b_factor
                  jksupvmnf(m,n,:) = orthonorm(m,n)*jksupvmnf(m,n,:)/b_factor
               END DO
            END DO

            CALL cdf_write(ncid, vn_jksupsmnc, jksupsmnf)
            CALL cdf_write(ncid, vn_jksupumns, jksupumnf)
            CALL cdf_write(ncid, vn_jksupvmns, jksupvmnf)
         END IF

         DEALLOCATE(bsupsijh)
         DEALLOCATE(bsupuijh)
         DEALLOCATE(bsupvijh)
         DEALLOCATE(pijh)

         DEALLOCATE(bsubsmnh)
         DEALLOCATE(bsubumnh)
         DEALLOCATE(bsubvmnh)

         DEALLOCATE(bsubsmnf)
         DEALLOCATE(bsubumnf)
         DEALLOCATE(bsubvmnf)

         DEALLOCATE(jksupsmnf)
         DEALLOCATE(jksupumnf)
         DEALLOCATE(jksupvmnf)
         DEALLOCATE(dbuds)
         DEALLOCATE(dbvds)

      END IF RECON1

      DEALLOCATE(tempmn_w)

      CALL cdf_close(ncid)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Interpolate fourier quantites from the restart file.
!>
!>  Restart files can change the number of surfaces and modes.
!>
!>  @param[in]  aold     Value from the restart file.
!>  @param[out] anew     New interpolated value.
!>  @param[in]  ns_old   Radial grid size in the restart file.
!>  @param[in]  ns_new   New radial grid size.
!>  @param[in]  mpol_old Number of poloidal modes in the restart file.
!>  @param[in]  mpol_new New number of poloidal modes.
!>  @param[in]  ntor_old Number of totoidal modes in the restart file.
!>  @param[in]  ntor_new New number of totoidal modes.
!>  @param[in]  nfp_old  Number of field periods in the restart file.
!>  @param[in]  nfp_new  New number of field periods.
!-------------------------------------------------------------------------------
      SUBROUTINE interpit(aold,     anew,                                      &
     &                    ns_old,   ns_new,                                    &
     &                    mpol_old, mpol_new,                                  &
     &                    ntor_old, ntor_new,                                  &
     &                    nfp_old,  nfp_new, lhalf)
      USE stel_constants, ONLY: one, zero

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), INTENT(in)  :: aold(0:mpol_old,-ntor_old:ntor_old,ns_old)
      REAL (dp), INTENT(out) :: anew(0:mpol_new,-ntor_new:ntor_new,ns_new)
      INTEGER, INTENT(in)    :: ns_old
      INTEGER, INTENT(in)    :: ns_new
      INTEGER, INTENT(in)    :: mpol_old
      INTEGER, INTENT(in)    :: mpol_new
      INTEGER, INTENT(in)    :: ntor_old
      INTEGER, INTENT(in)    :: ntor_new
      INTEGER, INTENT(in)    :: nfp_old
      INTEGER, INTENT(in)    :: nfp_new
      LOGICAL, INTENT(in)    :: lhalf

!  local variables
      INTEGER                :: real_ntor_min
      INTEGER                :: i_n
      INTEGER                :: i_m
      REAL (dp)              :: ds
      INTEGER                :: parity
      INTEGER                :: mpol_min
      INTEGER                :: ntor_min

!  Start of executable code
      mpol_min = MIN(mpol_new, mpol_old)
      real_ntor_min = MIN(ntor_new*nfp_new, ntor_old*nfp_old)

!  If the Fourier dimensions of the both arrays are the same just copy them.
!  Otherwise copy only the matcing modes.
      IF (ns_old   .eq. ns_new   .and.                                         &
          mpol_old .eq. mpol_new .and.                                         &
          ntor_old .eq. ntor_new .and.                                         &
          nfp_old  .eq. nfp_new) THEN
         anew = aold
      ELSE
!  Only some of the toroidal modes will match if nfp_new and nfp_old are
!  different.
         DO i_n = -real_ntor_min, real_ntor_min
            IF (i_n/nfp_new .eq. i_n/nfp_old) THEN
               DO i_m = 0, mpol_min
                  IF (MOD(i_m, 2) .EQ. 0) THEN
                     parity = -1
                  ELSE
                     parity = 1
                  END IF

                  CALL interpit_1d(aold(i_m, i_n, :), anew(i_m, i_n, :),       &
                                   ns_old, ns_new, lhalf, parity)
               END DO
            ELSE
               anew(:, i_n, :) = 0
            END IF
         END DO
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Interpolate radial quantites from the restart file.
!>
!>  Restart files can change the number of surfaces.
!>
!>  @param[in]  aold     Value from the restart file.
!>  @param[out] anew     New interpolated value.
!>  @param[in]  ns_old   Radial grid size in the restart file.
!>  @param[in]  ns_new   New radial grid size.
!-------------------------------------------------------------------------------
      SUBROUTINE interpit_1d(aold, anew, ns_old, ns_new, lhalf, parity)
      USE stel_constants, ONLY: one, zero

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), INTENT(in)                :: aold(ns_old)
      REAL (dp), INTENT(out)               :: anew(ns_new)
      INTEGER, INTENT(in)                  :: ns_old
      INTEGER, INTENT(in)                  :: ns_new
      LOGICAL, INTENT(in)                  :: lhalf
      INTEGER, INTENT(in)                  :: parity

!  local variables
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_y
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_x
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_new
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp
      REAL (dp)                            :: ds
      INTEGER                              :: i_s

!  local parameters
      REAL (dp), PARAMETER                 :: lower_b = -1.e30_dp
      REAL (dp), PARAMETER                 :: upper_b = -1.e30_dp

!  Start of executable code
      anew = 0
      IF (ns_old .EQ. ns_new) THEN
         anew = aold
         RETURN
      END IF

!  Need to radially interpolate. Extend the arrays from -ns to ns to get the
!  correct behavior at s=0.
      ds = 1.0/(ns_old - 1)
      IF (lhalf) THEN
         ALLOCATE(temp_y(2*(ns_old - 1)))
         ALLOCATE(temp_x(2*(ns_old - 1)))
         ALLOCATE(temp(2*(ns_old - 1)))

         DO i_s = 2, ns_old
            temp_x(ns_old - i_s + 1) = -ds*(i_s - 1.5)
            temp_y(ns_old + i_s - 2) = ds*(i_s - 1.5)
            temp_y(ns_old - i_s + 1) = parity*aold(i_s)
            temp_y(ns_old + i_s - 2) = aold(i_s)
         END DO
      ELSE
         ALLOCATE(temp_y(2*ns_old - 1))
         ALLOCATE(temp_x(2*ns_old - 1))
         ALLOCATE(temp(2*ns_old - 1))

         DO i_s = 1, ns_old
            temp_x(ns_old - i_s + 1) = -ds*(i_s - 1)
            temp_x(ns_old + i_s - 1) = ds*(i_s - 1)
            temp_y(ns_old - i_s + 1) = parity*aold(i_s)
            temp_y(ns_old + i_s - 1) = aold(i_s)
         END DO
      END IF

      CALL spline(temp_x, temp_y, SIZE(temp_x), lower_b, upper_b, temp)

      ALLOCATE(temp_new(ns_new))
      ds = 1.0/(ns_new - 1)
      IF (lhalf) THEN
         DO i_s = 2, ns_new
            temp_new(i_s) = ds*(i_s - 1.5)
         END DO
      ELSE
         DO i_s = 1, ns_new
            temp_new(i_s) = ds*(i_s - 1)
         END DO
      END IF

      CALL splINT(temp_x, temp_y, temp, SIZE(temp_x), temp_new, anew)

      DEALLOCATE(temp_y)
      DEALLOCATE(temp_x)
      DEALLOCATE(temp)

      END SUBROUTINE

      END MODULE
