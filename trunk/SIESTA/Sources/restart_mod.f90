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

      IMPLICIT NONE

!*******************************************************************************
!  restart module parameters
!*******************************************************************************
!>  Version number.
      INTEGER, PARAMETER :: restart_version = 0

!  Flag will be placed in the last bits. This way the version number and flags 
!  can fit in the same memory.
!>  Bit position for the l_asym flag.
      INTEGER, PARAMETER           :: restart_l_asym = 31
!>  Bit position for the l_recon flag.
      INTEGER, PARAMETER           :: restart_l_recon = 30

!>  Dimension names.
      CHARACTER (LEN=*), DIMENSION(3), PARAMETER ::                            &
     &   restart_dims = (/'m-mode','n-mode','radius'/)

!>  Name for the restart file number of radial points.
      CHARACTER (len=*), PARAMETER :: vn_nsin = 'nrad'
!>  Name for the restart file number of poloidal modes.
      CHARACTER (len=*), PARAMETER :: vn_mpolin = 'mpol'
!>  Name for the restart file number of toroidal modes.
      CHARACTER (len=*), PARAMETER :: vn_ntorin = 'ntor'
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
!>  @param[in]    restart_file_name Name of the restart file.
!>  @param[in]    mpolin            Namelist number of polodal modes.
!>  @param[in]    ntorin            Namelist number of toroidal modes.
!>  @param[in]    nsin              Namelist number of radial grid points.
!>  @param[inout] wout_file         Name of the wout file.
!>  @returns Error status of the read.
!-------------------------------------------------------------------------------
      FUNCTION restart_read(restart_file_name, mpolin, ntorin, nsin, wout_file)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch
      USE shared_data, ONLY: l_asym, nprecon

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                     :: restart_read
      CHARACTER (len=*), INTENT(in)               :: restart_file_name
      INTEGER, INTENT(in)                         :: mpolin
      INTEGER, INTENT(in)                         :: ntorin
      INTEGER, INTENT(in)                         :: nsin
      CHARACTER (len=*), INTENT(inout)            :: wout_file
!  local variables
      INTEGER                                     :: flags
      INTEGER                                     :: ncid
      INTEGER                                     :: ns
      INTEGER                                     :: mpol
      INTEGER                                     :: ntor
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: tempmn
      INTEGER                                     :: nmin

!  local parameters
      INTEGER, PARAMETER                          :: m0=0
      INTEGER, PARAMETER                          :: m1=1

!  Start of executable code
      CALL cdf_open(ncid, TRIM(restart_file_name), 'r', restart_read)
      IF (restart_read .ne. 0) RETURN

      CALL cdf_read(ncid, vn_flags, flags)

!  The namelist inut file can change the size of the radial grid and the 
!  poloidal and toroidal modes.
      CALL cdf_read(ncid, vn_nsin, ns)
      CALL cdf_read(ncid, vn_mpolin, mpol)
      CALL cdf_read(ncid, vn_ntorin, ntor)

      nmin = MIN(ntor, ntorin)

      CALL cdf_read(ncid, vn_wout, wout_file)

!  The namelist input file may turn the asymmetric terms on and off.
      ALLOCATE(tempmn(0:mpol,-ntor:ntor, ns))

      CALL cdf_read(ncid, vn_jbsupss, tempmn)
      jbsupsmnsh(:,:,1) = 0
      CALL interpit(tempmn, jbsupsmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupsmnsh(m1,-nmin:nmin,1) = jbsupsmnsh(m1,-nmin:nmin,2)

      CALL cdf_read(ncid, vn_jbsupuc, tempmn)
      CALL interpit(tempmn, jbsupumnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupumnch(m1,-nmin:nmin,1) = jbsupumnch(m1,-nmin:nmin,2)

      jbsupvmnch(:,:,1) = 0
      CALL cdf_read(ncid, vn_jbsupvc, tempmn)
      CALL interpit(tempmn, jbsupvmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupvmnch(m0,-nmin:nmin,1) = jbsupvmnch(m0,-nmin:nmin,2)

      jpmnch(:,:,1) = 0
      CALL cdf_read(ncid, vn_jpresc,  tempmn)
      CALL interpit(tempmn, jpmnch,     ns, nsin, mpol, mpolin, ntor, ntorin)
      jpmnch(m0,-nmin:nmin,1) = jpmnch(m0,-nmin:nmin,2)

      IF (BTEST(flags, restart_l_asym) .and. l_asym) THEN
         jbsupsmnch(:,:,1) = 0
         CALL cdf_read(ncid, vn_jbsupsc, tempmn)
         CALL interpit(tempmn, jbsupsmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupsmnch(m1,-nmin:nmin,1) = jbsupsmnch(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupus, tempmn)
         CALL interpit(tempmn, jbsupumnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupumnsh(m1,-nmin:nmin,1) = jbsupumnsh(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupvs, tempmn)
         CALL interpit(tempmn, jbsupvmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupvmnsh(m0,-nmin:nmin,1) = jbsupvmnsh(m0,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jpress,  tempmn)
         CALL interpit(tempmn, jpmnsh,     ns, nsin, mpol, mpolin, ntor, ntorin)
         jpmnsh(m0,-nmin:nmin,1) = jpmnsh(m0,-nmin:nmin,2)
      ELSE
         jbsupsmnch = 0
         jbsupumnsh = 0
         jbsupvmnsh = 0
         jpmnsh = 0
      END IF

      CALL cdf_close(ncid)

      DEALLOCATE(tempmn)

      WRITE (33, 400) mpol, ntor, ns
      nprecon = 2   !USED TO SKIP DIAGONAL PRECONDITIONER IN EVOLVE

 400  FORMAT(/,' RESTARTED FROM RUN PARAMETERS M: ',i3,' N: ',i3,' NS: ', i3,/)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write the restart file.
!>
!>  Writes the restart information.
!>
!>  @param[in] restart_file_name Name of the restart file.
!>  @param[in] wout_file         Name of the wout file.
!>  @returns Error status of the write.
!-------------------------------------------------------------------------------
      FUNCTION restart_write(restart_file_name, wout_file)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch,                                &
                            b_factor,   p_factor, jacobh
      USE fourier
      USE shared_data, ONLY: l_asym, l_recon
      USE island_params, ONLY: ohs => ohs_i

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                     :: restart_write
      CHARACTER (len=*), INTENT(in)               :: restart_file_name
      CHARACTER (len=*), INTENT(in)               :: wout_file

!  local variables
      INTEGER                                     :: flags
      INTEGER                                     :: ncid
      INTEGER                                     :: nblock
      INTEGER                                     :: m, n
      REAL (rprec)                                :: r0
      REAL (rprec)                                :: tempscalar
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsupsijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsupvijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubsijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubuijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubvijh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmnh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubumnh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubvmnh
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubumnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsubvmnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupsmnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupumnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: jksupvmnf
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: dbuds
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: dbvds
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: tempmn
      REAL (rprec), DIMENSION(:,:,:), ALLOCATABLE :: pijh

!  local parameters
      INTEGER, PARAMETER                          :: m0=0
      INTEGER, PARAMETER                          :: m1=1
      INTEGER, PARAMETER                          :: m2=2

!  Start of executable code
      CALL cdf_open(ncid, TRIM(restart_file_name), 'w', restart_write)
      IF (restart_write .NE. 0) RETURN

      flags = restart_version
      IF (l_asym) THEN
         flags = IBSET(flags, restart_l_asym)
      END IF
      IF (l_recon) THEN
         flags = IBSET(flags, restart_l_recon)
      END IF

      CALL cdf_define(ncid, vn_flags, flags)
      CALL cdf_define(ncid, vn_nsin, ns)
      CALL cdf_define(ncid, vn_mpolin, mpol)
      CALL cdf_define(ncid, vn_ntorin, ntor)
      CALL cdf_define(ncid, vn_wout, wout_file)

      CALL cdf_define(ncid, vn_jbsupss, jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupuc, jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupvc, jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jpresc,  jpmnch,     dimname=restart_dims)

      IF (l_asym) THEN
         CALL cdf_define(ncid, vn_jbsupsc, jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupus, jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupvs, jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jpress,  jpmnsh,     dimname=restart_dims)
      END IF

      IF (l_recon) THEN
!  Using variable with jacobian only for the dimension sizes and types. Be for
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
         CALL cdf_define(ncid, vn_p_max,     b_factor)
         CALL cdf_define(ncid, vn_p_min,     p_factor)
         IF (l_asym) THEN
            CALL cdf_define(ncid, vn_bsupsmnc,  jbsupsmnch, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsupumns,  jbsupumnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsupvmns,  jbsupvmnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubsmnc,  jbsupsmnch, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubumns,  jbsupumnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_bsubvmns,  jbsupvmnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupsmnc, jbsupsmnch, dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupumns, jbsupumnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_jksupvmns, jbsupvmnsh, dimname=restart_dims)
            CALL cdf_define(ncid, vn_pmns,      jpmnsh,     dimname=restart_dims)
         END IF
      END IF

      CALL cdf_write(ncid, vn_flags, flags)
      CALL cdf_write(ncid, vn_nsin, ns)
      CALL cdf_write(ncid, vn_mpolin, mpol)
      CALL cdf_write(ncid, vn_ntorin, ntor)
      CALL cdf_write(ncid, vn_wout, wout_file)

      CALL cdf_write(ncid, vn_jbsupss, jbsupsmnsh)
      CALL cdf_write(ncid, vn_jbsupuc, jbsupumnch)
      CALL cdf_write(ncid, vn_jbsupvc, jbsupvmnch)
      CALL cdf_write(ncid, vn_jpresc,  jpmnch)

      IF (l_asym) THEN
         CALL cdf_write(ncid, vn_jbsupsc, jbsupsmnch)
         CALL cdf_write(ncid, vn_jbsupus, jbsupumnsh)
         CALL cdf_write(ncid, vn_jbsupvs, jbsupvmnsh)
         CALL cdf_write(ncid, vn_jpress,  jpmnsh)
      END IF

      IF (l_recon) THEN
         CALL cdf_write(ncid, vn_p_factor, p_factor)
         CALL cdf_write(ncid, vn_b_factor, b_factor)

         ALLOCATE(bsupsijh(ntheta,nzeta,ns))
         ALLOCATE(bsupuijh(ntheta,nzeta,ns))
         ALLOCATE(bsupvijh(ntheta,nzeta,ns))
         ALLOCATE(pijh(ntheta,nzeta,ns))

         bsupsijh = 0
         bsupuijh = 0
         bsupvijh = 0
         pijh = 0

         CALL toijsp(jbsupsmnsh, bsupsijh, f_none, f_sin, f_half)
         CALL toijsp(jbsupumnch, bsupuijh, f_none, f_cos, f_half)
         CALL toijsp(jbsupvmnch, bsupvijh, f_none, f_cos, f_half)
         CALL toijsp(jpmnch, pijh, f_none, f_cos, f_half)

         IF (l_asym) THEN
            CALL toijsp(jbsupsmnch, bsupsijh, f_sum, f_cos, f_half)
            CALL toijsp(jbsupumnsh, bsupuijh, f_sum, f_sin, f_half)
            CALL toijsp(jbsupvmnsh, bsupvijh, f_sum, f_sin, f_half)
            CALL toijsp(jpmnsh, pijh, f_sum, f_sin, f_half)
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

         ALLOCATE(tempmn(0:mpol,-ntor:ntor,ns))
         tempmn = 0

         CALL tomnsp(bsupsijh, tempmn, f_sin, f_half)
         CALL cdf_write(ncid, vn_bsupsmns, tempmn)
         CALL tomnsp(bsupuijh, tempmn, f_cos, f_half)
         CALL cdf_write(ncid, vn_bsupumnc, tempmn)
         CALL tomnsp(bsupvijh, tempmn, f_cos, f_half)
         CALL cdf_write(ncid, vn_bsupvmnc, tempmn)
         CALL tomnsp(pijh, tempmn, f_cos, f_half)
         CALL cdf_write(ncid, vn_pmnc, tempmn)

         ALLOCATE(bsubsijh(ntheta,nzeta,ns))
         ALLOCATE(bsubuijh(ntheta,nzeta,ns))
         ALLOCATE(bsubvijh(ntheta,nzeta,ns))

         CALL tolowerh(bsupsijh, bsupuijh, bsupvijh, bsubsijh, bsubuijh, bsubvijh)

!  Need these to compute the currents later.
         ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnh(0:mpol,-ntor:ntor,ns))

         CALL tomnsp(bsubsijh, bsubsmnh, f_sin, f_half)
         CALL cdf_write(ncid, vn_bsubsmns, bsubsmnh)
         CALL tomnsp(bsubuijh, bsubumnh, f_cos, f_half)
         CALL cdf_write(ncid, vn_bsubumnc, bsubumnh)
         CALL tomnsp(bsubvijh, bsubvmnh, f_cos, f_half)
         CALL cdf_write(ncid, vn_bsubvmnc, bsubvmnh)

!  Need full mesh values to compute currents on the full mesh.
         ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnf(0:mpol,-ntor:ntor,ns))

         nblock = (2*ntor + 1)*(mpol + 1)
         CALL to_full_mesh(bsubsmnh, bsubsmnf, nblock)
         CALL to_full_mesh(bsubumnh, bsubumnf, nblock)
         CALL to_full_mesh(bsubvmnh, bsubvmnf, nblock)

!         bsubsmnf(:,:,1) = 0
!         bsubumnf(:,:,1) = 0
!         bsubsmnf(m1,:,1) = bsubsmnh(m1,:,1)
!         bsubumnf(m1,:,1) = bsubumnh(m1,:,1)
!         bsubvmnf(:,:,1) = 0
!         bsubvmnf(m0:m1,:,1) = bsubvmnh(m0:m1,:,2)

         ALLOCATE(jksupsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupvmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbuds(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbvds(0:mpol,-ntor:ntor,ns))

!  Half grid radial derivatives become full grid.
         CALL GradientFull(dbuds, bsubumnh, nblock)
         CALL GradientFull(dbvds, bsubvmnh, nblock)

!         r0 = hs_i/2
!         dbuds(m0:m1,:,1) = 2*ohs*bsubumnh(m0:m1,:,2)
!         bsubvmnf(m1,:,1) = bsubvmnh(m1,:,2)
!         dbvds(m1,:,1)    = 2*ohs*bsubvmnh(m1,:,2)                      !--> 0
!         bsubumnf(m1,:,1) = r0*bsubsmnf(m1,:,1)

         DO n = -ntor, ntor
            DO m = 0, mpol
               jksupsmnf(m,n,:) = -m*    bsubvmnf(m,n,:) + n*nfp*bsubumnf(m,n,:)
               jksupumnf(m,n,:) =  n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
               jksupvmnf(m,n,:) =  dbuds(m,n,:)          - m*    bsubsmnf(m,n,:)
            END DO
         END DO

         CALL cdf_write(ncid, vn_jksupsmns, jksupsmnf)
         CALL cdf_write(ncid, vn_jksupumnc, jksupumnf)
         CALL cdf_write(ncid, vn_jksupvmnc, jksupvmnf)

         IF (l_asym) THEN
            CALL tomnsp(bsupsijh, tempmn, f_cos, f_half)
            CALL cdf_write(ncid, vn_bsupsmnc, tempmn)
            CALL tomnsp(bsupuijh, tempmn, f_sin, f_half)
            CALL cdf_write(ncid, vn_bsupumns, tempmn)
            CALL tomnsp(bsupvijh, tempmn, f_sin, f_half)
            CALL cdf_write(ncid, vn_bsupvmns, tempmn)
            CALL tomnsp(pijh, tempmn, f_sin, f_half)
            CALL cdf_write(ncid, vn_pmns, tempmn)

            CALL tomnsp(bsubsijh, bsubsmnh, f_cos, f_half)
            CALL cdf_write(ncid, vn_bsubsmnc, bsubsmnh)
            CALL tomnsp(bsubuijh, bsubumnh, f_sin, f_half)
            CALL cdf_write(ncid, vn_bsubumns, bsubumnh)
            CALL tomnsp(bsubvijh, bsubvmnh, f_sin, f_half)
            CALL cdf_write(ncid, vn_bsubvmns, bsubvmnh)

            CALL to_full_mesh(bsubsmnh, bsubsmnf, nblock)
            CALL to_full_mesh(bsubumnh, bsubumnf, nblock)
            CALL to_full_mesh(bsubvmnh, bsubvmnf, nblock)

!            bsubsmnf(:,:,1) = 0
!            bsubumnf(:,:,1) = 0
!            bsubsmnf(m1,:,1) = bsubsmnh(m1,:,1)
!            bsubumnf(m1,:,1) = bsubumnh(m1,:,1)
!            bsubvmnf(:,:,1) = 0
!            bsubvmnf(m0:m1,:,1) = bsubvmnh(m0:m1,:,2)

            CALL GradientFull(dbuds, bsubumnh, nblock)
            CALL GradientFull(dbvds, bsubvmnh, nblock)

!            dbuds(m0:m1,:,1) = 2*ohs*bsubumnh(m0:m1,:,2)
!            bsubvmnf(m1,:,1) = bsubvmnh(m1,:,2)
!            dbvds(m1,:,1)    = 2*ohs*bsubvmnh(m1,:,2)                      !--> 0
!            bsubumnf(m1,:,1) = r0*bsubsmnf(m1,:,1)

            DO n = -ntor, ntor
               DO m = 0, mpol
                  jksupsmnf(m,n,:) =  m*    bsubvmnf(m,n,:) - n*nfp*bsubumnf(m,n,:)
                  jksupumnf(m,n,:) = -n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
                  jksupvmnf(m,n,:) =  dbuds(m,n,:)          + m*    bsubsmnf(m,n,:)
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

         DEALLOCATE(tempmn)
      END IF

      CALL cdf_close(ncid)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Reads the restart file.
!>
!>  Reads the restart information and initalizes SIESTA quantities.
!>
!>  @param[in]  aold     Value from the restart file.
!>  @param[out] anew     New interpolated value.
!>  @param[in]  ns_old   Radial grid size in the restart file.
!>  @param[in]  ns_new   New radial grid size.
!>  @param[in]  mpol_old Number of poloidal modes in the restart file.
!>  @param[in]  mpol_new New number of poloidal modes.
!>  @param[in]  ntor_old Number of totoidal modes in the restart file.
!>  @param[in]  ntor_new New number of totoidal modes.
!-------------------------------------------------------------------------------
      SUBROUTINE interpit(aold, anew, ns_old, ns_new,                          &
     &                    mpol_old, mpol_new, ntor_old, ntor_new)
      USE stel_constants, ONLY: one, zero
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(in)  :: aold(0:mpol_old,-ntor_old:ntor_old,ns_old)
      REAL(rprec), INTENT(out) :: anew(0:mpol_new,-ntor_new:ntor_new,ns_new)
      INTEGER, INTENT(in)      :: ns_old
      INTEGER, INTENT(in)      :: ns_new
      INTEGER, INTENT(in)      :: mpol_old
      INTEGER, INTENT(in)      :: mpol_new
      INTEGER, INTENT(in)      :: ntor_old
      INTEGER, INTENT(in)      :: ntor_new
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER   :: p5 = 0.5_dp
      REAL(rprec), ALLOCATABLE :: temp1(:,:,:)
      INTEGER                  :: mpol_min, ntor_min
      INTEGER                  :: js_old, js_new, istat
      REAL(rprec)              :: hs_old, hs_new, snew, sold, x1
!-----------------------------------------------

      mpol_min = MIN(mpol_old, mpol_new)
      ntor_min = MIN(ntor_old, ntor_new)

      anew = 0

      IF (ns_old .EQ. ns_new) THEN
         anew(0:mpol_min,-ntor_min:ntor_min,:) = aold(0:mpol_min,-ntor_min:ntor_min,:)
         RETURN
      END IF

      ALLOCATE (temp1(0:mpol_min,-ntor_min:ntor_min,ns_old), stat=istat)
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
            anew(:,:,js_new) = (one-x1)*temp1(:,:,js_old)                      &
                             + x1*temp1(:,:,js_old+1)
         END IF
      END DO
!
!     AVERAGE BACK TO HALF-MESH
!
      DEALLOCATE (temp1)
      ALLOCATE (temp1(0:mpol_min,-ntor_min:ntor_min,ns_new), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error #2 in INTERPIT'
      temp1 = anew
      anew(:,:,1) = 0
      anew(:,:,2:ns_new) = p5*(temp1(:,:,2:ns_new) + temp1(:,:,1:ns_new-1))
      DEALLOCATE(temp1)

      END SUBROUTINE interpit

      END MODULE
