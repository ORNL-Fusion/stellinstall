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
      USE v3_utilities, ONLY: assert
      USE ezcdf
      USE stel_kinds
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE metrics, ONLY: tolowerh
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: lasym, lrecon, unit_out

      IMPLICIT NONE

!*******************************************************************************
!  restart module parameters
!*******************************************************************************
!>  Version number.
      INTEGER, PARAMETER :: restart_version = 0

!  Flag will be placed in the last bits. This way the version number and flags 
!  can fit in the same memory.
!>  Bit position for the lasym flag.
      INTEGER, PARAMETER           :: restart_lasym = 31
!>  Bit position for the lrecon flag.
      INTEGER, PARAMETER           :: restart_lrecon = 30

!>  Dimension names.
      CHARACTER (LEN=*), DIMENSION(3), PARAMETER ::                     &
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

!>  Name for the restart file rmnc.
      CHARACTER (len=*), PARAMETER :: vn_rmnc = 'rmnc(m,n,r)'
!>  Name for the restart file rmns.
      CHARACTER (len=*), PARAMETER :: vn_rmns = 'rmns(m,n,r)'
!>  Name for the restart file zmnc.
      CHARACTER (len=*), PARAMETER :: vn_zmnc = 'zmnc(m,n,r)'
!>  Name for the restart file zmns.
      CHARACTER (len=*), PARAMETER :: vn_zmns = 'zmns(m,n,r)'

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
      FUNCTION restart_read(restart_file_name, wout_file, mpolin, ntorin, nsin)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                     &
                            jbsupumnsh, jbsupumnch,                     &
                            jbsupvmnsh, jbsupvmnch,                     &
                            jpmnsh,     jpmnch

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
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE    :: tempmn_r
      INTEGER                                     :: nmin

!  local parameters
      INTEGER, PARAMETER                          :: m0=0
      INTEGER, PARAMETER                          :: m1=1
      CHARACTER (LEN=256)                         :: filename
      LOGICAL                                     :: lex

!  Start of executable code

      restart_read = 1

#if defined(NETCDF)
      filename = TRIM(restart_file_name) // ".nc"
      INQUIRE(file=TRIM(filename), exist=lex)
      IF (lex) THEN
         CALL cdf_open(ncid, TRIM(filename), 'r', restart_read)
      END IF          
#endif            
      IF (restart_read .NE. 0) THEN
          filename = TRIM(restart_file_name) // ".txt"
          restart_read = restart_read_txt(filename, mpolin, ntorin, nsin)
      END IF

      IF (restart_read .NE. 0) THEN
         IF (iam .EQ. 0) THEN
            filename = 'UNABLE TO READ SIESTA RESTART FILE: ' //        &
                        TRIM(restart_file_name)
            PRINT 50, TRIM(filename)
         END IF
         RETURN
      ELSE
         IF (iam .EQ. 0) THEN
            WRITE (unit_out, 400) mpol, ntor, ns
            filename = 'RESTART FILE WAS SUCCESSFULLY READ: ' //        &
                        TRIM(restart_file_name)
            PRINT 50, TRIM(filename)
            CALL FLUSH(unit_out)
         END IF
      END IF
 50   FORMAT(' ************',/,1x,a,/,' ************')
400   FORMAT(/,' RESTARTED FROM RUN PARAMETERS M: ',i3,' N: ',i3,' NS: ', i3,/)

#if defined(NETCDF)
      CALL cdf_read(ncid, vn_flags, flags)

!  The namelist inut file can change the size of the radial grid and the 
!  poloidal and toroidal modes.
      CALL cdf_read(ncid, vn_nsin, ns)
      CALL cdf_read(ncid, vn_mpolin, mpol)
      CALL cdf_read(ncid, vn_ntorin, ntor)

      nmin = MIN(ntor, ntorin)

      CALL cdf_read(ncid, vn_wout, wout_file)

!  The namelist input file may turn the asymmetric terms on and off.
      ALLOCATE(tempmn_r(0:mpol,-ntor:ntor, ns))

      CALL cdf_read(ncid, vn_jbsupss, tempmn_r)
      jbsupsmnsh(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupsmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupsmnsh(m1,-nmin:nmin,1) = jbsupsmnsh(m1,-nmin:nmin,2)
      
      CALL cdf_read(ncid, vn_jbsupuc, tempmn_r)
      CALL interpit(tempmn_r, jbsupumnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupumnch(m1,-nmin:nmin,1) = jbsupumnch(m1,-nmin:nmin,2)

      jbsupvmnch(:,:,1) = 0
      CALL cdf_read(ncid, vn_jbsupvc, tempmn_r)
      CALL interpit(tempmn_r, jbsupvmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupvmnch(m0,-nmin:nmin,1) = jbsupvmnch(m0,-nmin:nmin,2)

      jpmnch(:,:,1) = 0
      CALL cdf_read(ncid, vn_jpresc,  tempmn_r)
      CALL interpit(tempmn_r, jpmnch,     ns, nsin, mpol, mpolin, ntor, ntorin)
      jpmnch(m0,-nmin:nmin,1) = jpmnch(m0,-nmin:nmin,2)


      IF (BTEST(flags, restart_lasym) .and. lasym) THEN
         jbsupsmnch(:,:,1) = 0
         CALL cdf_read(ncid, vn_jbsupsc, tempmn_r)
         CALL interpit(tempmn_r, jbsupsmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupsmnch(m1,-nmin:nmin,1) = jbsupsmnch(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupus, tempmn_r)
         CALL interpit(tempmn_r, jbsupumnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupumnsh(m1,-nmin:nmin,1) = jbsupumnsh(m1,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jbsupvs, tempmn_r)
         CALL interpit(tempmn_r, jbsupvmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
         jbsupvmnsh(m0,-nmin:nmin,1) = jbsupvmnsh(m0,-nmin:nmin,2)

         CALL cdf_read(ncid, vn_jpress,  tempmn_r)
         CALL interpit(tempmn_r, jpmnsh,     ns, nsin, mpol, mpolin, ntor, ntorin)
         jpmnsh(m0,-nmin:nmin,1) = jpmnsh(m0,-nmin:nmin,2)
      ELSE IF (lasym) THEN
         jbsupsmnch = 0
         jbsupumnsh = 0
         jbsupvmnsh = 0
         jpmnsh = 0
      END IF

      CALL cdf_close(ncid)

      DEALLOCATE(tempmn_r)
#endif
      END FUNCTION restart_read


      FUNCTION restart_read_txt(filename, mpolin, ntorin, nsin)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                     &
                            jbsupumnsh, jbsupumnch,                     &
                            jbsupvmnsh, jbsupvmnch,                     &
                            jpmnsh,     jpmnch
      INTEGER, PARAMETER :: m0=0, m1=1
      INTEGER  :: restart_read_txt, mpolin, ntorin, nsin
      INTEGER  :: nmin, mpol, ntor, ns, istat, m, n, js
      CHARACTER (LEN=*), INTENT(IN)  :: filename
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE    ::                    &
                 bsupss, bsupuc, bsupvc, presc,                         &
                 bsupsc, bsupus, bsupvs, press
    
      OPEN (unit=15, file=TRIM(filename), status='old',                 &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) THEN 
         restart_read_txt = 1
         RETURN
      END IF

      restart_read_txt = 0
!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      READ (15,'(L1)') lasym
      READ (15,100) ns, mpol, ntor
      nmin = MIN(ntor, ntorin)

      ALLOCATE (bsupss(0:mpol,-ntor:ntor,ns),                           &
                bsupuc(0:mpol,-ntor:ntor,ns),                           &
                bsupvc(0:mpol,-ntor:ntor,ns),                           & 
                presc(0:mpol,-ntor:ntor,ns), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in restart_read_txt')

      DO js=1,ns
      DO m=0,mpol
      DO n=-ntor,ntor
      READ (15,300) bsupss(m,n,js), bsupuc(m,n,js),                     &
                    bsupvc(m,n,js), presc(m,n,js)
      ENDDO
      ENDDO
      ENDDO

      IF (lasym) THEN
      ALLOCATE (bsupsc(0:mpol,-ntor:ntor,ns),                           &
                bsupus(0:mpol,-ntor:ntor,ns),                           &
                bsupvs(0:mpol,-ntor:ntor,ns),                           & 
                press(0:mpol,-ntor:ntor,ns), stat=istat)
      DO js=1,ns
      DO m=0,mpol
      DO n=-ntor,ntor
      READ (15,300) bsupsc(m,n,js), bsupus(m,n,js),                     &
                    bsupvs(m,n,js), press(m,n,js)
      ENDDO
      ENDDO
      ENDDO
      END IF

      CLOSE (unit=15)

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

 100  FORMAT(3i8)
 300  FORMAT(1p4e24.16)
      
      END FUNCTION restart_read_txt

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
      USE island_params, ONLY: ohs => ohs_i
      USE shared_data, ONLY: r1_i, z1_i

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: restart_write
      CHARACTER (len=*), INTENT(in)            :: restart_file_name
      CHARACTER (len=*), INTENT(in)            :: wout_file

!  local variables
      INTEGER                                  :: flags
      INTEGER                                  :: ncid
      INTEGER                                  :: nblock
      INTEGER                                  :: s, m, n

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

!  local parameters
      INTEGER, PARAMETER                       :: m0=0
      INTEGER, PARAMETER                       :: m1=1
      INTEGER, PARAMETER                       :: m2=2
      CHARACTER (LEN=256)                      :: filename

!  Start of executable code
      restart_write = 1
#if defined(NETCDF)
      filename = TRIM(restart_file_name) // ".nc"
      CALL cdf_open(ncid, TRIM(filename), 'w', restart_write)
#endif            
      IF (restart_write .NE. 0) THEN
         filename = TRIM(restart_file_name) // ".txt"
         restart_write = restart_write_txt(filename)
         RETURN
      END IF

      flags = restart_version
      IF (lasym) THEN
         flags = IBSET(flags, restart_lasym)
      END IF
      IF (lrecon) THEN
         flags = IBSET(flags, restart_lrecon)
      END IF

      ALLOCATE(tempmn_w(0:mpol,-ntor:ntor,ns))

#if defined(NETCDF)
      CALL cdf_define(ncid, vn_flags, flags)
      CALL cdf_define(ncid, vn_nsin, ns)
      CALL cdf_define(ncid, vn_mpolin, mpol)
      CALL cdf_define(ncid, vn_ntorin, ntor)
      CALL cdf_define(ncid, vn_wout, wout_file)

      CALL cdf_define(ncid, vn_jbsupss, jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupuc, jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupvc, jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jpresc,  jpmnch,     dimname=restart_dims)

      CALL cdf_define(ncid, vn_rmnc, tempmn_w, dimname=restart_dims)
      CALL cdf_define(ncid, vn_zmns, tempmn_w, dimname=restart_dims)

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
         CALL cdf_define(ncid, vn_p_max,     b_factor)
         CALL cdf_define(ncid, vn_p_min,     p_factor)
         IF (lasym) THEN
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
      END IF RECON0

      CALL cdf_write(ncid, vn_flags, flags)
      CALL cdf_write(ncid, vn_nsin, ns)
      CALL cdf_write(ncid, vn_mpolin, mpol)
      CALL cdf_write(ncid, vn_ntorin, ntor)
      CALL cdf_write(ncid, vn_wout, wout_file)

      CALL cdf_write(ncid, vn_jbsupss, jbsupsmnsh)
      CALL cdf_write(ncid, vn_jbsupuc, jbsupumnch)
      CALL cdf_write(ncid, vn_jbsupvc, jbsupvmnch)
      CALL cdf_write(ncid, vn_jpresc,  jpmnch)

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

      RECON1: IF (lrecon) THEN

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

         CALL toijsp(jbsupsmnsh, bsupsijh, f_none, f_sin)
         CALL toijsp(jbsupumnch, bsupuijh, f_none, f_cos)
         CALL toijsp(jbsupvmnch, bsupvijh, f_none, f_cos)
         CALL toijsp(jpmnch, pijh, f_none, f_cos)

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

         CALL tomnsp(bsupsijh, tempmn_w, f_sin)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsupsmns, tempmn_w)
         CALL tomnsp(bsupuijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsupumnc, tempmn_w)
         CALL tomnsp(bsupvijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsupvmnc, tempmn_w)
         CALL tomnsp(pijh, tempmn_w, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_pmnc, tempmn_w)

         ALLOCATE(bsubsijh(ntheta,nzeta,ns))
         ALLOCATE(bsubuijh(ntheta,nzeta,ns))
         ALLOCATE(bsubvijh(ntheta,nzeta,ns))

         CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                    &
                       bsubsijh, bsubuijh, bsubvijh, 1, ns)

!  Need these to compute the currents later.
         ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnh(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnh(0:mpol,-ntor:ntor,ns))

         CALL tomnsp(bsubsijh, bsubsmnh, f_sin)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubsmnh(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsubsmns, tempmn_w)
         CALL tomnsp(bsubuijh, bsubumnh, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubumnh(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsubumnc, tempmn_w)
         CALL tomnsp(bsubvijh, bsubvmnh, f_cos)
         DO s = 1, ns
            tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubvmnh(:,:,s)
         END DO
         CALL cdf_write(ncid, vn_bsubvmnc, tempmn_w)

!  Need full mesh values to compute currents on the full mesh.
         ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(bsubvmnf(0:mpol,-ntor:ntor,ns))

         CALL to_full_mesh(bsubsmnh, bsubsmnf)
         CALL to_full_mesh(bsubumnh, bsubumnf)
         CALL to_full_mesh(bsubvmnh, bsubvmnf)

         ALLOCATE(jksupsmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupumnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jksupvmnf(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbuds(0:mpol,-ntor:ntor,ns))
         ALLOCATE(dbvds(0:mpol,-ntor:ntor,ns))

!  Half grid radial derivatives become full grid.
         CALL GradientFull(dbuds, bsubumnh)
         CALL GradientFull(dbvds, bsubvmnh)

         DO n = -ntor, ntor
            DO m = 0, mpol
               jksupsmnf(m,n,:) = -m*    bsubvmnf(m,n,:) + n*nfp*bsubumnf(m,n,:)
               jksupumnf(m,n,:) =  n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
               jksupvmnf(m,n,:) =  dbuds(m,n,:)          - m*    bsubsmnf(m,n,:)
               jksupsmnf(m,n,:) =  orthonorm(m,n)*jksupsmnf(m,n,:)
               jksupumnf(m,n,:) =  orthonorm(m,n)*jksupumnf(m,n,:)
               jksupvmnf(m,n,:) =  orthonorm(m,n)*jksupvmnf(m,n,:)
            END DO
         END DO

         CALL cdf_write(ncid, vn_jksupsmns, jksupsmnf)
         CALL cdf_write(ncid, vn_jksupumnc, jksupumnf)
         CALL cdf_write(ncid, vn_jksupvmnc, jksupvmnf)

         IF (lasym) THEN
            CALL tomnsp(bsupsijh, tempmn_w, f_cos)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsupsmnc, tempmn_w)
            CALL tomnsp(bsupuijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsupumns, tempmn_w)
            CALL tomnsp(bsupvijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsupvmns, tempmn_w)
            CALL tomnsp(pijh, tempmn_w, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*tempmn_w(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_pmns, tempmn_w)

            CALL tomnsp(bsubsijh, bsubsmnh, f_cos)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubsmnh(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsubsmnc, tempmn_w)
            CALL tomnsp(bsubuijh, bsubumnh, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubumnh(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsubumns, tempmn_w)
            CALL tomnsp(bsubvijh, bsubvmnh, f_sin)
            DO s = 1, ns
               tempmn_w(:,:,s) = orthonorm(0:mpol,-ntor:ntor)*bsubvmnh(:,:,s)
            END DO
            CALL cdf_write(ncid, vn_bsubvmns, tempmn_w)

            CALL to_full_mesh(bsubsmnh, bsubsmnf)
            CALL to_full_mesh(bsubumnh, bsubumnf)
            CALL to_full_mesh(bsubvmnh, bsubvmnf)

            CALL GradientFull(dbuds, bsubumnh)
            CALL GradientFull(dbvds, bsubvmnh)

            DO n = -ntor, ntor
               DO m = 0, mpol
                  jksupsmnf(m,n,:) =  m*    bsubvmnf(m,n,:) - n*nfp*bsubumnf(m,n,:)
                  jksupumnf(m,n,:) = -n*nfp*bsubsmnf(m,n,:) - dbvds(m,n,:)
                  jksupvmnf(m,n,:) =  dbuds(m,n,:)          + m*    bsubsmnf(m,n,:)
                  jksupsmnf(m,n,:) =  orthonorm(m,n)*jksupsmnf(m,n,:)
                  jksupumnf(m,n,:) =  orthonorm(m,n)*jksupumnf(m,n,:)
                  jksupvmnf(m,n,:) =  orthonorm(m,n)*jksupvmnf(m,n,:)
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
#endif
      END FUNCTION restart_write

      FUNCTION restart_write_txt(filename)
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,  &
                            jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh,  &
                            ns, mpol, ntor
      CHARACTER (LEN=*)     :: filename
      INTEGER               :: restart_write_txt, istat, js, m, n
      
      OPEN (unit=15, file=TRIM(filename), status='replace',              &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) THEN
         restart_write_txt = 1
         RETURN
      END IF
      
      restart_write_txt = 0

!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      WRITE (15,'(L1)') lasym
      WRITE (15,100) ns, mpol, ntor

      DO js=1,ns
         DO m =0,mpol
            DO n = -ntor,ntor
            WRITE (15,300) jbsupsmnsh(m,n,js), jbsupumnch(m,n,js),       &
                           jbsupvmnch(m,n,js), jpmnch(m,n,js)
            ENDDO
         ENDDO
      ENDDO

      IF (lasym) THEN

      DO js=1,ns
         DO m =0,mpol
            DO n = -ntor,ntor
            WRITE (15,300) jbsupsmnch(m,n,js), jbsupumnsh(m,n,js),       &
                           jbsupvmnsh(m,n,js), jpmnsh(m,n,js)
            ENDDO
         ENDDO
      ENDDO

      ENDIF

      CLOSE (UNIT=15)

 100  FORMAT(3i8)
 300  FORMAT(1p4e24.16)
      
      END FUNCTION restart_write_txt

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
      INTEGER, INTENT(in)      :: ns_old
      INTEGER, INTENT(in)      :: ns_new
      INTEGER, INTENT(in)      :: mpol_old
      INTEGER, INTENT(in)      :: mpol_new
      INTEGER, INTENT(in)      :: ntor_old
      INTEGER, INTENT(in)      :: ntor_new
      REAL(dp), INTENT(in)  :: aold(0:mpol_old,-ntor_old:ntor_old,ns_old)
      REAL(dp), INTENT(out) :: anew(0:mpol_new,-ntor_new:ntor_new,ns_new)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER   :: p5 = 0.5_dp
      REAL(dp), ALLOCATABLE :: temp1(:,:,:)
      INTEGER                  :: mpol_min, ntor_min
      INTEGER                  :: js_old, js_new, istat
      REAL(dp)              :: hs_old, hs_new, snew, sold, x1
!-----------------------------------------------

      mpol_min = MIN(mpol_old, mpol_new)
      ntor_min = MIN(ntor_old, ntor_new)

      anew = 0

      IF (ns_old .EQ. ns_new) THEN
         anew(0:mpol_min,-ntor_min:ntor_min,:) = aold(0:mpol_min,-ntor_min:ntor_min,:)
         RETURN
      END IF

      ALLOCATE (temp1(0:mpol_min,-ntor_min:ntor_min,ns_old), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error #1 in INTERPIT')
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
         CALL ASSERT(js_old.le.ns_old,'js_old >= ns_old in INTERPIT')
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
      ALLOCATE (temp1(0:mpol_min,-ntor_min:ntor_min,ns_new), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error #2 in INTERPIT')
      temp1 = anew
      anew(:,:,1) = 0
      anew(:,:,2:ns_new) = p5*(temp1(:,:,2:ns_new) + temp1(:,:,1:ns_new-1))
      DEALLOCATE(temp1)

      END SUBROUTINE interpit

      END MODULE
