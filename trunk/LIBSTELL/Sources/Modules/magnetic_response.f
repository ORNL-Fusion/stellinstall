!*******************************************************************************
!>  @file magnetic_response.f
!>  @brief Contains module @ref magnetic_response
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref magnetic_response_class.
!*******************************************************************************

      MODULE magnetic_response
      USE stel_kinds, only: rprec
#if defined (MPI_OPT)
      USE mpi
#endif
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  magnetic_response module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: magnetic_response_len = 80

!>  Bit position for the use coil response flag.
      INTEGER, PARAMETER :: magnetic_response_use_coil_flag   = 1
!>  Bit position for the force coil response flag.
      INTEGER, PARAMETER :: magnetic_response_use_plasma_flag = 2
!>  Bit position for the stellerator symmetry flag.
      INTEGER, PARAMETER :: magnetic_response_stell_sym_flag  = 3
!>  Bit position for the use conducting shell flag.
      INTEGER, PARAMETER :: magnetic_response_use_shell_flag  = 4

!  NETCDF file variable names --------------------------------------------------
!  Identification Variables ----------------------------------------------------
!>  NETCDF code name variable.
      CHARACTER (len=*), PARAMETER :: nc_name = 'mddc_mrf_code_name'
!>  NETCDF version variable.
      CHARACTER (len=*), PARAMETER :: nc_version =                             &
     &   'mddc_mrf_code_version'
!>  NETCDF date run variable.
      CHARACTER (len=*), PARAMETER :: nc_date = 'mddc_mrf_date_run'
!>  NETCDF coil identifier variable.
      CHARACTER (len=*), PARAMETER :: nc_coil_id =                             &
     &   'mddc_mrf_field_coils_id'

!  Coil Responce Function Variables --------------------------------------------
!>  NETCDF number of field coils variable.
      CHARACTER (len=*), PARAMETER :: nc_n_field_cg =                          &
     &   'mddc_mrf_n_field_cg'
!>  NETCDF mutual inductance variable.
      CHARACTER (len=*), PARAMETER :: nc_inductance =                          &
     &   'mddc_mrf_rdiag_coilg_1'
!>  NETCDF current scale variable.
      CHARACTER (len=*), PARAMETER :: nc_current_scale =                       &
     &   'mddc_mrf_extcur_mg'

!  Plasma Response Grid Variables ----------------------------------------------
!>  NETCDF number torodial planes variable.
      CHARACTER (len=*), PARAMETER :: nc_num_t = 'mddc_mrf_kp'
!>  NETCDF number radial grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_r = 'mddc_mrf_ir'
!>  NETCDF number vertical grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_z = 'mddc_mrf_jz'
!>  NETCDF minimum radial grid variable.
      CHARACTER (len=*), PARAMETER :: nc_rmin = 'mddc_mrf_rmin'
!>  NETCDF maximum radial grid variable.
      CHARACTER (len=*), PARAMETER :: nc_rmax = 'mddc_mrf_rmax'
!>  NETCDF minimum vertical grid variable.
      CHARACTER (len=*), PARAMETER :: nc_zmin = 'mddc_mrf_zmin'
!>  NETCDF maximum vertical grid variable.
      CHARACTER (len=*), PARAMETER :: nc_zmax = 'mddc_mrf_zmax'
!>  NETCDF number of field periods variable.
      CHARACTER (len=*), PARAMETER :: nc_n_field_periods =                     &
     &   'mddc_mrf_n_field_periods'

!>  NETCDF stell symmetric variable.
      CHARACTER (len=*), PARAMETER :: nc_stell_sym =                           &
     &   'mddc_mrf_lstell_sym'

!>  NETCDF radial response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_r = 'mddc_mrf_a_r'
!>  NETCDF toroidal response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_f = 'mddc_mrf_a_f'
!>  NETCDF vertical response grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_z = 'mddc_mrf_a_z'

!  Conducting Shell Response Arrays --------------------------------------------
!>
      CHARACTER (len=*), PARAMETER :: nc_use_shell =                           &
     &   'mddc_mrf_use_con_shell'
!>  NETCDF number torodial shell grid points variable.
      CHARACTER (len=*), PARAMETER :: nc_num_t_shell =                         &
     &   'mddc_mrf_kp_shell'

!>  NETCDF radial response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_r = 'mddc_mrf_a_s_r'
!>  NETCDF toroidal response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_f = 'mddc_mrf_a_s_f'
!>  NETCDF vertical response shell grid variable.
      CHARACTER (len=*), PARAMETER :: nc_a_s_z = 'mddc_mrf_a_s_z'

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a magnetic signal response function.
!-------------------------------------------------------------------------------
      TYPE magnetic_response_class
!  Identification Variables ----------------------------------------------------
!>  Name of the code which computed the responses.
         CHARACTER (len=magnetic_response_len) :: name
!>  Version of the code.
         CHARACTER (len=magnetic_response_len) :: version
!>  Data and time the code was run.
         CHARACTER (len=magnetic_response_len) :: date
!>  Idenifier of field-coils.
         CHARACTER (len=magnetic_response_len) :: coil_id

!  Coil Responce Function Variables --------------------------------------------
!>  Number of field-coil groups.
         INTEGER                             :: n_field_cg
!>  Array of diagnostic field-coil-group responses.
         REAL (rprec), DIMENSION(:), POINTER :: inductance => null()
!>  Array of external current scales.
         REAL (rprec), DIMENSION(:), POINTER :: current_scale => null()

!  Plasma Response Grid Variables ----------------------------------------------
!>  Number of grid points in R, plasma grid.
         INTEGER                            :: num_r
!>  Number of grid points in z, plasma grid.
         INTEGER                            :: num_z
!>  Number of phi planes per field period in plasma grid.
         INTEGER                            :: num_t
!>  Number of phi planes actually stored in plasma grid.
         INTEGER                            :: num_t_store
!>  Minimum R for plasma grid.
         REAL (rprec)                       :: rmin
!>  Maximum R for plasma grid.
         REAL (rprec)                       :: rmax
!>  Minimum z for plasma grid.
         REAL (rprec)                       :: zmin
!>  Maximum z for plasma grid.
         REAL (rprec)                       :: zmax

!>  Number of field periods.
         INTEGER                            :: n_field_periods
!>  Use stellarator symmetry.
         LOGICAL                            :: stell_sym

!>  Number of grid points in u, conducting shell grid.
         INTEGER                            :: num_u
!>  Number of phi planes per field period in the conducting shell grid.
         INTEGER                            :: num_t_shell
!>  Number of phi planes actually stored in the conducting shell grid.
         INTEGER                            :: num_t_shell_store

!  Plasma Response Arrays ------------------------------------------------------
!>  R component of plasma response function.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
!>  Phi component of plasma response function.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_f => null()
!>  Z component of plasma response function.
         REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z => null()

!  Conducting Shell Response Arrays --------------------------------------------
!>  R component of the conducting shell response function.
         REAL (rprec), DIMENSION(:,:), POINTER :: a_s_r => null()
!>  Phi component of the conducting shell response function.
         REAL (rprec), DIMENSION(:,:), POINTER :: a_s_f => null()
!>  Z component of the conducting shell response function.
         REAL (rprec), DIMENSION(:,:), POINTER :: a_s_z => null()

      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref magnetic_response_class types using
!>  @ref magnetic_response_construct_netcdf
!>
!>  @note @ref magnetic_construct_mpi cannot be included since netcdf handel and
!>  mpi comm have the dsame base type. Keep magnetic_construct for style.
!-------------------------------------------------------------------------------
      INTERFACE magnetic_response_construct
         MODULE PROCEDURE magnetic_response_construct_new,                     &
     &                    magnetic_response_construct_netcdf
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class object.
!>
!>  Allocates memory and initializes a @ref magnetic_response_class object.
!>
!>  @param[in] mdsig_iou       An instance of a the netcdf id of the open mdsig
!>                             file.
!>  @param[in] name            Name of the code which computed the responses.
!>  @param[in] version         Version of the code.
!>  @param[in] date            Data and time the code was run.
!>  @param[in] coil_id         Idenifier of field-coils.
!>  @param[in] inductance      Array of diagnostic field-coil-group responses.
!>  @param[in] current_scale   Array of external current scales.
!>  @param[in] num_t           Number of phi planes per field period in plasma
!>                             grid.
!>  @param[in] num_t_shell     Number of phi planes per field period in the
!>                             conducting shell grid.
!>  @param[in] rmin            Minimum R for plasma grid.
!>  @param[in] rmax            Maximum R for plasma grid.
!>  @param[in] zmin            Minimum z for plasma grid.
!>  @param[in] zmax            Maximum z for plasma grid.
!>  @param[in] n_field_periods Number of field periods.
!>  @param[in] stell_sym       Use stellarator symmetry.
!>  @param[in] a_r             R component of plasma response function.
!>  @param[in] a_f             Phi component of plasma response function.
!>  @param[in] a_z             Z component of plasma response function.
!>  @param[in] a_s_r           R component of the conducting shell response
!>                             function.
!>  @param[in] a_s_f           Phi component of the conducting shell response
!>                             function.
!>  @param[in] a_s_z           Z component of the conducting shell response
!>                             function.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_new(name, version, date,            &
     &                                         coil_id, inductance,            &
     &                                         current_scale,                  &
     &                                         num_t, num_t_shell,             &
     &                                         rmin, rmax, zmin, zmax,         &
     &                                         n_field_periods,                &
     &                                         stell_sym, a_r, a_f, a_z,       &
     &                                         a_s_r, a_s_f, a_s_z)
      USE ezcdf
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER                                  &
     &   :: magnetic_response_construct_new
      CHARACTER (len=*), INTENT(in)           :: name
      CHARACTER (len=*), INTENT(in)           :: version
      CHARACTER (len=*), INTENT(in)           :: date
      CHARACTER (len=*), INTENT(in)           :: coil_id
      REAL (rprec), DIMENSION(:), POINTER     :: inductance
      REAL (rprec), DIMENSION(:), POINTER     :: current_scale
      INTEGER, INTENT(in)                     :: num_t
      INTEGER, INTENT(in)                     :: num_t_shell
      REAL (rprec), INTENT(in)                :: rmin
      REAL (rprec), INTENT(in)                :: rmax
      REAL (rprec), INTENT(in)                :: zmin
      REAL (rprec), INTENT(in)                :: zmax
      INTEGER, INTENT(in)                     :: n_field_periods
      LOGICAL, INTENT(in)                     :: stell_sym
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_f
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_r
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_f
      REAL (rprec), DIMENSION(:,:), POINTER   :: a_s_z

!  local variables
      INTEGER, DIMENSION(3)                   :: dim_lengths
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_response_construct_new)

!  Identification Variables ----------------------------------------------------
      magnetic_response_construct_new%name = name
      magnetic_response_construct_new%version = version
      magnetic_response_construct_new%date = date
      magnetic_response_construct_new%coil_id = coil_id

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(inductance)) THEN
         ALLOCATE(magnetic_response_construct_new%inductance(                  &
     &               SIZE(inductance)))
         magnetic_response_construct_new%inductance = inductance
         magnetic_response_construct_new%n_field_cg = SIZE(inductance)
      END IF

      IF (ASSOCIATED(current_scale)) THEN
         ALLOCATE(magnetic_response_construct_new%current_scale(               &
     &               SIZE(current_scale)))
         magnetic_response_construct_new%current_scale = current_scale
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      magnetic_response_construct_new%num_t = num_t
      magnetic_response_construct_new%num_t_shell = num_t_shell
      magnetic_response_construct_new%rmin = rmin
      magnetic_response_construct_new%rmax = rmax
      magnetic_response_construct_new%zmin = zmin
      magnetic_response_construct_new%zmax = zmax
      magnetic_response_construct_new%n_field_periods = n_field_periods
      magnetic_response_construct_new%stell_sym = stell_sym

!  Plasma Response Arrays ------------------------------------------------------
      IF (ASSOCIATED(a_r)) THEN
         CALL assert(ASSOCIATED(a_f), 'a_f response function not ' //          &
     &                                'allocated')
         CALL assert(ASSOCIATED(a_z), 'a_z response function not ' //          &
     &                                'allocated')

         dim_lengths = SHAPE(a_r)
         magnetic_response_construct_new%num_r = dim_lengths(1)
         magnetic_response_construct_new%num_z = dim_lengths(2)
         magnetic_response_construct_new%num_t_store = dim_lengths(3)

         ALLOCATE(magnetic_response_construct_new%a_r(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         ALLOCATE(magnetic_response_construct_new%a_f(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         ALLOCATE(magnetic_response_construct_new%a_z(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         magnetic_response_construct_new%a_r = a_r
         magnetic_response_construct_new%a_f = a_f
         magnetic_response_construct_new%a_z = a_z
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (ASSOCIATED(a_s_r)) THEN
         CALL assert(ASSOCIATED(a_s_f), 'a_s_f response function ' //          &
     &                                  'not allocated')
         CALL assert(ASSOCIATED(a_s_z), 'a_s_z response function ' //          &
     &                                  'not allocated')

         dim_lengths(1:2) = SHAPE(a_s_r)
!         magnetic_response_construct_new%num_u = dim_lengths(1)
         magnetic_response_construct_new%num_t_shell_store =                   &
     &      dim_lengths(3)

         ALLOCATE(magnetic_response_construct_new%a_s_r(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         ALLOCATE(magnetic_response_construct_new%a_s_f(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         ALLOCATE(magnetic_response_construct_new%a_s_z(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         magnetic_response_construct_new%a_s_r = a_s_r
         magnetic_response_construct_new%a_s_f = a_s_f
         magnetic_response_construct_new%a_s_z = a_s_z
      END IF

      CALL profiler_set_stop_time('magnetic_response_construct_new',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class object.
!>
!>  Allocates memory and initializes a @ref magnetic_class object from a mdsig
!>  file.
!>
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!>  @returns A pointer to a constructed @ref magnetic_response_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_netcdf(mdsig_iou)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER ::                               &
     &   magnetic_response_construct_netcdf
      INTEGER, INTENT(in)                     :: mdsig_iou

!  local variables
      LOGICAL                                 :: use_inductance
      LOGICAL                                 :: use_response
      LOGICAL                                 :: use_shell
      INTEGER, DIMENSION(3)                   :: dim_lengths
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_response_construct_netcdf)

!  Identification Variables ----------------------------------------------------
      CALL cdf_read(mdsig_iou, nc_name,                                        &
     &              magnetic_response_construct_netcdf%name)
      CALL cdf_read(mdsig_iou, nc_version,                                     &
     &              magnetic_response_construct_netcdf%version)
      CALL cdf_read(mdsig_iou, nc_date,                                        &
     &              magnetic_response_construct_netcdf%date)
      CALL cdf_read(mdsig_iou, nc_coil_id,                                     &
     &              magnetic_response_construct_netcdf%coil_id)

!  Check the version code verison. If the version is 'MRC 2014-09-28' or later,
!  there is a flag variable that indicates the contents of the file.
      SELECT CASE (TRIM(magnetic_response_construct_netcdf%version))

         CASE ('MRC 2014-09-28')
            CALL cdf_read(mdsig_iou, 'mddc_mrf_use_con_shell',                 &
     &                    use_shell)
            use_inductance = .true.
            use_response = .true.

         CASE DEFAULT
            use_inductance = .true.
            use_response = .true.
            use_shell = .false.

      END SELECT

!  Coil Responce Function Variables --------------------------------------------
      IF (use_inductance) THEN
         CALL cdf_read(mdsig_iou, nc_n_field_cg,                               &
     &                 magnetic_response_construct_netcdf%n_field_cg)

         ALLOCATE(magnetic_response_construct_netcdf%inductance(               &
     &               magnetic_response_construct_netcdf%n_field_cg))
         CALL cdf_read(mdsig_iou, nc_inductance,                               &
     &                 magnetic_response_construct_netcdf%inductance)

         ALLOCATE(magnetic_response_construct_netcdf%current_scale(            &
     &               magnetic_response_construct_netcdf%n_field_cg))
         CALL cdf_read(mdsig_iou, nc_current_scale,                            &
     &                 magnetic_response_construct_netcdf%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (use_response) THEN
         CALL cdf_read(mdsig_iou, nc_num_t,                                    &
     &                 magnetic_response_construct_netcdf%num_t)
         CALL cdf_read(mdsig_iou, nc_rmin,                                     &
     &                 magnetic_response_construct_netcdf%rmin)
         CALL cdf_read(mdsig_iou, nc_rmax,                                     &
     &                 magnetic_response_construct_netcdf%rmax)
         CALL cdf_read(mdsig_iou, nc_zmin,                                     &
     &                 magnetic_response_construct_netcdf%zmin)
         CALL cdf_read(mdsig_iou, nc_zmax,                                     &
     &                 magnetic_response_construct_netcdf%zmax)
         CALL cdf_read(mdsig_iou, nc_n_field_periods,                          &
     &           magnetic_response_construct_netcdf%n_field_periods)
         CALL cdf_read(mdsig_iou, nc_stell_sym,                                &
     &                 magnetic_response_construct_netcdf%stell_sym)

         CALL cdf_inquire(mdsig_iou, nc_a_r, dim_lengths)
         magnetic_response_construct_netcdf%num_r = dim_lengths(1)
         magnetic_response_construct_netcdf%num_z = dim_lengths(2)
         magnetic_response_construct_netcdf%num_t_store = dim_lengths(3)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_r(dim_lengths(1),             &
     &                                             dim_lengths(2),             &
     &                                             dim_lengths(3)))
         CALL cdf_read(mdsig_iou, nc_a_r,                                      &
     &                 magnetic_response_construct_netcdf%a_r)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_f(dim_lengths(1),             &
     &                                             dim_lengths(2),             &
     &                                             dim_lengths(3)))
         CALL cdf_read(mdsig_iou, nc_a_f,                                      &
     &                 magnetic_response_construct_netcdf%a_f)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_z(dim_lengths(1),             &
     &                                             dim_lengths(2),             &
     &                                             dim_lengths(3)))
         CALL cdf_read(mdsig_iou, nc_a_z,                                      &
     &                 magnetic_response_construct_netcdf%a_z)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (use_shell) THEN
         IF (.not.use_response) THEN
            CALL cdf_read(mdsig_iou, nc_stell_sym,                             &
     &                    magnetic_response_construct_netcdf%stell_sym)
         END IF

         CALL cdf_read(mdsig_iou, nc_num_t_shell,                              &
     &                 magnetic_response_construct_netcdf%num_t_shell)

         CALL cdf_inquire(mdsig_iou, nc_a_s_r, dim_lengths(1:2))
         magnetic_response_construct_netcdf%num_u = dim_lengths(1)
         magnetic_response_construct_netcdf%num_t_shell_store =                &
     &      dim_lengths(2)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_s_r(dim_lengths(1),           &
     &                                               dim_lengths(2)))
         CALL cdf_read(mdsig_iou, nc_a_s_r,                                    &
     &                 magnetic_response_construct_netcdf%a_s_r)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_s_f(dim_lengths(1),           &
     &                                               dim_lengths(2)))
         CALL cdf_read(mdsig_iou, nc_a_s_f,                                    &
     &                 magnetic_response_construct_netcdf%a_s_f)

         ALLOCATE(                                                             &
     &      magnetic_response_construct_netcdf%a_s_z(dim_lengths(1),           &
     &                                               dim_lengths(2)))
         CALL cdf_read(mdsig_iou, nc_a_s_z,                                    &
     &                 magnetic_response_construct_netcdf%a_s_z)
      END IF

      CALL profiler_set_stop_time('magnetic_response_construct_netcdf',        &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_response_class.
!>
!>  Allocates memory and initializes a @ref magnetic_response_class object by
!>  syncing through mpi.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref magnetic_response_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_response_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER ::                               &
     &   magnetic_response_construct_mpi
      INTEGER, INTENT(in)            :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                        :: flags
      INTEGER                        :: error
      INTEGER, DIMENSION(3)          :: dim_lengths
      INTEGER                        :: total_length

!  Start of executable code
      ALLOCATE(magnetic_response_construct_mpi)

!  All state flags are bit packed into a single integer to reduce mpi
!  communication overhead.
      CALL MPI_BCAST(flags, 1, MPI_INTEGER, 0, intra_comm, error)
      magnetic_response_construct_mpi%stell_sym =                              &
     &   BTEST(flags, magnetic_response_stell_sym_flag)

!  Identification Variables ----------------------------------------------------
      CALL MPI_BCAST(magnetic_response_construct_mpi%name,                     &
     &               magnetic_response_len, MPI_CHARACTER, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_response_construct_mpi%version,                  &
     &               magnetic_response_len, MPI_CHARACTER, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_response_construct_mpi%date,                     &
     &               magnetic_response_len, MPI_CHARACTER, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_response_construct_mpi%coil_id,                  &
     &               magnetic_response_len, MPI_CHARACTER, 0,                  &
     &               intra_comm, error)

!  Coil Responce Function Variables --------------------------------------------
      IF (BTEST(flags, magnetic_response_use_coil_flag)) THEN
         CALL MPI_BCAST(magnetic_response_construct_mpi%n_field_cg, 1,         &
     &                  MPI_INTEGER, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%inductance(                  &
     &               magnetic_response_construct_mpi%n_field_cg))
         CALL MPI_BCAST(magnetic_response_construct_mpi%inductance,            &
     &                  magnetic_response_construct_mpi%n_field_cg,            &
     &                  MPI_REAL8, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%current_scale(               &
     &               magnetic_response_construct_mpi%n_field_cg))
         CALL MPI_BCAST(magnetic_response_construct_mpi%current_scale,         &
     &                  magnetic_response_construct_mpi%n_field_cg,            &
     &                  MPI_REAL8, 0, intra_comm, error)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (BTEST(flags, magnetic_response_use_plasma_flag)) THEN
         CALL MPI_BCAST(magnetic_response_construct_mpi%num_t, 1,              &
     &                  MPI_INTEGER, 0, intra_comm, error)
         CALL MPI_BCAST(magnetic_response_construct_mpi%rmin, 1,               &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(magnetic_response_construct_mpi%rmax, 1,               &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(magnetic_response_construct_mpi%zmin, 1,               &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(magnetic_response_construct_mpi%zmax, 1,               &
     &                  MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(magnetic_response_construct_mpi%n_field_periods,       &
     &                  1, MPI_INTEGER, 0, intra_comm, error)

         CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm,            &
     &                  error)
         magnetic_response_construct_mpi%num_r = dim_lengths(1)
         magnetic_response_construct_mpi%num_z = dim_lengths(2)
         magnetic_response_construct_mpi%num_t_store = dim_lengths(3)

         total_length = dim_lengths(1)*dim_lengths(2)*dim_lengths(3)

         ALLOCATE(magnetic_response_construct_mpi%a_r(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_r,                   &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%a_f(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_f,                   &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%a_z(dim_lengths(1),          &
     &                                                dim_lengths(2),          &
     &                                                dim_lengths(3)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_z,                   &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (BTEST(flags, magnetic_response_use_shell_flag)) THEN
         CALL MPI_BCAST(magnetic_response_construct_mpi%num_t_shell,           &
     &                  1, MPI_INTEGER, 0, intra_comm, error)

         CALL MPI_BCAST(dim_lengths(1:2), 2, MPI_INTEGER, 0, intra_comm,       &
     &                  error)
!         magnetic_response_construct_mpi%num_u = dim_lengths(1)
         magnetic_response_construct_mpi%num_t_shell_store =                   &
     &      dim_lengths(2)

         total_length = dim_lengths(1)*dim_lengths(2)

         ALLOCATE(magnetic_response_construct_mpi%a_s_r(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_s_r,                 &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%a_s_f(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_s_f,                 &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)

         ALLOCATE(magnetic_response_construct_mpi%a_s_z(dim_lengths(1),        &
     &                                                  dim_lengths(2)))
         CALL MPI_BCAST(magnetic_response_construct_mpi%a_s_z,                 &
     &                  total_length, MPI_REAL8, 0, intra_comm, error)
      END IF
#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref magnetic_response_class object.
!>
!>  Deallocates memory and uninitializes a @ref magnetic_response_class object.
!>
!>  @param[inout] this A @ref magnetic_response_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_destruct(this)

!  Declare Arguments
      TYPE (magnetic_response_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%inductance)) THEN
         DEALLOCATE(this%inductance)
         this%inductance => null()
      END IF

      IF (ASSOCIATED(this%current_scale)) THEN
         DEALLOCATE(this%current_scale)
         this%current_scale => null()
      END IF

      IF (ASSOCIATED(this%a_r)) THEN
         DEALLOCATE(this%a_r)
         this%a_r => null()
      END IF

      IF (ASSOCIATED(this%a_f)) THEN
         DEALLOCATE(this%a_f)
         this%a_f => null()
      END IF

      IF (ASSOCIATED(this%a_z)) THEN
         DEALLOCATE(this%a_z)
         this%a_z => null()
      END IF

      IF (ASSOCIATED(this%a_s_r)) THEN
         DEALLOCATE(this%a_s_r)
         this%a_s_r => null()
      END IF

      IF (ASSOCIATED(this%a_s_f)) THEN
         DEALLOCATE(this%a_s_f)
         this%a_s_f => null()
      END IF

      IF (ASSOCIATED(this%a_s_z)) THEN
         DEALLOCATE(this%a_s_z)
         this%a_s_z => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE
!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Defines the variables for the NETCDF file.
!>
!>  Write out the magnetic response to an mdsig file.
!>
!>  @param[in] this      A @ref magnetic_response_class instance.
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_define(this, mdsig_iou)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(in) :: this
      INTEGER, INTENT(in)                        :: mdsig_iou

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Define Variables ------------------------------------------------------------
!  Identification Variables ----------------------------------------------------
      CALL cdf_define(mdsig_iou, nc_name, this%name)
      CALL cdf_define(mdsig_iou, nc_version, this%version)
      CALL cdf_define(mdsig_iou, nc_date, this%date)
      CALL cdf_define(mdsig_iou, nc_coil_id, this%coil_id)

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(this%inductance)) THEN
         CALL cdf_define(mdsig_iou, nc_n_field_cg, this%n_field_cg)
         CALL cdf_define(mdsig_iou, nc_inductance, this%inductance)
         CALL cdf_define(mdsig_iou, nc_current_scale,                          &
     &                   this%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (ASSOCIATED(this%a_r)) THEN
         CALL cdf_define(mdsig_iou, nc_num_t, this%num_t)
         CALL cdf_define(mdsig_iou, nc_num_r, this%num_r)
         CALL cdf_define(mdsig_iou, nc_num_z, this%num_z)
         CALL cdf_define(mdsig_iou, nc_rmin, this%rmin)
         CALL cdf_define(mdsig_iou, nc_rmax, this%rmax)
         CALL cdf_define(mdsig_iou, nc_zmin, this%zmin)
         CALL cdf_define(mdsig_iou, nc_zmax, this%zmax)
         CALL cdf_define(mdsig_iou, nc_n_field_periods,                        &
     &                   this%n_field_periods)
         CALL cdf_define(mdsig_iou, nc_stell_sym, this%stell_sym)
         CALL cdf_define(mdsig_iou, nc_a_r, this%a_r)
         CALL cdf_define(mdsig_iou, nc_a_f, this%a_f)
         CALL cdf_define(mdsig_iou, nc_a_z, this%a_z)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      CALL cdf_define(mdsig_iou, nc_use_shell, ASSOCIATED(this%a_s_r))
      IF (ASSOCIATED(this%a_s_r)) THEN
         IF (.not.ASSOCIATED(this%a_r)) THEN
            CALL cdf_define(mdsig_iou, nc_stell_sym, this%stell_sym)
         END IF
         CALL cdf_define(mdsig_iou, nc_use_shell, .true.)
         CALL cdf_define(mdsig_iou, nc_num_t_shell, this%num_t_shell)
         CALL cdf_define(mdsig_iou, nc_a_s_r, this%a_s_r)
         CALL cdf_define(mdsig_iou, nc_a_s_f, this%a_s_f)
         CALL cdf_define(mdsig_iou, nc_a_s_z, this%a_s_z)
      END IF

      CALL profiler_set_stop_time('magnetic_response_define',                  &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write variables to the NETCDF file.
!>
!>  Write out the magnetic response to an mdsig file.
!>
!>  @param[in] this      A @ref magnetic_response_class instance.
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_write(this, mdsig_iou)
      USE ezcdf
      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(in) :: this
      INTEGER, INTENT(in)                        :: mdsig_iou

!  local variables
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Define Variables ------------------------------------------------------------
!  Identification Variables ----------------------------------------------------
      CALL cdf_write(mdsig_iou, nc_name, this%name)
      CALL cdf_write(mdsig_iou, nc_version, this%version)
      CALL cdf_write(mdsig_iou, nc_date, this%date)
      CALL cdf_write(mdsig_iou, nc_coil_id, this%coil_id)

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(this%inductance)) THEN
         CALL cdf_write(mdsig_iou, nc_n_field_cg, this%n_field_cg)
         CALL cdf_write(mdsig_iou, nc_inductance, this%inductance)
         CALL cdf_write(mdsig_iou, nc_current_scale, this%current_scale)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (ASSOCIATED(this%a_r)) THEN
         CALL cdf_write(mdsig_iou, nc_num_t, this%num_t)
         CALL cdf_write(mdsig_iou, nc_num_r, this%num_r)
         CALL cdf_write(mdsig_iou, nc_num_z, this%num_z)
         CALL cdf_write(mdsig_iou, nc_rmin, this%rmin)
         CALL cdf_write(mdsig_iou, nc_rmax, this%rmax)
         CALL cdf_write(mdsig_iou, nc_zmin, this%zmin)
         CALL cdf_write(mdsig_iou, nc_zmax, this%zmax)
         CALL cdf_write(mdsig_iou, nc_n_field_periods,                         &
     &                  this%n_field_periods)
         CALL cdf_write(mdsig_iou, nc_stell_sym, this%stell_sym)
         CALL cdf_write(mdsig_iou, nc_a_r, this%a_r)
         CALL cdf_write(mdsig_iou, nc_a_f, this%a_f)
         CALL cdf_write(mdsig_iou, nc_a_z, this%a_z)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      CALL cdf_write(mdsig_iou, nc_use_shell, ASSOCIATED(this%a_s_r))
      IF (ASSOCIATED(this%a_s_r)) THEN
         IF (.not.ASSOCIATED(this%a_r)) THEN
            CALL cdf_write(mdsig_iou, nc_stell_sym, this%stell_sym)
         END IF
         CALL cdf_write(mdsig_iou, nc_use_shell, .true.)
         CALL cdf_write(mdsig_iou, nc_num_t_shell, this%num_t_shell)
         CALL cdf_write(mdsig_iou, nc_a_s_r, this%a_s_r)
         CALL cdf_write(mdsig_iou, nc_a_s_f, this%a_s_f)
         CALL cdf_write(mdsig_iou, nc_a_s_z, this%a_s_z)
      END IF

      CALL profiler_set_stop_time('magnetic_response_write', start_time)

      END SUBROUTINE

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref magnetic_response_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_response_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_response_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                           :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                                       :: flags
      INTEGER                                       :: error
      INTEGER, DIMENSION(3)                         :: dim_lengths
      INTEGER                                       :: total_length
      REAL (rprec)                                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  All state flags are bit packed into a single integer to reduce mpi
!  communication overhead.
      flags = 0
      IF (ASSOCIATED(this%inductance)) THEN
         flags = IBSET(flags, magnetic_response_use_coil_flag)
      END IF
      IF (ASSOCIATED(this%a_r)) THEN
         flags = IBSET(flags, magnetic_response_use_plasma_flag)
      END IF
      IF (this%stell_sym) THEN
         flags = IBSET(flags, magnetic_response_stell_sym_flag)
      END IF
      IF (ASSOCIATED(this%a_s_r)) THEN
         flags = IBSET(flags, magnetic_response_use_shell_flag)
      END IF
      CALL MPI_BCAST(flags, 1, MPI_INTEGER, 0, intra_comm, error)

!  Identification Variables ----------------------------------------------------
      CALL MPI_BCAST(this%name, magnetic_response_len, MPI_CHARACTER, 0,       &
     &               intra_comm, error)
      CALL MPI_BCAST(this%version, magnetic_response_len, MPI_CHARACTER,       &
     &               0, intra_comm, error)
      CALL MPI_BCAST(this%date, magnetic_response_len, MPI_CHARACTER, 0,       &
     &               intra_comm, error)
      CALL MPI_BCAST(this%coil_id, magnetic_response_len, MPI_CHARACTER,       &
     &               0, intra_comm, error)

!  Coil Responce Function Variables --------------------------------------------
      IF (ASSOCIATED(this%inductance)) THEN
         CALL MPI_BCAST(this%n_field_cg, 1, MPI_INTEGER, 0, intra_comm,        &
     &                  error)
         CALL MPI_BCAST(this%inductance, this%n_field_cg, MPI_REAL8, 0,        &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%current_scale, this%n_field_cg, MPI_REAL8,        &
     &                  0, intra_comm, error)
      END IF

!  Plasma Response Grid Variables ----------------------------------------------
      IF (ASSOCIATED(this%a_r)) THEN
         CALL MPI_BCAST(this%num_t, 1, MPI_INTEGER, 0, intra_comm,             &
     &                  error)
         CALL MPI_BCAST(this%rmin, 1, MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%rmax, 1, MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%zmin, 1, MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%zmax, 1, MPI_REAL8, 0, intra_comm, error)
         CALL MPI_BCAST(this%n_field_periods, 1, MPI_INTEGER, 0,               &
     &                  intra_comm, error)

         dim_lengths = SHAPE(this%a_r)
         CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm,            &
     &                  error)

         total_length = dim_lengths(1)*dim_lengths(2)*dim_lengths(3)

         CALL MPI_BCAST(this%a_r, total_length, MPI_REAL8, 0,                  &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%a_f, total_length, MPI_REAL8, 0,                  &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%a_z, total_length, MPI_REAL8, 0,                  &
     &                  intra_comm, error)
      END IF

!  Conducting Shell Response Arrays --------------------------------------------
      IF (ASSOCIATED(this%a_s_r)) THEN
         CALL MPI_BCAST(this%num_t_shell, 1, MPI_INTEGER, 0, intra_comm,       &
     &                  error)

         dim_lengths(1:2) = SHAPE(this%a_s_r)
         CALL MPI_BCAST(dim_lengths(1:2), 2, MPI_INTEGER, 0, intra_comm,       &
     &                  error)

         total_length = dim_lengths(1)*dim_lengths(2)

         CALL MPI_BCAST(this%a_s_r, total_length, MPI_REAL8, 0,                &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%a_s_f, total_length, MPI_REAL8, 0,                &
     &                  intra_comm, error)
         CALL MPI_BCAST(this%a_s_z, total_length, MPI_REAL8, 0,                &
     &                  intra_comm, error)
      END IF

      CALL profiler_set_stop_time('magnetic_response_sync_children',           &
     &                            start_time)
#endif

      END SUBROUTINE

      END MODULE
