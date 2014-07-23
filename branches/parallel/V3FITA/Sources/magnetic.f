!*******************************************************************************
!>  @file magnetic.f
!>  @brief Contains module @ref magnetic
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref magnetic_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE magnetic
      USE stel_kinds, only: rprec
      USE data_parameters
#if defined (MPI_OPT)
      USE mpi
#endif

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a magnetic signal.
!>  @par Super Class:
!>  @ref diagnostic
!>
!>  @note Copied from @ref mddc_mrf.
!-------------------------------------------------------------------------------
      TYPE magnetic_class
!  Identification Variables ----------------------------------------------------
!>  Name of the code which computed the responses.
         CHARACTER(len=data_name_length)    :: code_name
!>  Version of the code.
         CHARACTER(len=data_name_length)    :: code_version
!>  Data and time the code was run.
         CHARACTER(len=data_name_length)    :: date_run
!>  Idenifier of field-coilds.
         CHARACTER(len=data_name_length)    :: field_coils_id

!  Coil Responce Function Variables --------------------------------------------
!>  Number of field-coil groups.
         INTEGER                            :: n_field_cg
!>  Array of diagnostic field-coil-group responses.
         REAL(rprec), DIMENSION(:), POINTER :: rdiag_coilg_1 => null()
!>  Array of external currents.
         REAL(rprec), DIMENSION(:), POINTER :: extcur_mg => null()

!  Plasma Response Grid Variables ----------------------------------------------
!>  Number of grid points in R, plasma grid.
         INTEGER                            :: ir
!>  Number of grid points in z, plasma grid.
         INTEGER                            :: jz
!>  Number of phi planes per field period in plasma grid.
         INTEGER                            :: kp
!>  Number of phi planes actually stored in plasma grid.
         INTEGER                            :: kp_store
!>  Minimum R for plasma grid.
         REAL(rprec)                        :: rmin
!>  Maximum R for plasma grid.
         REAL(rprec)                        :: rmax
!>  Minimum z for plasma grid.
         REAL(rprec)                        :: zmin
!>  Maximum z for plasma grid.
         REAL(rprec)                        :: zmax
!>  Number of field periods.
         INTEGER                            :: n_field_periods
!>  Toggles if the response is stellarator symmetric.
         LOGICAL                            :: lstell_sym

!  Plasma Response Arrays ------------------------------------------------------
!>  R component of plasma response function.
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
!>  Phi component of plasma response function.
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_f => null()
!>  Z component of plasma response function.
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_z => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref magnetic_class types using
!>  @ref magnetic_construct_netcdf
!>
!>  @note @ref magnetic_construct_mpi cannot be included since netcdf handel and
!>  mpi comm have the dsame base type. Keep magnetic_construct for style.
!-------------------------------------------------------------------------------
      INTERFACE magnetic_construct
         MODULE PROCEDURE magnetic_construct_netcdf
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_class object containing a @ref mddc_mrf
!>  object.
!>
!>  Allocates memory and initializes a @ref magnetic_class object.
!>
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_netcdf(mdsig_iou)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: magnetic_construct_netcdf
      INTEGER, INTENT(in)            :: mdsig_iou

!  local variables
      INTEGER, DIMENSION(3)          :: dim_lengths

!  Start of executable code
      ALLOCATE(magnetic_construct_netcdf)

!  Identification Variables ----------------------------------------------------
      CALL cdf_read(mdsig_iou, 'mddc_mrf_code_name',                           &
     &              magnetic_construct_netcdf%code_name)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_code_version',                        &
     &              magnetic_construct_netcdf%code_version)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_date_run',                            &
     &              magnetic_construct_netcdf%date_run)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_field_coils_id',                      &
     &              magnetic_construct_netcdf%field_coils_id)

!  Coil Responce Function Variables --------------------------------------------
      CALL cdf_read(mdsig_iou, 'mddc_mrf_n_field_cg',                          &
     &              magnetic_construct_netcdf%n_field_cg)

      ALLOCATE(magnetic_construct_netcdf%rdiag_coilg_1(                        &
     &            magnetic_construct_netcdf%n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rdiag_coilg_1',                       &
     &              magnetic_construct_netcdf%rdiag_coilg_1)

      ALLOCATE(magnetic_construct_netcdf%extcur_mg(                            &
     &            magnetic_construct_netcdf%n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_extcur_mg',                           &
     &              magnetic_construct_netcdf%extcur_mg)

!  Plasma Response Grid Variables ----------------------------------------------
      CALL cdf_read(mdsig_iou, 'mddc_mrf_kp',                                  &
     &              magnetic_construct_netcdf%kp)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmin',                                &
     &              magnetic_construct_netcdf%rmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmax',                                &
     &              magnetic_construct_netcdf%rmax)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmin',                                &
     &              magnetic_construct_netcdf%zmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmax',                                &
     &              magnetic_construct_netcdf%zmax)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_n_field_periods',                     &
     &              magnetic_construct_netcdf%n_field_periods)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_lstell_sym',                          &
     &              magnetic_construct_netcdf%lstell_sym)

!  Plasma Response Arrays ------------------------------------------------------
      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_r', dim_lengths)
      magnetic_construct_netcdf%ir = dim_lengths(1)
      magnetic_construct_netcdf%jz = dim_lengths(2)
      magnetic_construct_netcdf%kp_store = dim_lengths(3)

      ALLOCATE(magnetic_construct_netcdf%a_r(dim_lengths(1),                   &
     &                                       dim_lengths(2),                   &
     &                                       dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_r',                                 &
     &              magnetic_construct_netcdf%a_r)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_f', dim_lengths)
      ALLOCATE(magnetic_construct_netcdf%a_f(dim_lengths(1),                   &
     &                                       dim_lengths(2),                   &
     &                                       dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_f',                                 &
     &              magnetic_construct_netcdf%a_f)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_z', dim_lengths)
      ALLOCATE(magnetic_construct_netcdf%a_z(dim_lengths(1),                   &
     &                                       dim_lengths(2),                   &
     &                                       dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_z',                                 &
     &              magnetic_construct_netcdf%a_z)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_class by syncing through mpi.
!>
!>  Allocates memory and initializes a @ref equilibrium_class object.
!>
!>  @param[in] intra_comm A MPI intra_comm handle.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_mpi(intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: magnetic_construct_mpi
      INTEGER, INTENT(in)            :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                        :: error
      INTEGER, DIMENSION(3)          :: dim_lengths
      INTEGER                        :: total_length

!  Start of executable code
      ALLOCATE(magnetic_construct_mpi)

!  Identification Variables ----------------------------------------------------
      CALL MPI_BCAST(magnetic_construct_mpi%code_name,                         &
     &               data_name_length, MPI_CHARACTER, 0, intra_comm,           &
     &               error)
      CALL MPI_BCAST(magnetic_construct_mpi%code_version,                      &
     &               data_name_length, MPI_CHARACTER, 0, intra_comm,           &
     &               error)
      CALL MPI_BCAST(magnetic_construct_mpi%date_run,                          &
     &               data_name_length, MPI_CHARACTER, 0, intra_comm,           &
     &               error)
      CALL MPI_BCAST(magnetic_construct_mpi%field_coils_id,                    &
     &               data_name_length, MPI_CHARACTER, 0, intra_comm,           &
     &               error)

!  Coil Responce Function Variables --------------------------------------------
      CALL MPI_BCAST(magnetic_construct_mpi%n_field_cg, 1,                     &
     &               MPI_INTEGER, 0, intra_comm, error)

      ALLOCATE(magnetic_construct_mpi%rdiag_coilg_1(                           &
     &            magnetic_construct_mpi%n_field_cg))
      CALL MPI_BCAST(magnetic_construct_mpi%rdiag_coilg_1,                     &
     &               magnetic_construct_mpi%n_field_cg, MPI_REAL8, 0,          &
     &               intra_comm, error)

      ALLOCATE(magnetic_construct_mpi%extcur_mg(                               &
     &            magnetic_construct_mpi%n_field_cg))
      CALL MPI_BCAST(magnetic_construct_mpi%extcur_mg,                         &
     &               magnetic_construct_mpi%n_field_cg, MPI_REAL8, 0,          &
     &               intra_comm, error)

!  Plasma Response Grid Variables ----------------------------------------------
      CALL MPI_BCAST(magnetic_construct_mpi%kp, 1, MPI_INTEGER, 0,             &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%rmin, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%rmax, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%zmin, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%rmax, 1, MPI_INTEGER, 0,           &
     &               intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%n_field_periods, 1,                &
     &               MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(magnetic_construct_mpi%lstell_sym, 1, MPI_LOGICAL,        &
     &               0, intra_comm, error)

!  Plasma Response Arrays ------------------------------------------------------
      CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm, error)
      magnetic_construct_mpi%ir = dim_lengths(1)
      magnetic_construct_mpi%jz = dim_lengths(2)
      magnetic_construct_mpi%kp_store = dim_lengths(3)

      total_length = dim_lengths(1)*dim_lengths(2)*dim_lengths(3)

      ALLOCATE(magnetic_construct_mpi%a_r(dim_lengths(1),                      &
     &                                    dim_lengths(2),                      &
     &                                    dim_lengths(3)))
      CALL MPI_BCAST(magnetic_construct_mpi%a_r, total_length,                 &
     &               MPI_REAL8, 0, intra_comm, error)

      ALLOCATE(magnetic_construct_mpi%a_f(dim_lengths(1),                      &
     &                                    dim_lengths(2),                      &
     &                                    dim_lengths(3)))
      CALL MPI_BCAST(magnetic_construct_mpi%a_f, total_length,                 &
     &               MPI_REAL8, 0, intra_comm, error)

      ALLOCATE(magnetic_construct_mpi%a_z(dim_lengths(1),                      &
     &                                    dim_lengths(2),                      &
     &                                    dim_lengths(3)))
      CALL MPI_BCAST(magnetic_construct_mpi%a_z, total_length,                 &
     &               MPI_REAL8, 0, intra_comm, error)

#endif

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref magnetic_class object.
!>
!>  Deallocates memory and uninitializes a @ref magnetic_class object.
!>
!>  @param[inout] this A @ref magnetic_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%rdiag_coilg_1)) THEN
         DEALLOCATE(this%rdiag_coilg_1)
         this%rdiag_coilg_1 => null()
      END IF

      IF (ASSOCIATED(this%extcur_mg)) THEN
         DEALLOCATE(this%extcur_mg)
         this%extcur_mg => null()
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

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  The plasma signal is the volume integral of the current denisty dotted with
!>  the diagnostic response function vector. The response_function contains a r z
!>  grid of reponce function vectors. These grid points are then interpolated
!>  into the equilibrium grid positions for a single toroidal plane.
!>  Hishman et. al. Phys. Plasmas 11, 595 (2004); doi: 10.1063/1.1637347
!>
!>  phi = Int[dV J_plasma * A_diag]
!>
!>  J_plasma and dV are obtained from the @ref equilibrium. A_diag is obtained
!>  from the @ref mddc_mrf
!-------------------------------------------------------------------------------
!>
!>  External current contribrutions are included by the dot product of the
!>  external currents with mutual inductance values. These currents are obtained
!>  from @ref equilibrium if they exist.
!-------------------------------------------------------------------------------
!>
!>  The plasma only signal is stored in the 2nd index.
!>  The external current signal is stored in the 3rd index.
!>  Total signal is the sum of plasma signal and external signal.
!>
!>  @param[in]    this    A @ref magnetic_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_modeled_signal(this, a_model, sigma)
      USE model
      USE bivariate
      USE v3_utilities, only: err_warn

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal
      TYPE (magnetic_class), INTENT(in)         :: this
      TYPE (model_class), INTENT(in)            :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out)   :: sigma

!  Local variables
      INTEGER                                   :: phi, num_phi, r, z
      INTEGER                                   :: numrvol, numzvol
      TYPE (bivariate_type)                     :: bivariate_object
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: rgrid
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: zgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: rvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: zvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jrvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jphivolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jzvolgrid
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_rvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_phivol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_zvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sumrz
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sumphi

      LOGICAL                                   :: scale_current
      INTEGER                                   :: num_currents
      REAL (rprec), DIMENSION(:), POINTER       :: extcur

!  Start of executable code
      magnetic_get_modeled_signal = 0.0

!  Calculate the plasma only signal from the response function.
!  Get the volume integration grid points.
!  NOTE: Do not deallocate any of these arrays. They are owned by the
!  equilibrium
      rvolgrid =>                                                              &
     &   equilibrium_get_magnetic_volume_rgrid(a_model%equilibrium)
      zvolgrid =>                                                              &
     &   equilibrium_get_magnetic_volume_zgrid(a_model%equilibrium)
      jrvolgrid =>                                                             &
     &   equilibrium_get_magnetic_volume_jrgrid(a_model%equilibrium)
      jphivolgrid =>                                                           &
     &   equilibrium_get_magnetic_volume_jphigrid(a_model%equilibrium)
      jzvolgrid =>                                                             &
     &   equilibrium_get_magnetic_volume_jzgrid(a_model%equilibrium)

!  The equilibrium may return null pointers indicating the lack of a plasma.
      IF (ASSOCIATED(rvolgrid) .and. ASSOCIATED(zvolgrid) .and.                &
     &    ASSOCIATED(jrvolgrid) .and. ASSOCIATED(jphivolgrid) .and.            &
     &    ASSOCIATED(jzvolgrid)) THEN

!  Allocate space for the interpolated response functions. Only need two
!  dimensions because the interpolation is performed one phi position at a time.
         numrvol = SIZE(rvolgrid, 1)
         numzvol = SIZE(rvolgrid, 2)
         ALLOCATE(response_rvol(numrvol, numzvol))
         ALLOCATE(response_phivol(numrvol, numzvol))
         ALLOCATE(response_zvol(numrvol, numzvol))

!  Allocate work arrays to hold the results of the integration.
         num_phi = this%kp_store
         ALLOCATE(sumrz(numrvol, numzvol))
         ALLOCATE(sumphi(num_phi))

!  Create the r and z gridpoints for the magnetic response function.
!>  @todo FIXME: For a single diagnostic these never change. Allocate these once
!>  and cache these in the @ref magnetic_class object.
         ALLOCATE(rgrid(this%ir))
         DO r = 1, this%ir
            rgrid(r) = this%rmin + (r - 1)*(this%rmax - this%rmin)             &
     &               / (this%ir - 1)
         END DO
         ALLOCATE(zgrid(this%jz))
         DO z = 1, this%jz
            zgrid(z) = this%zmin + (z - 1)*(this%zmax - this%zmin)             &
     &               / (this%jz - 1)
         END DO

!  Create a bivariate object to do the interpolation.
!>  @todo FIXME: The biavariate interpolation is performed one plane at a time.
!>  Remove the torodial array index and make this interpolation 2D.
         CALL createbivariate(bivariate_object, numrvol, numzvol,              &
     &                        SIZE(rvolgrid, 3), 1)

!  Sum over grid points for each toroidal plane
         DO phi = 1, num_phi
            CALL setbivariate(bivariate_object, rvolgrid(:,:,phi),             &
     &                        zvolgrid(:,:,phi), rgrid, zgrid, phi, 1)

!  Interpolat from the diagnotic grid points to the volume grid points.
            CALL bivariate4pt(bivariate_object, this%a_r(:,:,phi),             &
     &                        response_rvol, phi)
            CALL bivariate4pt(bivariate_object, this%a_f(:,:,phi),             &
     &                        response_phivol, phi)
            CALL bivariate4pt(bivariate_object, this%a_z(:,:,phi),             &
     &                        response_zvol, phi)

!  J_plasma * A_diag
            sumrz = jrvolgrid(:,:,phi)*response_rvol                           &
     &            + jphivolgrid(:,:,phi)*response_phivol                       &
     &            + jzvolgrid(:,:,phi)*response_zvol
            sumphi(phi) = SUM(sumrz(2:numrvol - 1,:))                          &
     &                  + SUM(sumrz(1,:) + sumrz(numrvol,:))/2.0
         END DO

         CALL clear_bivar(bivariate_object)

!  Correct for stellarator symetry.
         magnetic_get_modeled_signal(2) = SUM(sumphi)
         IF (this%lstell_sym) THEN
            magnetic_get_modeled_signal(2)                                     &
     &         = magnetic_get_modeled_signal(2)                                & 
     &         - (sumphi(1) + sumphi(num_phi))/2.0
         END IF

!  Finish integration with the volume element.
         magnetic_get_modeled_signal(2) = magnetic_get_modeled_signal(2)       &
     &      * equilibrium_get_volume_int_element(a_model%equilibrium)

!  Deallocate work arrays.
!  NOTE: Do not deallocate any of the vol grid arrays that were obtained from
!  the equilibrium_get_magnetic_volume_*grid routines. These arrays are owned by
!  by the equilibrium.
         DEALLOCATE(response_rvol, response_phivol, response_zvol)
         DEALLOCATE(sumrz, sumphi)
         DEALLOCATE(rgrid, zgrid)

      END IF

!-------------------------------------------------------------------------------
!  Calculate the induced signal from external currents.
      extcur => equilibrium_get_ext_currents(a_model%equilibrium,              &
     &                                       scale_current)

!  If there are no external current extcur will be a null pointer
      IF (ASSOCIATED(extcur)) THEN

!  Make sure we don't go out of bounds on any array.
         num_currents = MIN(SIZE(extcur),                                      &
     &                      SIZE(this%rdiag_coilg_1))
         IF (num_currents .lt. this%n_field_cg) THEN
            CALL err_warn('magnetic_get_modeled_signal: ' //                   &
     &                    'mdsig_file expected more currents')
         END IF

         IF (scale_current) THEN
            magnetic_get_modeled_signal(3) = DOT_PRODUCT(                      &
     &         this%rdiag_coilg_1(1:num_currents),                             &
     &         extcur(1:num_currents)/this%extcur_mg)
         ELSE
            magnetic_get_modeled_signal(3) = DOT_PRODUCT(                      &
     &         this%rdiag_coilg_1(1:num_currents),                             &
     &         extcur(1:num_currents))
         END IF
      END IF

!-------------------------------------------------------------------------------

!  Add plasma signals and induced signal to create total signal.
      magnetic_get_modeled_signal(1) =                                         &
     &   SUM(magnetic_get_modeled_signal(2:3))

      sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the magnetic type.
!>
!>  Returns a description of the magnetic type for use when writting output
!>  files.
!>
!>  @param[in] this A @ref magnetic_class instance.
!>  @returns A string describing the magnetic type.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: magnetic_get_signal_type
      TYPE (magnetic_class), INTENT(in) :: this

!  Start of executable code
      magnetic_get_signal_type = 'mddc'

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref magnetic_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   magnetic_get_header
      TYPE (magnetic_class), INTENT(in) :: this

!  Start of executable code
      magnetic_get_header(1) = 'plasma only'
      magnetic_get_header(2) = 'eq currents'
      magnetic_get_header(3) = 'N/A'
      magnetic_get_header(4) = 'model_sig(1)'
      magnetic_get_header(5) = 'model_sig(2)'
      magnetic_get_header(6) = 'model_sig(3)'
      magnetic_get_header(7) = 'model_sig(4)'

      END FUNCTION

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize data to children.
!>
!>  Syncs data between the parent and child processes. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref magnetic_class instance.
!>  @param[in]    intra_comm A MPI intra_comm handle.
!-------------------------------------------------------------------------------
      SUBROUTINE magnetic_sync_children(this, intra_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                  :: intra_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                              :: error
      INTEGER, DIMENSION(3)                :: dim_lengths
      INTEGER                              :: total_length

!  Start of executable code
!  Identification Variables ----------------------------------------------------
      CALL MPI_BCAST(this%code_name, data_name_length, MPI_CHARACTER, 0,       &
     &               intra_comm, error)
      CALL MPI_BCAST(this%code_version, data_name_length, MPI_CHARACTER,       &
     &               0, intra_comm, error)
      CALL MPI_BCAST(this%date_run, data_name_length, MPI_CHARACTER, 0,        &
     &               intra_comm, error)
      CALL MPI_BCAST(this%field_coils_id, data_name_length,                    &
     &               MPI_CHARACTER, 0, intra_comm, error)

!  Coil Responce Function Variables --------------------------------------------
      CALL MPI_BCAST(this%n_field_cg, 1,  MPI_INTEGER, 0, intra_comm,          &
     &               error)
      CALL MPI_BCAST(this%rdiag_coilg_1, this%n_field_cg, MPI_REAL8, 0,        &
     &               intra_comm, error)
      CALL MPI_BCAST(this%extcur_mg, this%n_field_cg, MPI_REAL8, 0,            &
     &               intra_comm, error)

!  Plasma Response Grid Variables ----------------------------------------------
      CALL MPI_BCAST(this%kp, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%rmin, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%rmax, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%zmin, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%rmax, 1, MPI_INTEGER, 0, intra_comm, error)
      CALL MPI_BCAST(this%n_field_periods, 1, MPI_INTEGER, 0,                  &
     &               intra_comm, error)
      CALL MPI_BCAST(this%lstell_sym, 1, MPI_LOGICAL, 0, intra_comm,           &
     &               error)

!  Plasma Response Arrays ------------------------------------------------------
      CALL MPI_BCAST(dim_lengths, 3, MPI_INTEGER, 0, intra_comm, error)
      this%ir = dim_lengths(1)
      this%jz = dim_lengths(2)
      this%kp_store = dim_lengths(3)

      total_length = dim_lengths(1)*dim_lengths(2)*dim_lengths(3)
      CALL MPI_BCAST(this%a_r, total_length, MPI_REAL8, 0, intra_comm,         &
     &               error)
      CALL MPI_BCAST(this%a_f, total_length, MPI_REAL8, 0, intra_comm,         &
     &               error)
      CALL MPI_BCAST(this%a_z, total_length, MPI_REAL8, 0, intra_comm,         &
     &               error)

#endif

      END SUBROUTINE

      END MODULE
