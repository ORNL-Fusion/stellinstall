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
      USE magnetic_response
#if defined (MPI_OPT)
      USE mpi
#endif
      USE profiler

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
!>  Setting to force the computing of the coil constribution.
         LOGICAL                                 :: force_coil_reponse
!>  Magnetic response function object.
         TYPE (magnetic_response_class), POINTER :: response => null()
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
         MODULE PROCEDURE magnetic_construct_netcdf,                           &
     &                    magnetic_construct_mpi
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
!>  @param[in] mdsig_iou          An instance of a the netcdf id of the open
!>                                mdsig file.
!>  @param[in] use_coil_response  Toggle if the total signal contains the
!>                                induced signal in addition to the plasma only
!>                                signal.
!>  @param[in] force_coil_reponse Force the coil response to be computed.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_netcdf(mdsig_iou, use_coil_response,         &
     &                                   force_coil_reponse)
      USE ezcdf
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: magnetic_construct_netcdf
      INTEGER, INTENT(in)            :: mdsig_iou
      LOGICAL, INTENT(in)            :: use_coil_response
      LOGICAL, INTENT(in)            :: force_coil_reponse

!  local variables
      INTEGER, DIMENSION(3)          :: dim_lengths
      LOGICAL                        :: temp_flag
      REAL (rprec)                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(magnetic_construct_netcdf)

      magnetic_construct_netcdf%response =>                                    &
     &   magnetic_response_construct(mdsig_iou)

      IF (force_coil_reponse) THEN
         CALL assert(magnetic_response_use_coil(                               &
     &                  magnetic_construct_netcdf%response),                   &
     &               'Cannot force coil response. Coil responses ' //          &
     &               'not found in mdsig file.')
      END IF
      magnetic_construct_netcdf%force_coil_reponse = force_coil_reponse

!  Check if using the coil response for this diagnostic. If not using it clear
!  the responce flag bit.
      IF (.not.use_coil_response) THEN
         CALL magnetic_response_clr_use_coil(                                  &
     &           magnetic_construct_netcdf%response)
      END IF

      CALL profiler_set_stop_time('magnetic_construct_netcdf',                 &
     &                            start_time)

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
      INTEGER                        :: total_length

!  Start of executable code
      ALLOCATE(magnetic_construct_mpi)

      CALL MPI_BCAST(magnetic_construct_mpi%force_coil_reponse, 1,             &
     &               MPI_LOGICAL, 0, intra_comm, error)
      magnetic_construct_mpi%response =>                                       &
     &   magnetic_response_construct_mpi(intra_comm)
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
      IF (ASSOCIATED(this%response)) THEN
         CALL magnetic_response_destruct(this%response)
         this%response => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  The modeled signal can be computed in one of two ways. Modeled as a full
!>  coil or as a point measurement.
!>
!>  @param[in]    this    A @ref magnetic_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal
      TYPE (magnetic_class), INTENT(in)         :: this
      TYPE (model_class), INTENT(in)            :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out)   :: sigma

!  local variables
      REAL (rprec)                              :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (magnetic_response_is_point(this%response)) THEN
         magnetic_get_modeled_signal =                                         &
     &      magnetic_get_modeled_signal_point(this, a_model, sigma)
      ELSE
         magnetic_get_modeled_signal =                                         &
     &      magnetic_get_modeled_signal_coil(this, a_model, sigma)
      END IF

      CALL profiler_set_stop_time('magnetic_get_modeled_signal',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the point modeled signal.
!>
!>  The plasma signal is the surface integral of the virtual current denisty.
!>  This is computed using the virtual casing routines in LIBSTELL.
!-------------------------------------------------------------------------------
!>
!>  External current contribrutions are included by the dot product of the
!>  external currents with a linear combination of the external fields.
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
      FUNCTION magnetic_get_modeled_signal_point(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal_point
      TYPE (magnetic_class), INTENT(in)       :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      REAL (rprec), DIMENSION(3)              :: b_field
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_get_modeled_signal_point = 0.0

      magnetic_get_modeled_signal_point(2) =                                   &
     &   DOT_PRODUCT(b_field, this%response%direction)

      CALL profiler_set_stop_time('magnetic_get_modeled_signal_point',         &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the coil modeled signal.
!>
!>  The plasma signal is the volume integral of the current denisty dotted with
!>  the diagnostic response function vector. The response_function contains a rz
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
      FUNCTION magnetic_get_modeled_signal_coil(this, a_model, sigma)
      USE model
      USE bivariate
      USE v3_utilities, only: err_warn

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal_coil
      TYPE (magnetic_class), INTENT(in)         :: this
      TYPE (model_class), INTENT(in)            :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out)   :: sigma

!  local variables
      INTEGER                                   :: phi, num_phi, r, z
      INTEGER                                   :: numrvol, numzvol
      TYPE (bivariate_type), POINTER            :: bivariate_object
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: rgrid
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: zgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: rvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: zvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jrvolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jphivolgrid
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jzvolgrid
      REAL (rprec), DIMENSION(:,:), POINTER     :: krgrid
      REAL (rprec), DIMENSION(:,:), POINTER     :: kphigrid
      REAL (rprec), DIMENSION(:,:), POINTER     :: kzgrid
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_rvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_phivol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: response_zvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sumrz
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sumphi

      LOGICAL                                   :: scale_current
      INTEGER                                   :: num_currents
      REAL (rprec), DIMENSION(:), POINTER       :: extcur

      REAL (rprec)                              :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_get_modeled_signal_coil = 0.0

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
         num_phi = this%response%num_t_store
         ALLOCATE(sumrz(numrvol, numzvol))
         ALLOCATE(sumphi(num_phi))

!  Create the r and z gridpoints for the magnetic response function.
!>  @todo FIXME: For a single diagnostic these never change. Allocate these once
!>  and cache these in the @ref magnetic_class object.
         ALLOCATE(rgrid(this%response%num_r))
         DO r = 1, this%response%num_r
            rgrid(r) = this%response%rmin                                      &
     &               + (r - 1)*(this%response%rmax - this%response%rmin)       &
     &               / (this%response%num_r - 1)
         END DO
         ALLOCATE(zgrid(this%response%num_z))
         DO z = 1, this%response%num_z
            zgrid(z) = this%response%zmin                                      &
     &               + (z - 1)*(this%response%zmax - this%response%zmin)       &
     &               / (this%response%num_z - 1)
         END DO

!  Create a bivariate object to do the interpolation.
         bivariate_object => bivariate_construct(numrvol, numzvol)

!  Sum over grid points for each toroidal plane
         DO phi = 1, num_phi
            CALL bivariate_set_grids(bivariate_object,                         &
     &                               rvolgrid(:,:,phi),                        &
     &                               zvolgrid(:,:,phi),                        &
     &                               rgrid, zgrid)

!  Interpolate from the diagnotic grid points to the volume grid points.
            CALL bivariate_get_4pt(bivariate_object,                           &
     &                             this%response%a_r(:,:,phi),                 &
     &                             response_rvol)
            CALL bivariate_get_4pt(bivariate_object,                           &
     &                             this%response%a_f(:,:,phi),                 &
     &                             response_phivol)
            CALL bivariate_get_4pt(bivariate_object,                           &
     &                             this%response%a_z(:,:,phi),                 &
     &                             response_zvol)

!  J_plasma * A_diag
!  First compute the rz grid signal. Then compute the signal for a single phi
!  plane. Over the radial component, the ends need to be taken as half.
            sumrz = jrvolgrid(:,:,phi)*response_rvol                           &
     &            + jphivolgrid(:,:,phi)*response_phivol                       &
     &            + jzvolgrid(:,:,phi)*response_zvol
            sumphi(phi) = SUM(sumrz(2:numrvol - 1,:))                          &
     &                  + SUM(sumrz(1,:) + sumrz(numrvol,:))/2.0
         END DO

         CALL bivariate_destruct(bivariate_object)

!  Correct for stellarator symetry.
         magnetic_get_modeled_signal_coil(2) = SUM(sumphi)
         IF (magnetic_response_is_stell_sym(this%response)) THEN
            magnetic_get_modeled_signal_coil(2)                                &
     &         = magnetic_get_modeled_signal_coil(2)                           &
     &         - (sumphi(1) + sumphi(num_phi))/2.0
         END IF

!  Finish integration with the volume element.
         magnetic_get_modeled_signal_coil(2)                                   &
     &      = magnetic_get_modeled_signal_coil(2)                              &
     &      * equilibrium_get_volume_int_element(a_model%equilibrium)

!  Deallocate work arrays.
!  NOTE: Do not deallocate any of the vol grid arrays that were obtained from
!  the equilibrium_get_magnetic_volume_*grid routines. These arrays are owned by
!  by the equilibrium.
         DEALLOCATE(response_rvol, response_phivol, response_zvol)
         DEALLOCATE(sumrz, sumphi)
         DEALLOCATE(rgrid, zgrid)

      END IF

!  Conducting Shell ------------------------------------------------------------
      IF (magnetic_response_use_shell(this%response)) THEN

         num_phi = this%response%num_t_shell_store

         krgrid =>                                                             &
     &      equilibrium_get_con_surface_krgrid(a_model%equilibrium)
         kphigrid =>                                                           &
     &      equilibrium_get_con_surface_kphigrid(a_model%equilibrium)
         kzgrid =>                                                             &
     &      equilibrium_get_con_surface_kzgrid(a_model%equilibrium)

         IF (ASSOCIATED(krgrid) .and. ASSOCIATED(kphigrid) .and.               &
     &       ASSOCIATED(kzgrid)) THEN
            ALLOCATE(sumphi(num_phi))

            DO phi = 1, num_phi
               sumphi(phi)                                                     &
     &            = SUM(krgrid(:,phi)*this%response%a_s_r(:,phi))              &
     &            + SUM(kphigrid(:,phi)*this%response%a_s_f(:,phi))            &
     &            + SUM(kzgrid(:,phi)*this%response%a_s_z(:,phi))
            END DO

!  Correct for stellarator symetry.
            magnetic_get_modeled_signal_coil(4) = SUM(sumphi)
            IF (magnetic_response_is_stell_sym(this%response)) THEN
               magnetic_get_modeled_signal_coil(4)                             &
     &            = magnetic_get_modeled_signal_coil(4)                        &
     &            - (sumphi(1) + sumphi(num_phi))/2.0
            END IF

!  Finish integration with the area element.
            magnetic_get_modeled_signal_coil(4)                                &
     &         = magnetic_get_modeled_signal_coil(4)                           &
     &         * equilibrium_get_area_int_element(a_model%equilibrium)

            DEALLOCATE(sumphi)
         END IF

      END IF

!  Coil Pickup -----------------------------------------------------------------
      IF (this%force_coil_reponse) THEN
         num_currents = SIZE(this%response%inductance)
      ELSE
         num_currents = -1
      END IF

!  Calculate the induced signal from external currents.
      extcur => equilibrium_get_ext_currents(a_model%equilibrium,              &
     &                                       num_currents,                     &
     &                                       scale_current)

!  If there are no external current extcur will be a null pointer
      IF (ASSOCIATED(extcur)) THEN

!  Make sure we don't go out of bounds on any array.
         num_currents = MIN(SIZE(extcur),                                      &
     &                      SIZE(this%response%inductance))
         IF (num_currents .lt. this%response%n_field_cg) THEN
            CALL err_warn('magnetic_get_modeled_signal_coil: ' //              &
     &                    'mdsig_file expected more currents')
         END IF

         IF (scale_current) THEN
            magnetic_get_modeled_signal_coil(3) = DOT_PRODUCT(                 &
     &         this%response%inductance(1:num_currents),                       &
     &         extcur(1:num_currents)/this%response%current_scale)
         ELSE
            magnetic_get_modeled_signal_coil(3) = DOT_PRODUCT(                 &
     &         this%response%inductance(1:num_currents),                       &
     &         extcur(1:num_currents))
         END IF
      END IF

!  Total Signal ----------------------------------------------------------------

!  Add plasma signals and induced signal to create total signal.
      IF (magnetic_response_use_coil(this%response)) THEN
         magnetic_get_modeled_signal_coil(1) =                                 &
     &      SUM(magnetic_get_modeled_signal_coil(2:3))
      ELSE
         magnetic_get_modeled_signal_coil(1) =                                 &
     &      magnetic_get_modeled_signal_coil(2)
      END IF

!  Add shell currents if used.
      IF (magnetic_response_use_shell(this%response)) THEN
         magnetic_get_modeled_signal_coil(1)                                   &
     &      = magnetic_get_modeled_signal_coil(1)                              &
     &      + magnetic_get_modeled_signal_coil(4)
      END IF

      sigma = 0.0

      CALL profiler_set_stop_time('magnetic_get_modeled_signal_coil',          &
     &                            start_time)

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

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_get_signal_type = 'mddc'

      CALL profiler_set_stop_time('magnetic_get_signal_type',                  &
     &                            start_time)

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
      TYPE (magnetic_class), INTENT(in)              :: this

!  local variables
      REAL (rprec)                                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      magnetic_get_header(1) = 'plasma only'
      magnetic_get_header(2) = 'eq currents'

      IF (magnetic_response_use_shell(this%response)) THEN
         magnetic_get_header(3) = 'shell currents'
      ELSE
         magnetic_get_header(3) = 'N/A'
      END IF

      magnetic_get_header(4) = 'model_sig(1)'
      magnetic_get_header(5) = 'model_sig(2)'
      magnetic_get_header(6) = 'model_sig(3)'
      magnetic_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('magnetic_get_header', start_time)

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
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_BCAST(this%force_coil_reponse, 1, MPI_LOGICAL, 0,               &
     &               intra_comm, error)
      CALL magnetic_response_sync_children(this%response, intra_comm)

      CALL profiler_set_stop_time('magnetic_sync_children', start_time)
#endif

      END SUBROUTINE

      END MODULE
