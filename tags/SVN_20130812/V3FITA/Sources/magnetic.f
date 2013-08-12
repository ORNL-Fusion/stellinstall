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
      USE mddc_T, only: mddc_mrf, mddc_mrf_destroy, mddc_construct

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
!-------------------------------------------------------------------------------
      TYPE magnetic_class
         TYPE (mddc_mrf), POINTER :: responce_function => null()
      END TYPE magnetic_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the constriction of @ref magnetic_class types using
!>  @ref magnetic_construct_netcdf or @ref magnetic_construct_responce
!-------------------------------------------------------------------------------
      INTERFACE magnetic_construct
      MODULE PROCEDURE magnetic_construct_netcdf,                              &
     &                 magnetic_construct_responce
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
!>  @param[in] responce_function An instance of a @ref mddc_mrf object.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_responce(responce_function)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER :: magnetic_construct_responce
      TYPE (mddc_mrf), POINTER, INTENT(in) :: responce_function

!  Start of executable code
      ALLOCATE(magnetic_construct_responce)

      magnetic_construct_responce%responce_function => responce_function

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref magnetic_class from a mdsig file.
!>
!>  Reads the description of a magnetic diagnostic from an mdsig file.
!>  Constructs an instance of a @ref mddc_mrf and return the constructed
!>  magnetic_class from @ref magnetic_construct_responce
!>  @param[in] mdsig_iou An instance of a the netcdf id of the open mdsig file.
!>  @returns A pointer to a constructed @ref magnetic_class object.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct_netcdf(mdsig_iou)
      USE data_parameters
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER  :: magnetic_construct_netcdf
      INTEGER, INTENT(in)             :: mdsig_iou

!  local variables
      TYPE (mddc_mrf), POINTER        :: responce_function

      CHARACTER (len=data_name_length) :: code_name
      CHARACTER (len=data_name_length) :: code_version
      CHARACTER (len=data_name_length) :: date_run
      CHARACTER (len=data_name_length) :: field_coils_id
      INTEGER                          :: n_field_cg
      REAL (rprec), DIMENSION(:), POINTER :: rdiag_coilg_1
      REAL (rprec), DIMENSION(:), POINTER :: extcur_mg
      INTEGER                          :: kp
      REAL (rprec)                     :: rmin
      REAL (rprec)                     :: rmax
      REAL (rprec)                     :: zmin
      REAL (rprec)                     :: zmax
      INTEGER                          :: n_field_periods
      LOGICAL                          :: lstell_sym
      INTEGER, DIMENSION(3)            :: dim_lengths
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_r
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_f
      REAL (rprec), DIMENSION(:,:,:), POINTER :: a_z

!  Start of executable code
      ALLOCATE(responce_function)

      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_code_name', code_name)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_code_version', code_version)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_date_run', date_run)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_field_coils_id', field_coils_id)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_n_field_cg', n_field_cg)

      ALLOCATE(rdiag_coilg_1(n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rdiag_coilg_1', rdiag_coilg_1)

      ALLOCATE(extcur_mg(n_field_cg))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_extcur_mg', extcur_mg)

      CALL cdf_read(mdsig_iou, 'mddc_mrf_kp', kp)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmin', rmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_rmax', rmax)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmin', zmin)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_zmax', zmax)
      CALL cdf_read(mdsig_iou,                                                 &
     &              'mddc_mrf_n_field_periods', n_field_periods)
      CALL cdf_read(mdsig_iou, 'mddc_mrf_lstell_sym', lstell_sym)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_r', dim_lengths)
      ALLOCATE(a_r(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_r', a_r)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_f', dim_lengths)
      ALLOCATE(a_f(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_f', a_f)

      CALL cdf_inquire(mdsig_iou, 'mddc_mrf_a_z', dim_lengths)
      ALLOCATE(a_z(dim_lengths(1),dim_lengths(2),dim_lengths(3)))
      CALL cdf_read(mdsig_iou, 'mddc_mrf_a_z', a_z)

      CALL mddc_construct(responce_function, code_name, code_version,          &
     &                    date_run, field_coils_id, rdiag_coilg_1,             &
     &                    extcur_mg, kp, rmin, rmax, zmin, zmax,               &
     &                    n_field_periods, lstell_sym, a_r, a_f, a_z)

      magnetic_construct_netcdf => magnetic_construct(responce_function)

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
      TYPE (magnetic_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%responce_function)) THEN
         CALL mddc_mrf_destroy(this%responce_function)
         this%responce_function => null()
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
!>  the diagnostic responce function vector. The responce_function contains a r z
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
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: rvolgrid => null()
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: zvolgrid => null()
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jrvolgrid => null()
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jphivolgrid => null()
      REAL (rprec), DIMENSION(:,:,:), POINTER   :: jzvolgrid => null()
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: responce_rvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: responce_phivol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: responce_zvol
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: sumrz
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sumphi

      LOGICAL                                   :: scale_current
      INTEGER                                   :: num_currents
      REAL (rprec), DIMENSION(:), POINTER       :: extcur

!  Start of executable code
      magnetic_get_modeled_signal = 0.0

!  Calculate the plasma only signal from the responce function.
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

!  Allocate space for the interpolated responce functions. Only need two
!  dimensions because the interpolation is performed one phi position at a time.
      numrvol = SIZE(rvolgrid, 1)
      numzvol = SIZE(rvolgrid, 2)
      ALLOCATE(responce_rvol(numrvol, numzvol))
      ALLOCATE(responce_phivol(numrvol, numzvol))
      ALLOCATE(responce_zvol(numrvol, numzvol))

!  Allocate work arrays to hold the results of the integration.
      num_phi = this%responce_function%kp_store
      ALLOCATE(sumrz(numrvol, numzvol))
      ALLOCATE(sumphi(num_phi))

!  Create the r and z gridpoints for the magnetic responce function.
!>  @todo FIXME: For a single diagnostic these never change. Allocate these once
!>  and cache these in the @ref magnetic_class object.
      ALLOCATE(rgrid(this%responce_function%ir))
      DO r = 1, this%responce_function%ir
         rgrid(r) = this%responce_function%rmin                                &
     &            + (r - 1)*(this%responce_function%rmax -                     &
     &                       this%responce_function%rmin)                      &
     &            / (this%responce_function%ir - 1)
      END DO
      ALLOCATE(zgrid(this%responce_function%jz))
      DO z = 1, this%responce_function%jz
         zgrid(z) = this%responce_function%zmin                                &
     &            + (z - 1)*(this%responce_function%zmax -                     &
     &                       this%responce_function%zmin)                      &
     &            / (this%responce_function%jz - 1)
      END DO

!  Create a bivariate object to do the interpolation.
!>  @todo FIXME: The biavariate interpolation is performed one plane at a time.
!>  Remove the torodial array index and make this interpolation 2D.
      CALL createbivariate(bivariate_object, numrvol, numzvol,                 &
     &                     SIZE(rvolgrid, 3), 1)

!  Sum over grid points for each toroidal plane
      DO phi = 1, num_phi
         CALL setbivariate(bivariate_object, rvolgrid(:,:,phi),                &
     &                     zvolgrid(:,:,phi), rgrid, zgrid, phi, 1)

!  Interpolat from the diagnotic grid points to the volume grid points.
         CALL bivariate4pt(bivariate_object,                                   &
     &                     this%responce_function%a_r(:,:,phi),                &
     &                     responce_rvol, phi)
         CALL bivariate4pt(bivariate_object,                                   &
     &                     this%responce_function%a_f(:,:,phi),                &
     &                     responce_phivol, phi)
         CALL bivariate4pt(bivariate_object,                                   &
     &                     this%responce_function%a_z(:,:,phi),                &
     &                     responce_zvol, phi)

!  J_plasma * A_diag
         sumrz = jrvolgrid(:,:,phi)*responce_rvol                              &
     &         + jphivolgrid(:,:,phi)*responce_phivol                          &
     &         + jzvolgrid(:,:,phi)*responce_zvol
         sumphi(phi) = SUM(sumrz(2:numrvol - 1,:))                             &
     &               + SUM(sumrz(1,:) + sumrz(numrvol,:))/2.0
      END DO

      CALL clear_bivar(bivariate_object)

!  Correct for stellarator symetry.
      magnetic_get_modeled_signal(2) = SUM(sumphi)
      IF (this%responce_function%lstell_sym) THEN
         magnetic_get_modeled_signal(2) = magnetic_get_modeled_signal(2)       &
     &      - (sumphi(1) + sumphi(num_phi))/2.0
      END IF

!  Finish integration with the volume element.
      magnetic_get_modeled_signal(2) = magnetic_get_modeled_signal(2)          &
     &   * equilibrium_get_volume_int_element(a_model%equilibrium)

!  Deallocate work arrays.
!  NOTE: Do not deallocate any of the vol grid arrays that were obtained from
!  the equilibrium_get_magnetic_volume_*grid routines. These arrays are owned by
!  by the equilibrium.
      DEALLOCATE(responce_rvol, responce_phivol, responce_zvol)
      DEALLOCATE(sumrz, sumphi)
      DEALLOCATE(rgrid, zgrid)

!-------------------------------------------------------------------------------
!  Calculate the induced signal from external currents.
      extcur => equilibrium_get_ext_currents(a_model%equilibrium,              &
     &                                       scale_current)

!  If there are no external current extcur will be a null pointer
      IF (ASSOCIATED(extcur)) THEN

!  Make sure we don't go out of bounds on any array.
         num_currents = MIN(SIZE(extcur),                                      &
     &                      SIZE(this%responce_function%rdiag_coilg_1))
         IF (num_currents .lt. this%responce_function%n_field_cg) THEN
            CALL err_warn('magnetic_get_modeled_signal: ' //                   &
     &                    'mdsig_file expected more currents')
         END IF

         IF (scale_current) THEN
            magnetic_get_modeled_signal(3) = DOT_PRODUCT(                      &
     &         this%responce_function%rdiag_coilg_1(1:num_currents),           &
     &         extcur(1:num_currents)/this%responce_function%extcur_mg)
         ELSE
            magnetic_get_modeled_signal(3) = DOT_PRODUCT(                      &
     &         this%responce_function%rdiag_coilg_1(1:num_currents),           &
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

      END MODULE
