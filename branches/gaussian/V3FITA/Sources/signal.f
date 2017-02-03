!*******************************************************************************
!>  @file signal.f
!>  @brief Contains module @ref signal.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref signal_class and the base class of
!>  @ref combination_class and @ref gaussp_class.
!>  @par
!>  signal Sub Classes:
!>  @ref diagnostic, @ref geometric, @ref prior, @ref combination_class,
!>  @ref gaussp_class
!>
!>  @par
!>  combination and gaussian process Super Class:
!>  @ref signal
!>
!>  @note
!>  Combination and Gaussian process are defined in the signal module to
!>  avoid a circular module dependency.
!*******************************************************************************
      MODULE signal
      USE stel_kinds, only: rprec
      USE diagnostic
      USE geometric
      USE prior
      USE feedback
      USE data_parameters

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
!>  Type descriptor for the @ref signal_class no subclass.
      INTEGER, PARAMETER :: signal_no_type          = -1
!>  Type descriptor for the @ref signal_class subclass @ref diagnostic
      INTEGER, PARAMETER :: signal_diagnostic_type  = 0
!>  Type descriptor for the @ref signal_class subclass @ref geometric
      INTEGER, PARAMETER :: signal_geometric_type   = 1
!>  Type descriptor for the @ref signal_class subclass @ref prior
      INTEGER, PARAMETER :: signal_prior_type       = 2
!>  Type descriptor for the @ref signal_class subclass @ref feedback
      INTEGER, PARAMETER :: signal_feedback_type    = 3
!>  Type descriptor for the @ref signal_class subclass @ref combination_class
!  Use combination_class since combination lacks it's own module.
      INTEGER, PARAMETER :: signal_combination_type = 4
!>  Type descriptor for the @ref signal_class subclass @ref gaussp_class
!  Use gaussp_class since gaussian process lacks it's own module.
      INTEGER, PARAMETER :: signal_gaussp_type = 5
!*******************************************************************************
!  combination module parameters
!*******************************************************************************
!>  Maximum length of the combination type descritpion string.
      INTEGER, PARAMETER :: combination_type_length = 3

!>  Type descriptor for combination type no type.
      INTEGER, PARAMETER :: combination_none = -1
!>  Type descriptor for combination type sum.
      INTEGER, PARAMETER :: combination_sum  = 0
!>  Type descriptor for combination type max.
      INTEGER, PARAMETER :: combination_max  = 1
!>  Type descriptor for combination type min.
      INTEGER, PARAMETER :: combination_min  = 2
!*******************************************************************************
!  gaussian process module parameters
!*******************************************************************************
!>  Maximum length of the gaussp type descritpion string.
      INTEGER, PARAMETER :: gaussp_type_length = 10

!>  Type descriptor for gaussp type no type.
      INTEGER, PARAMETER :: gaussp_no_type    =  -1
!>  Type descriptor for the gaussp type sxrem type.
      INTEGER, PARAMETER :: gaussp_sxrem_type = 0

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) signal base class
!  2) signal pointer type
!  3) combination base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a signal.
!>  @par Sub Classes:
!>  @ref diagnostic, @ref geometric, @ref prior, @ref combination_class,
!>  @ref gaussp_class
!-------------------------------------------------------------------------------
      TYPE signal_class
!>  Type descirptor of the signal subclass type.
!>  @par Possible values are:
!>  * @ref signal_no_type
!>  * @ref signal_diagnostic_type,
!>  * @ref signal_geometric_type
!>  * @ref signal_combination_type
!>  * @ref signal_gaussp_type
         INTEGER                            :: type = signal_no_type
!>  Short name of the signal
         CHARACTER (len=data_short_name_length) :: s_name
!>  Long name of the signal
         CHARACTER (len=data_name_length)       :: l_name
         CHARACTER (len=data_short_name_length) :: units
!>  Eperimentally measured signal value.
         REAL (rprec)                       :: observed
!>  Eperimentally measured signal uncertainty.
         REAL (rprec)                       :: observed_sigma
!>  Weighting parameter of the signal.
         REAL (rprec)                       :: weight
!>  Cached value of the modeled signal.
         REAL (rprec), DIMENSION(4)         :: modeled
!>  Cached value of the modeled sigma.
         REAL (rprec), DIMENSION(4)         :: modeled_sigma
!>  An instance of a @ref diagnostic object.
         TYPE (diagnostic_class), POINTER   :: diagnostic => null()
!>  An instance of a @ref geometric object.
         TYPE (geometric_class), POINTER    :: geometric => null()
!>  An instance of a @ref prior object.
         TYPE (prior_class), POINTER        :: prior => null()
!>  An instance of a @ref feedback object.
         TYPE (feedback_class), POINTER     :: feedback => null()
!>  An instance of a @ref combination object.
         TYPE (combination_class), POINTER  :: combination => null()
!>  An instance of a @ref gaussp object.
         TYPE (gaussp_class), POINTER       :: gaussp => null()
      END TYPE signal_class

!-------------------------------------------------------------------------------
!>  Pointer to a signal object. Used for creating arrays of signal pointers.
!>  This is needed because fortran does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE signal_pointer
         TYPE (signal_class), POINTER :: p => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Base class representing a combination signal.
!>  @par Super Class:
!>  @ref signal
!-------------------------------------------------------------------------------
      TYPE combination_class
!>  Type descirptor of the cobination type.
!>  @par Possible values are:
!>  * @ref combination_none
!>  * @ref combination_sum
!>  * @ref combination_max
!>  * @ref combination_min
         INTEGER :: type = combination_none
!>  Array of @ref signal_class pointer presenting the sigals to combine.
         TYPE (signal_pointer), DIMENSION(:), POINTER :: signals
!>  Array of coeffients to for the combined signals.
         REAL(rprec), DIMENSION(:), POINTER           :: a
      END TYPE combination_class
!-------------------------------------------------------------------------------
!>  Base class representing a gaussian process signal.
!>  @par Super Class:
!>  @ref signal
!-------------------------------------------------------------------------------
      TYPE gaussp_class
!> Type descriptor of the gaussian procees type
!> @par Possible values are:
!> * @ref gaussp_no_type
!> * @ref gaussp_sxrem_type
! Keep for now, but this might not be needed later
          INTEGER :: type = gaussp_no_type
!> Array of @ref diagnostic_class 
          TYPE (signal_pointer), DIMENSION(:), POINTER :: signals
!> Array of weights for the signals used in the gaussian process.
          REAL(rprec), DIMENSION(:), POINTER           :: a
! ERIC GP TODO (more detailed description)
!> Matrices used in the calculations of the gaussian process
          REAL(rprec), ALLOCATABLE, DIMENSION(:) :: matrix_sigmay
          REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: matrix_kss
          REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: matrix_kls
          REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: matrix_kll
          REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: matrix_sigmaf
!> Lower Triangular matrix found using a Cholesky Factorization
          REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: cholesky_l
!> Positive factor to add to the diagonal entries of the matrix before factoring
          REAL(rprec) :: cholesky_fac =0.0
!> Array of radial locations to evaluate the underlying function
          REAL(rprec), ALLOCATABLE, DIMENSION(:) :: r_array
!> Array of the calculated most likely model values at r_array
          REAL(rprec), ALLOCATABLE, DIMENSION(:) :: mfs_array
!> Array of observed signals
          REAL(rprec), ALLOCATABLE, DIMENSION(:) :: y_array
          INTEGER  :: n_signals = 0 
          INTEGER  :: n_points = 0
          INTEGER  :: profile_id = -1
          INTEGER  :: kernel_id = -1
!> ID of the radial grid of the profile set by GP
          INTEGER  :: as_id = -1
!> ID of the value of the profile set by GP
          INTEGER  :: af_id = -1
!> Flag to check if gaussp_get_modeled_signal has been called for this gp signal
          LOGICAL  :: first_get_call=.TRUE.

      END TYPE gaussp_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref signal_class types using
!>  @ref signal_construct_diagnostic, @ref signal_construct_geometric,
!>  @ref signal_construct_combination, @ref signal_construct_gaussp,
!>  or @ref signal_construct_diagnostic_netcdf
!-------------------------------------------------------------------------------
      INTERFACE signal_construct
         MODULE PROCEDURE signal_construct_new,                                &
     &                    signal_construct_diagnostic,                         &
     &                    signal_construct_geometric,                          &
     &                    signal_construct_prior,                              &
     &                    signal_construct_feedback,                           &
     &                    signal_construct_combination,                        &
     &                    signal_construct_gaussp,                             &
     &                    signal_construct_diagnostic_netcdf
      END INTERFACE

      CONTAINS 
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct new @ref signal_class object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] signal_type An instance of a @ref diagnostic subclass.
!>  @param[in] s_name      The short name of the signal.
!>  @param[in] l_name      The short name of the signal.
!>  @param[in] units       The units of the signal.
!>  @param[in] observed    Observed value of the signal.
!>  @param[in] sigma       Observed sigma of the signal.
!>  @param[in] weight            Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_new(signal_type, s_name, l_name, units,        &
     &                              observed, sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_new
      INTEGER                      :: signal_type
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight

!  local variables
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(signal_construct_new)

!  Warn the user that the signal is missing a non zero sigma.
      IF (sigma .eq. 0.0) THEN
         WRITE (*,1000) TRIM(s_name)
      END IF

      signal_construct_new%type = signal_type
      signal_construct_new%s_name = s_name
      signal_construct_new%l_name = l_name
      signal_construct_new%units = units
      signal_construct_new%observed = observed
      signal_construct_new%observed_sigma = sigma
      signal_construct_new%weight = weight
      signal_construct_new%modeled = 0.0
      signal_construct_new%modeled_sigma = 0.0

      CALL profiler_set_stop_time('signal_construct_new', start_time)

1000  FORMAT('Warning: sigma for ',a,' is zero.')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref diagnostic object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] diagnostic_object An instance of a @ref diagnostic subclass.
!>  @param[in] s_name            The short name of the signal.
!>  @param[in] l_name            The short name of the signal.
!>  @param[in] units             The units of the signal.
!>  @param[in] observed          Observed value of the signal.
!>  @param[in] sigma             Observed sigma of the signal.
!>  @param[in] weight            Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_diagnostic(diagnostic_object, s_name,          &
     &                                     l_name, units, observed,            &
     &                                     sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER        :: signal_construct_diagnostic
      TYPE (diagnostic_class), POINTER    :: diagnostic_object
      CHARACTER(len=*), INTENT(in)        :: s_name
      CHARACTER(len=*), INTENT(in)        :: l_name
      CHARACTER(len=*), INTENT(in)        :: units
      REAL(rprec), INTENT(in)             :: observed
      REAL(rprec), INTENT(in)             :: sigma
      REAL(rprec), INTENT(in)             :: weight

!  local variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_diagnostic =>                                           &
     &   signal_construct(signal_diagnostic_type, s_name, l_name, units,       &
     &                    observed, sigma, weight)
      signal_construct_diagnostic%diagnostic => diagnostic_object

      CALL profiler_set_stop_time('signal_construct_diagnostic',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref geometric object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] geometric_object An instance of a @ref geometric subclass.
!>  @param[in] s_name           The short name of the signal.
!>  @param[in] l_name           The short name of the signal.
!>  @param[in] units            The units of the signal.
!>  @param[in] observed         Observed value of the signal.
!>  @param[in] sigma            Observed sigma of the signal.
!>  @param[in] weight           Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_geometric(geometric_object, s_name,            &
     &                                    l_name, units, observed,             &
     &                                    sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER       :: signal_construct_geometric
      TYPE (geometric_class), POINTER    :: geometric_object
      CHARACTER(len=*), INTENT(in)       :: s_name
      CHARACTER(len=*), INTENT(in)       :: l_name
      CHARACTER(len=*), INTENT(in)       :: units
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  local variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_geometric =>                                            &
     &   signal_construct(signal_geometric_type, s_name, l_name, units,        &
     &                    observed, sigma, weight)
      signal_construct_geometric%geometric => geometric_object

      CALL profiler_set_stop_time('signal_construct_geometric',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref prior object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] prior_object An instance of a @ref prior subclass.
!>  @param[in] s_name       The short name of the signal.
!>  @param[in] l_name       The short name of the signal.
!>  @param[in] units        The units of the signal.
!>  @param[in] observed     Observed value of the signal.
!>  @param[in] sigma        Observed sigma of the signal.
!>  @param[in] weight       Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_prior(prior_object, s_name, l_name,            &
     &                                units, observed, sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER       :: signal_construct_prior
      TYPE (prior_class), POINTER        :: prior_object
      CHARACTER(len=*), INTENT(in)       :: s_name
      CHARACTER(len=*), INTENT(in)       :: l_name
      CHARACTER(len=*), INTENT(in)       :: units
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  local variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_prior =>                                                &
     &   signal_construct(signal_prior_type, s_name, l_name, units,            &
     &                    observed, sigma, weight)
      signal_construct_prior%prior => prior_object

      CALL profiler_set_stop_time('signal_construct_prior', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref prior object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] prior_object An instance of a @ref prior subclass.
!>  @param[in] s_name       The short name of the signal.
!>  @param[in] l_name       The short name of the signal.
!>  @param[in] units        The units of the signal.
!>  @param[in] observed     Observed value of the signal.
!>  @param[in] sigma        Observed sigma of the signal.
!>  @param[in] weight       Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_feedback(feedback_object, s_name,              &
     &                                   l_name, units, observed, sigma,       &
     &                                   weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER       :: signal_construct_feedback
      TYPE (feedback_class), POINTER     :: feedback_object
      CHARACTER(len=*), INTENT(in)       :: s_name
      CHARACTER(len=*), INTENT(in)       :: l_name
      CHARACTER(len=*), INTENT(in)       :: units
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  local variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_feedback =>                                             &
     &   signal_construct(signal_feedback_type, s_name, l_name, units,         &
     &                    observed, sigma, weight)
      signal_construct_feedback%feedback => feedback_object

      CALL profiler_set_stop_time('signal_construct_feedback',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref combination_class
!>  object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] combination_object An instance of a @ref combination_class
!>                                subclass.
!>  @param[in] s_name             The short name of the signal.
!>  @param[in] l_name             The short name of the signal.
!>  @param[in] units              The units of the signal.
!>  @param[in] observed           Observed value of the signal.
!>  @param[in] sigma              Observed sigma of the signal.
!>  @param[in] weight             Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_combination(combination_object, s_name,        &
     &                                      l_name, units, observed,           &
     &                                      sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_combination
      TYPE (combination_class), POINTER :: combination_object
      CHARACTER(len=*), INTENT(in)      :: s_name
      CHARACTER(len=*), INTENT(in)      :: l_name
      CHARACTER(len=*), INTENT(in)      :: units
      REAL(rprec), INTENT(in)           :: observed
      REAL(rprec), INTENT(in)           :: sigma
      REAL(rprec), INTENT(in)           :: weight

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_combination =>                                          &
     &   signal_construct(signal_combination_type, s_name, l_name,             &
     &                    units, observed, sigma, weight)
      signal_construct_combination%combination => combination_object

      CALL profiler_set_stop_time('signal_construct_combination',              &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref gaussp_class
!>  object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] combination_object An instance of a @ref combination_class
!>                                subclass.
!>  @param[in] s_name             The short name of the signal.
!>  @param[in] l_name             The short name of the signal.
!>  @param[in] units              The units of the signal.
!>  @param[in] observed           Observed value of the signal.
!>  @param[in] sigma              Observed sigma of the signal.
!>  @param[in] weight             Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_gaussp(gaussp_object, s_name,                  &
     &                                      l_name, units, observed,           &
     &                                      sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_gaussp
      TYPE (gaussp_class), POINTER :: gaussp_object
      CHARACTER(len=*), INTENT(in)      :: s_name
      CHARACTER(len=*), INTENT(in)      :: l_name
      CHARACTER(len=*), INTENT(in)      :: units
      REAL(rprec), INTENT(inout)        :: observed
      REAL(rprec), INTENT(inout)        :: sigma
      REAL(rprec), INTENT(inout)        :: weight

!  local variables
      REAL (rprec)                      :: start_time
      REAL(rprec), PARAMETER            :: tol =1.0e-8
!  Start of executable code
      start_time = profiler_get_start_time()

!  Warn the user that the observed signal for the Gaussian process is non-zero.
      IF ( ABS(observed)>tol ) THEN
         WRITE (*,1001) TRIM(s_name)
      END IF
! The Gaussian process should always have an observed value of 0.0
      observed = 0.0
! The Gaussian process should always have an sigma value of 0.0 to avoid 
! double counting. Note that the sigmas are included in the calculation of 
! the gp profile and the return gp is the log evidence which accounts for
! these errors
      IF ( ABS(sigma)>tol ) THEN
         WRITE (*,1002) TRIM(s_name)
      END IF
      sigma = 0.0

      IF (ABS(weight) > tol) THEN
         WRITE (*,1003) TRIM(s_name)
         WRITE (*,*) "The ability to calculate chi^2 from the Gaussian "
         WRITE (*,*) "process is still under development."
      END IF
      weight = 0.0

      signal_construct_gaussp =>                                               &
     &   signal_construct(signal_gaussp_type, s_name, l_name,                  &
     &                    units, observed, sigma, weight)
      signal_construct_gaussp%gaussp => gaussp_object

      CALL profiler_set_stop_time('signal_construct_combination',              &
     &                            start_time)

1001  FORMAT('Warning: Setting the observed value of ', a, ' to zero.')
1002  FORMAT('Warning: Setting the sigma value of ', a, ' to zero.')
1003  FORMAT('Warning: Setting the weight of ', a, ' to zero.')  
      END FUNCTION


!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref diagnostic object.
!>
!>  Reads the s_name, l_name and units from a netcdf file. Used for constructing
!>  the @ref magnetic objects.
!>
!>  @param[in] diagnostic_object An instance of a @ref diagnostic subclass.
!>  @param[in] mdsig_iou         NetCDF id of the opened mdsig file.
!>  @param[in] observed          Observed value of the signal.
!>  @param[in] sigma             Observed sigma of the signal.
!>  @param[in] weight            Weight factor of the signal.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_diagnostic_netcdf(diagnostic_object,           &
     &                                            mdsig_iou, observed,         &
     &                                            sigma, weight)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_diagnostic_netcdf
      TYPE (diagnostic_class), POINTER :: diagnostic_object
      INTEGER, INTENT(in)                    :: mdsig_iou
      REAL(rprec), INTENT(in)                :: observed
      REAL(rprec), INTENT(in)                :: sigma
      REAL(rprec), INTENT(in)                :: weight

!  local variables
      CHARACTER (len=data_short_name_length) :: s_name
      CHARACTER (len=data_name_length)       :: l_name
      CHARACTER (len=data_short_name_length) :: units
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL cdf_read(mdsig_iou, 'diagnostic_desc_s_name', s_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_l_name', l_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_units', units)

      signal_construct_diagnostic_netcdf =>                                    &
     &   signal_construct(diagnostic_object, s_name, l_name, units,            &
     &                    observed, sigma, weight)

      CALL profiler_set_stop_time('signal_construct_diagnostic_netcdf',        &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref signal_class object.
!>
!>  Deallocates memory and uninitializes a @ref signal_class object.
!>
!>  @param[inout] this A @ref signal_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%diagnostic)) THEN
         CALL diagnostic_destruct(this%diagnostic)
         this%diagnostic => null()
      END IF

      IF (ASSOCIATED(this%geometric)) THEN
         CALL geometric_destruct(this%geometric)
         this%geometric => null()
      END IF

      IF (ASSOCIATED(this%prior)) THEN
         CALL prior_destruct(this%prior)
         this%prior => null()
      END IF

      IF (ASSOCIATED(this%feedback)) THEN
         CALL feedback_destruct(this%feedback)
         this%feedback => null()
      END IF

      IF (ASSOCIATED(this%combination)) THEN
         CALL combination_destruct(this%combination)
         this%combination => null()
      END IF

      IF (ASSOCIATED(this%gaussp)) THEN
         CALL gaussp_destruct(this%gaussp)
         this%combination => null()
      END IF

      this%type = signal_no_type
      this%s_name = ''
      this%l_name = ''
      this%units = ''
      this%observed = 0.0
      this%observed_sigma = 0.0
      this%weight = 0.0
      this%modeled = 0.0
      this%modeled_sigma = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_modeled_signal method of the subclass instance
!>  @ref signal_class was constructed with unless the signal is read from the
!>  cached value. Read from the cached value when the signal doesn't need to be
!>  recalulated. If the sigal is to be recalculated, write these to the cache.
!>
!>  @see diagnostic::diagnostic_get_modeled_signal
!>  @see geometric::geometric_get_modeled_signal
!>  @see prior::prior_get_modeled_signal
!>  @see combination_get_modeled_signal
!>  @see gaussp_get_modeled_signal
!>
!>  @param[in]    this      A @ref diagnostic_class instance.
!>  @param[inout] a_model   A @ref model instance.
!>  @param[out]   sigma     The modeled sigma.
!>  @param[in]    use_cache If false calculate the signal and sigma and write
!>                          them to the cache. Otherwise read from the cache
!>                          without recalculating the signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_modeled_signal(this, a_model, sigma,                 &
     &                                   use_cache)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: signal_get_modeled_signal
      TYPE (signal_class), INTENT(inout)      :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma
      LOGICAL, INTENT(in)                     :: use_cache

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (use_cache) THEN
         sigma = this%modeled_sigma
         signal_get_modeled_signal = this%modeled
      ELSE

         SELECT CASE(this%type)
            CASE (signal_diagnostic_type)
               signal_get_modeled_signal =                                     &
     &            diagnostic_get_modeled_signal(this%diagnostic,               &
     &                                          a_model, sigma)

            CASE (signal_geometric_type)
               signal_get_modeled_signal =                                     &
     &            geometric_get_modeled_signal(this%geometric,                 &
     &                                         a_model, sigma)

            CASE (signal_prior_type)
               signal_get_modeled_signal =                                     &
     &            prior_get_modeled_signal(this%prior, a_model, sigma)

            CASE (signal_feedback_type)
               signal_get_modeled_signal =                                     &
     &            feedback_get_modeled_signal(this%feedback, a_model,          &
     &                                        sigma)

            CASE (signal_combination_type)
               signal_get_modeled_signal =                                     &
     &            combination_get_modeled_signal(this%combination,             &
     &                                           a_model, sigma)

            CASE (signal_gaussp_type)
               signal_get_modeled_signal =                                     &
     &            gaussp_get_modeled_signal(this%gaussp,                       &
     &                                           a_model, sigma)
         END SELECT

         this%modeled = signal_get_modeled_signal
         this%modeled_sigma = sigma

      END IF

      CALL profiler_set_stop_time('signal_get_modeled_signal',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the observed signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_modeled_signal method of the subclass instance
!>  @ref signal_class was constructed with. All signals expect feed back signals
!>  use a fixed observed value.
!>
!>  @see feedback::feedback_get_modeled_signal
!>
!>  @param[in]    this      A @ref diagnostic_class instance.
!>  @param[inout] a_model   A @ref model instance.
!>  @returns The observed value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_observed_signal(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: signal_get_observed_signal
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model

!  local variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)
         CASE (signal_feedback_type)
            signal_get_observed_signal =                                       &
     &         feedback_get_observed_signal(this%feedback, a_model)

         CASE  DEFAULT
            signal_get_observed_signal = this%observed

      END SELECT

      CALL profiler_set_stop_time('signal_get_observed_signal',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the g^2 contribution of a signal.
!>
!>  Gets the g^2 contribution of an individual signal.
!>  g^2 = W*(O - M)^2/(sigma_o^2 + sigma_m^2) = e*e
!>
!>  @param[in]    this      A @ref diagnostic_class instance.
!>  @param[inout] a_model   A @ref model instance.
!>  @param[in]    use_cache If the signals have been already calculated this
!>                          should be false. Otherwise set to true.
!>  @returns The g^2 value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_g2(this, a_model, use_cache)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                       :: signal_get_g2
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model
      LOGICAL, INTENT(in)                :: use_cache

!  local Variables
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  g^2 = W*(O - M)^2/(sigma_o^2 + sigma_m^2) = e*e
      signal_get_g2 = signal_get_e(this, a_model, use_cache)**2.0

      CALL profiler_set_stop_time('signal_get_g2', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the e contribution of a signal.
!>
!>  Gets the e contribution of an individual signal.
!>  e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
!>
!>  @param[in]    this      A @ref diagnostic_class instance.
!>  @param[inout] a_model   A @ref model instance.
!>  @param[in]    use_cache If the signals have been already calculated this
!>                          should be false. Otherwise set to true.
!>  @returns The e value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_e(this, a_model, use_cache)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                       :: signal_get_e
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model
      LOGICAL, INTENT(in)                :: use_cache

!  local Variables
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma
      REAL (rprec)                       :: observed_signal
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      observed_signal = signal_get_observed_signal(this, a_model)

      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma,                &
     &                                           use_cache)

!  e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
!  Users are can be lax when specifiying the sigma for signals weighted to zero.
!  Avoid a possible divide by zero by error checking the weight first. If the
!  weight is zero, set the error vector explicitly to zero.
      IF (this%weight .eq. 0.0) THEN
         signal_get_e = 0.0
      ELSE
         signal_get_e = (observed_signal - modeled_signal(1))                  &
     &                * SQRT(this%weight/signal_get_sigma2(this))
      END IF

      CALL profiler_set_stop_time('signal_get_e', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the total sigma^2 of a signal.
!>
!>  Gets the total sigma^2 of an individual signal.
!>  sigma^2 = sigma_o^2 + sigma_m^2
!>
!>  @param[in] this A @ref diagnostic_class instance.
!>  @returns The e value.
!>  @note this assumes that modeled sigma value has already been cached by a
!>  call to @ref signal_get_modeled_signal
!-------------------------------------------------------------------------------
      FUNCTION signal_get_sigma2(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                    :: signal_get_sigma2
      TYPE (signal_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  sigma_o^2 + sigma_m^2)
      signal_get_sigma2 = this%observed_sigma**2.0                             &
     &                  + this%modeled_sigma(1)**2.0

      CALL profiler_set_stop_time('signal_get_sigma2', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the signal type.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_signal_type method of the subclass instance
!>  @ref signal_class was constructed with.
!>
!>  @see diagnostic::diagnostic_get_signal_type
!>  @see geometric::geometric_get_signal_type
!>  @see prior::prior_get_signal_type
!>  @see combination_get_signal_type
!>  @see gaussp_get_signal_type
!>
!>  @param[in] this A @ref signal_class instance.
!>  @returns A string describing the signal type.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: signal_get_signal_type
      TYPE (signal_class), INTENT(in)  :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_signal_type =                                           &
     &         TRIM(diagnostic_get_signal_type(this%diagnostic))

         CASE (signal_geometric_type)
            signal_get_signal_type =                                           &
     &         TRIM(geometric_get_signal_type(this%geometric))

         CASE (signal_prior_type)
            signal_get_signal_type =                                           &
     &         TRIM(prior_get_signal_type(this%prior))

         CASE (signal_feedback_type)
            signal_get_signal_type =                                           &
     &         TRIM(feedback_get_signal_type(this%feedback))

         CASE (signal_combination_type)
            signal_get_signal_type =                                           &
     &         TRIM(combination_get_signal_type(this%combination))

         CASE (signal_gaussp_type)
            signal_get_signal_type =                                           &
     &         TRIM(gaussp_get_signal_type(this%gaussp))

      END SELECT

      CALL profiler_set_stop_time('signal_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_header method of the subclass instance
!>  @ref diagnostic_class was constructed with.
!>
!>  @see diagnostic::diagnostic_get_header
!>  @see geometric::geometric_get_header
!>  @see prior::prior_get_header
!>  @see combination_get_header
!>  @see gaussp_get_header
!>
!>  @param[in] this A @ref signal_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_header(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   signal_get_header
      TYPE (signal_class), INTENT(in)                :: this

!  local variables
      REAL (rprec)                                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_header = diagnostic_get_header(this%diagnostic)

         CASE (signal_geometric_type)
            signal_get_header = geometric_get_header(this%geometric)

         CASE (signal_prior_type)
            signal_get_header = prior_get_header(this%prior)

         CASE (signal_feedback_type)
            signal_get_header = feedback_get_header(this%feedback)

         CASE (signal_combination_type)
            signal_get_header = combination_get_header(this%combination)

         CASE (signal_gaussp_type)
            signal_get_header = gaussp_get_header(this%gaussp)

      END SELECT

      CALL profiler_set_stop_time('signal_get_header', start_time)

      END FUNCTION
 
!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Generate a short name by appending an index.
!>
!>  Returns an appended string with no leading spaces. This check the number of
!>  digits in the index then chooses the appropriate formatting.
!>
!>  @param[in] name  A name of a signal to append to.
!>  @param[in] index A index of a signal to append.
!>  @returns A string an appended string.
!-------------------------------------------------------------------------------
      FUNCTION signal_make_short_name(name, index)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_short_name_length) :: signal_make_short_name
      CHARACTER (len=*), INTENT(in)          :: name
      INTEGER, INTENT(in)                    :: index

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Find the number of digits in the index
      IF (index .lt. 10) THEN
         WRITE(signal_make_short_name,'(a,i1)') TRIM(name), index
      ELSE IF (index .lt. 100) THEN
         WRITE(signal_make_short_name,'(a,i2)') TRIM(name), index
      ELSE IF (index .lt. 1000) THEN
         WRITE(signal_make_short_name,'(a,i3)') TRIM(name), index
      ELSE
         WRITE(signal_make_short_name,'(a,i4)') TRIM(name), index
      END IF

      CALL profiler_set_stop_time('signal_make_short_name', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write out the signal information to an output file.
!>
!>  Signals information is formated as the index, signal type, s_name, g^2,
!>  weight, observed value, observed sigma, modeled value and modeled sigma.
!>
!>  @param[in] this    A @ref diagnostic_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] index   A index of a signal.
!>  @param[in] a_model The equilibrium model.
!>  @note This assumes that signals have already been calculated. The model
!>  isn't used by signal_get_g2 since it is reading from the cache but must be
!>  passed in anyway.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write(this, iou, index, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                :: iou
      INTEGER, INTENT(in)                :: index
      TYPE (model_class), INTENT(inout)  :: a_model

!  local variables
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .true.)

      WRITE (iou,1000) index, signal_get_signal_type(this), this%s_name,       &
     &                 signal_get_g2(this, a_model, .true.),                   &
     &                 this%weight, signal_get_observed_signal(this,           &
     &                                                         a_model),       &
     &                 this%observed_sigma, modeled_signal,                    &
     &                 modeled_sigma

      CALL profiler_set_stop_time('signal_write', start_time)

1000  FORMAT(i4,2x,a23,2x,a20,12(2x,es12.5))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out the signal header information to an output file.
!>
!>  Signals information is formated as a the index, signal type, s_name, g^2,
!>  weight, observed value, observed sigma, modeled value and modeled sigma.
!>
!>  @param[in] this A @ref diagnostic_class instance.
!>  @param[in] iou  A input/output representing the file to write to.
!>  @note This should only be called when the type of signal changes.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write_header(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou

!  local variables
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (iou,*)
      WRITE (iou,1000) signal_get_header(this)

      CALL profiler_set_stop_time('signal_write_header', start_time)

1000  FORMAT(3x,'#',2x,'type',21x,'s_name',16x,'g^2',11x,'weight',8x,          &
     &       'observed',6x,'sigma', 9x,'model',7x,2(2x,a12),7(2x,a12))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out any auxiliary signal information to an output file.
!>
!>  This method is meant to be overwritten by a subclass method if need be. As a
!>  result, returns the write_auxiliary method of the subclass instance
!>  @ref signal_class was constructed with. Otherwise uses the default. nothing
!>  is written.
!>
!>  @see combination_write_auxiliary
!>  @see gaussp_write_auxiliary
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] index   A index of a signal.
!>  @param[in] a_model The equilibrium model.
!-------------------------------------------------------------------------------
      SUBROUTINE signals_write_auxiliary(this, iou, index, a_model)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou
      INTEGER, INTENT(in)             :: index
      TYPE (model_class), INTENT(in)  :: a_model

!  local variables
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)

         CASE (signal_prior_type)
            CALL prior_write_auxiliary(this%prior, iou, index,                 &
     &                                 a_model)

         CASE (signal_feedback_type)
            CALL feedback_write_auxiliary(this%feedback, iou, index)

         CASE (signal_combination_type)
            CALL combination_write_auxiliary(this%combination,                 &
     &                                       iou, index)

         CASE (signal_gaussp_type)
            CALL gaussp_write_auxiliary(this%gaussp,                           &
     &                                       iou, index)


      CALL profiler_set_stop_time('signals_write_auxiliary', start_time)

      END SELECT

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the signal data for a step to the result netcdf file.
!>
!>  Writes out the observed, modeled and total sigma value of the signal.
!>
!>  @param[in] this                  A @ref diagnostic_class instance.
!>  @param[in] a_model               The equilibrium model.
!>  @param[in] result_ncid           A NetCDF id of the result file.
!>  @param[in] current_step          The surrent reconstruction step.
!>  @param[in] index                 A index of a signal.
!>  @param[in] signal_model_value_id NetCDF variable id of the model value.
!>  @param[in] signal_sigma_value_id NetCDF variable id of the total sigma.
!>  @note This assumes that signals have already been calculated. The model
!>  isn't used by signal_get_g2 since it is reading from the cache but must be
!>  passed in anyway.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write_step_data(this, a_model, result_ncid,            &
     &                                  current_step, index,                   &
     &                                  signal_model_value_id,                 &
     &                                  signal_sigma_value_id)
      USE model
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model
      INTEGER, INTENT(in)                :: result_ncid
      INTEGER, INTENT(in)                :: current_step
      INTEGER, INTENT(in)                :: index
      INTEGER, INTENT(in)                :: signal_model_value_id
      INTEGER, INTENT(in)                :: signal_sigma_value_id

!  local variables
      INTEGER                            :: status
      INTEGER                            :: varid
      REAL (rprec), DIMENSION(4)         :: modeled_sigma
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .true.)

      status = nf_put_vara_double(result_ncid, signal_model_value_id,          &
     &                            (/ 1, index, current_step /),                &
     &                            (/ 4, 1, 1 /), modeled_signal)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var1_double(result_ncid, signal_sigma_value_id,          &
     &                            (/ index, current_step /),                   &
     &                            SQRT(signal_get_sigma2(this)))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

!  Feedback signals change the observed signal. Write the current observed
!  value at eash step.
      SELECT CASE(this%type)

         CASE (signal_feedback_type)
            status = nf_inq_varid(result_ncid, 'signal_observed_value',        &
     &                            varid)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

            status =                                                           &
     &         nf_put_var1_double(result_ncid, varid, index,                   &
     &                            signal_get_observed_signal(this,             &
     &                                                       a_model))
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

      END SELECT

      CALL profiler_set_stop_time('signals_write_auxiliary', start_time)

      END SUBROUTINE

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!  Combination Signals
!
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref combination_class.
!>
!>  Allocates memory and initializes a @ref combination_class object. This just
!>  allocates the @ref combination_class::signals and @ref combination_class::a
!>  array. The indices of these arrays are set by @ref combination_set_signal
!>
!>  @param[in] n_signals        Number of signals in this combination.
!>  @param[in] combination_type The method to combine the signals. Possible
!>                              value are @ref combination_sum,
!>                              @ref combination_max and @ ref combination_min
!>  @returns A pointer to a constructed @ref combination_class object.
!-------------------------------------------------------------------------------
      FUNCTION combination_construct(n_signals, combination_type)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), POINTER :: combination_construct
      INTEGER, INTENT(in)               :: n_signals
      CHARACTER (len=combination_type_length), INTENT(in)                      &
     &                                  :: combination_type

!  Local Variables
      INTEGER                           :: status
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(combination_construct)

      ALLOCATE(combination_construct%signals(n_signals),                       &
     &         combination_construct%a(n_signals), STAT=status)

      SELECT CASE (combination_type)
         CASE ('sum')
            combination_construct%type = combination_sum

         CASE ('max')
            combination_construct%type = combination_max

         CASE ('min')
            combination_construct%type = combination_min

      END SELECT

      CALL profiler_set_stop_time('combination_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref combination_class object.
!>
!>  Deallocates memory and uninitializes a @ref combination_class object.
!>
!>  @param[inout] this A @ref combination_class instance.
!>  @note The destructors of the signals in the signal array are not called her
!>  because the combination signals do not own the memory for those signals.
!-------------------------------------------------------------------------------
      SUBROUTINE combination_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), POINTER :: this

!  Local Variables
      INTEGER                           :: i

!  Start of executable code

!  Null all pointers in the signals array. Do not deallocate the signals in the
!  array since this object is not the owner.
      DO i = 1, SIZE(this%signals)
         this%signals(i)%p => null()
      END DO

      DEALLOCATE(this%signals)
      this%type = combination_none

      DEALLOCATE(this)

      END SUBROUTINE
!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Set the object and coefficient for an index.
!>
!>  @param[inout] this   A @ref combination_class instance.
!>  @param[in]    signal A @ref signal_class instance.
!>  @param[in]    a      A coefficient for the signal
!>  @param[in]    index  The index to place the combined @ref signal.
!>  @note Only signals that were created before this combination can be added.
!>  Otherwise there it will creat an infinitely recursive loop in
!>  @ref combination_get_modeled_signal
!-------------------------------------------------------------------------------
      SUBROUTINE combination_set_signal(this, signal, a, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(inout)  :: this
      TYPE (signal_class), POINTER             :: signal
      INTEGER, INTENT(in)                      :: index
      REAL (rprec), INTENT(in)                 :: a

!  local Variables
      REAL (rprec)                             :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      this%a(index) = a
      this%signals(index)%p => signal

      CALL profiler_set_stop_time('combination_set_signal', start_time)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  This method combines the signals pointed in the
!>  @ref combination_class::signals array. 
!>
!>  @param[in]    this    A @ref combination_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION combination_get_modeled_signal(this, a_model, sigma)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)      :: combination_get_modeled_signal
      TYPE (combination_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Local Variables
      INTEGER                                 :: i
      REAL (rprec), DIMENSION(4)              :: sigma_local
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      combination_get_modeled_signal =                                         &
     &   this%a(1)*signal_get_modeled_signal(this%signals(1)%p, a_model,       &
     &                                       sigma, .true.)

      SELECT CASE(this%type)
         CASE (combination_sum)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal =                                &
     &            combination_get_modeled_signal +                             &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.)
               sigma = sigma + sigma_local
            END DO

         CASE (combination_max)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal = MAX(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.))
               sigma = MAX(sigma, sigma_local)
            END DO

         CASE (combination_min)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal = MIN(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.))
               sigma = MIN(sigma, sigma_local)
            END DO

      END SELECT

      CALL profiler_set_stop_time('combination_get_modeled_signal',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the combination type.
!>
!>  Returns a description of the combination type for use when writting output
!>  files.
!>
!>  @param[in] this A @ref combination_class instance.
!>  @returns A string describing the magnetic type.
!-------------------------------------------------------------------------------
      FUNCTION combination_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: combination_get_signal_type
      TYPE (combination_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)
         CASE (combination_sum)
            combination_get_signal_type = 'coosig sum'

         CASE (combination_max)
            combination_get_signal_type = 'coosig max'

         CASE (combination_min)
            combination_get_signal_type = 'coosig min'

      END SELECT

      CALL profiler_set_stop_time('combination_get_signal_type',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref combination_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION combination_get_header(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   combination_get_header
      TYPE (combination_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      combination_get_header(1) = 'model(2)'
      combination_get_header(2) = 'model(3)'
      combination_get_header(3) = 'model(4)'
      combination_get_header(4) = 'model_sig(1)'
      combination_get_header(5) = 'model_sig(2)'
      combination_get_header(6) = 'model_sig(3)'
      combination_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('combination_get_header', start_time)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out auxiliary signal information to an output file.
!>
!>  Writes out the s_name and coefficient of the combined signals.
!>
!>  @param[in] this  A @ref combination_class instance.
!>  @param[in] iou   A input/output representing the file to write to.
!>  @param[in] index A index of this signal.
!-------------------------------------------------------------------------------
      SUBROUTINE combination_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou
      INTEGER, INTENT(in)                  :: index

!  local variables
      INTEGER                              :: i
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (iou,*)
      WRITE (iou,1000) index, combination_get_signal_type(this)
      WRITE (iou,1001)

      DO i = 1, SIZE(this%signals)
         WRITE (iou,1002) i, this%signals(i)%p%s_name, this%a(i)
      END DO

      CALL profiler_set_stop_time('combination_write_auxiliary',               &
     &                            start_time)

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a combination of other signals, type: ',a)
1001  FORMAT('term #',2x,'s_name',16x,'Coefficient')
1002  FORMAT(2x,i4,2x,a20,2x,es12.5)

      END SUBROUTINE


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!  Gaussian Process Signals
!
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref gaussp_class.
!>
!>  Allocates memory and initializes a @ref gaussp_class object. This just
!>  allocates the @ref gaussp_class::signals, @ref gaussp_class::a
!>  array, and matrixices for the GP computataion.  The indices of these arrays
!>  are set by @ref gaussp_set_signal, and the matrices are calculated in
!>  @ref gaussp_get_modeled_signal
!>
!>  @param[in] n_signals        Number of signals in this gaussian process.
!>  @param[in] n_points         Number of points to evaluate the model proifle
!>  @param[in] kernel_id        Index of the kernel used for this GP Signal
!>  @param[in] gaussp_type      The type of profile modeled by the gaussian process
!>  @param[in] profile_id       Index of the profile if needed (sxrem(i) for example)
!>  @returns A pointer to a constructed @ref gaussp_class object.
!-------------------------------------------------------------------------------
      FUNCTION gaussp_construct(n_signals, n_points, kernel_id,           &
     &                          gaussp_type, profile_id, cholesky_fac)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class), POINTER :: gaussp_construct
      INTEGER, INTENT(in)               :: n_signals
      INTEGER, INTENT(in)               :: n_points      
      INTEGER, INTENT(in)               :: kernel_id
      INTEGER, INTENT(in)               :: profile_id
      CHARACTER (len=gaussp_type_length), INTENT(in)                      &
     &                                  :: gaussp_type
      REAL (rprec)                      :: cholesky_fac

!  Local Variables
      INTEGER                           :: status
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(gaussp_construct)

      ALLOCATE(gaussp_construct%signals(n_signals),                       &
     &         gaussp_construct%a(n_signals), STAT=status)

      ALLOCATE (gaussp_construct%matrix_sigmay(n_signals),                &
     &         gaussp_construct%matrix_kss(n_points,n_points),            &
     &         gaussp_construct%matrix_sigmaf(n_points,n_points),         &
     &         gaussp_construct%matrix_kls(n_signals,n_points),           &
     &         gaussp_construct%matrix_kll(n_signals,n_signals),          &
     &         gaussp_construct%cholesky_l(n_signals,n_signals),          &
     &         STAT=status)

      ALLOCATE( gaussp_construct%r_array(n_points),STAT=status)
      ALLOCATE( gaussp_construct%mfs_array(n_points),STAT=status)
      ALLOCATE( gaussp_construct%y_array(n_signals),STAT=status)
      gaussp_construct%n_signals = n_signals
      gaussp_construct%n_points = n_points
      gaussp_construct%profile_id = profile_id
      gaussp_construct%kernel_id = kernel_id
      gaussp_construct%first_get_call = .TRUE.
      gaussp_construct%cholesky_fac=cholesky_fac
      SELECT CASE (gaussp_type)
         CASE ('sxrem')
            gaussp_construct%type = gaussp_sxrem_type

      END SELECT

      CALL profiler_set_stop_time('gaussp_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref gaussp_class object.
!>
!>  Deallocates memory and uninitializes a @ref gaussp_class object.
!>
!>  @param[inout] this A @ref gaussp_class instance.
!>  @note The destructors of the signals in the signal array are not called her
!>  because the gaussian process signals do not own the memory for those signals.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class), POINTER :: this

!  Local Variables
      INTEGER                           :: i

!  Start of executable code

!  Null all pointers in the signals array. Do not deallocate the signals in the
!  array since this object is not the owner.
      DO i = 1, SIZE(this%signals)
         this%signals(i)%p => null()
      END DO

      DEALLOCATE(this%signals)
      DEALLOCATE(this%a)
      this%type = gaussp_no_type

      DEALLOCATE(this%matrix_sigmay)
      DEALLOCATE(this%matrix_sigmaf)
      DEALLOCATE(this%matrix_kss)
      DEALLOCATE(this%matrix_kls)
      DEALLOCATE(this%matrix_kll)
      DEALLOCATE(this%cholesky_l)

      DEALLOCATE(this%r_array)
      DEALLOCATE(this%mfs_array)
      DEALLOCATE(this%y_array)

      this%n_signals = 0
      this%n_points = 0
      this%profile_id = -1
      this%kernel_id = -1
      this%first_get_call = .TRUE.

      DEALLOCATE(this)

      END SUBROUTINE
!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Set the object and coefficient for an index.
!>
!>  @param[inout] this   A @ref gaussp_class instance.
!>  @param[in]    signal A @ref signal_class instance.
!>  @param[in]    a      A weight for the signal
!>  @param[in]    index  The index to place the gaussp @ref signal.
!>  @note Only signals that were created before this gaussian process
!>  can be added. Otherwise it will creat an infinitely recursive loop in
!>  @ref gaussp_get_modeled_signal
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_set_signal(this, signal, a, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class), INTENT(inout)  :: this
      TYPE (signal_class), POINTER             :: signal
      INTEGER, INTENT(in)                      :: index
      REAL (rprec), INTENT(in)                 :: a

!  local Variables
      REAL (rprec)                             :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      this%a(index) = a
      this%signals(index)%p => signal
      this%y_array(index) = signal%observed

      CALL profiler_set_stop_time('gaussp_set_signal', start_time)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  This method uses a gaussian process to set a model using the
!>  the signals pointed in the  @ref gaussp_class::signals array. It then 
!>  returns the evidence. 
!>
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION gaussp_get_modeled_signal(this, a_model, sigma)
      USE model
      USE stel_constants, ONLY: twopi
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)      :: gaussp_get_modeled_signal
      TYPE (gaussp_class), INTENT(inout)        :: this
      TYPE (model_class), INTENT(inout)         :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out)   :: sigma

!  Local Variables
      INTEGER                                  :: i, j, k
      REAL (rprec), DIMENSION(4)               :: sigma_local
      REAL (rprec)                             :: start_time
      INTEGER                                  :: ierr
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work_matrix
      REAL(rprec)                              :: log_evidence
!  Start of executable code
      start_time = profiler_get_start_time()


      ALLOCATE(work_matrix(this%n_signals,this%n_points+1))


! generate matrices
! The following matricies only need to be generated the first time
! gaussp_get_modeled_signal is called for this gpsig.

! Note that there can be multiple GP signals so I cannot use a local
! variable with the SAVE attribute. Also I want to avoid duplicating
! information in the input file to help cut down on user input errors. To
! do this I use some information stored in a_model to construct this matrices.
! This is why the matrices are not generated in the construct and set subroutines

      IF (this%first_get_call) THEN
        
         SELECT CASE (this%type)
           CASE (gaussp_sxrem_type)
             this%as_id = model_get_param_id(a_model,'pp_sxrem_as_a')
             this%af_id = model_get_param_id(a_model,'pp_sxrem_af_a')
         END SELECT

          CALL gaussp_set_rarray(this,a_model)
! r_array looks okay
          CALL gaussp_generate_sigmay(this)
! sigma _y looks okay
         CALL gaussp_generate_kss(this, a_model)
! kss looks okay
         this%first_get_call = .FALSE.
      END IF

      CALL gaussp_generate_kls(this, a_model)
! kls has the correct dimesnion and does not change between here
! and the end of the subroutine

      CALL gaussp_generate_kll(this, a_model)
! kls has the correct dimesnion and does not change between here
! and the end of the subroutine


      this%cholesky_l = this%matrix_kll
! This copy does what I want it to do

      DO i = 1, this%n_signals
         this%cholesky_l(i,i) = this%cholesky_l(i,i) +                  &
     &                           this%matrix_sigmay(i) +                &
     &                           this%cholesky_fac
         work_matrix(i,1) = this%y_array(i)
         DO j = 1 ,this%n_points
            work_matrix(i,j+1) = this%matrix_kls(i,j)
         END DO
      END DO

! THE DO LOOP looks to corecly modify the cholesky matrix
! The work matrix looks to be properly initialized



! CALL LAPACK to do Cholesky factorization of the A=(K_LL + Sigma_y)
      CALL DPOTRF ('L', this%n_signals , this%cholesky_l ,              &
     &             this%n_signals, ierr)


!Cholesky l looks resonable
!l does not change after call to DROTS
!l(i,j) is not importnat of j>i
!ie setting l(i,j) =0 for j>i does not change evidence

      IF (ierr /= 0) THEN
        IF (ierr < 0) THEN
          CALL err_fatal('gaussp_get_modeled_signal: DROTRF cannot '//  & 
     &                   'factor matrix. ierr is negative' )
        ELSE IF (ierr > 0) THEN
           CALL err_fatal('gaussp_get_modeled_signal: DROTRF cannot '//  & 
     &                     'factor matrix. ierr is positive. Try'   //   &
     &                     ' increasing gp_cholesky_fact' )
         ENDIF
      END IF

      
! IF statement correctly throws an error if I cannot factor the matrix!

! CALL LAPACK to solve A^-1y and A^-1 K_LS simultaneously
! The work stores the RHS on input, and returns the solution matrix
      CALL DPOTRS ('L', this%n_signals, this%n_points+1,                 &
     &             this%cholesky_l, this%n_signals, work_matrix,         &
     &             this%n_signals, ierr)

! INPUT TO DPOTRS looks ok

      IF (ierr < 0) THEN
          CALL err_fatal('gaussp_get_modeled_signal: DROTRS failed ' //  & 
     &                   'to solve the equation' )
      END IF
  
! This error state correctly catches LAPACK ERRORS
            


! Calculate mfs, sigma_f, and evidence
! mfs = K_SL * A^-1 y
! Log evidence = -1/2 y^TA^-1 y - sum log L_ii - n/2 log 2pi 

      log_evidence = -0.5 * this%n_signals*log(twopi)
      this%mfs_array = 0.0
      DO i = 1, this%n_signals
         log_evidence = log_evidence - log(this%cholesky_l(i,i)) -      &
     &   0.5*this%y_array(i) * work_matrix(i,1)
         DO j = 1, this%n_points
           this%mfs_array(j) = this%mfs_array(j) +                      &
     &                         this%matrix_kls(i,j) * work_matrix(i,1)
         ENDDO
      END DO
      
! Sigma_f = K_ss - k_sl * A^-1 K_ls
      DO i = 1, this%n_points
        DO j = 1, this%n_points
          this%matrix_sigmaf(i,j) = this%matrix_kss(i,j)
          DO k =1, this%n_signals
            this%matrix_sigmaf(i,j) = this%matrix_sigmaf(i,j) -         &
     &      this%matrix_kls(k,i)* work_matrix(k, j+1)
          END DO
        END DO
      END DO

      CALL gaussp_set_model(this, a_model)
 
      !gaussp_get_modeled_signal(1) = SQRT(-2.0* log_evidence)
      gaussp_get_modeled_signal(1) = 0
      sigma = 1.0
      DEALLOCATE(work_matrix)

      CALL profiler_set_stop_time('gaussp_get_modeled_signal',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the gaussian process type.
!>
!>  Returns a description of the gaussian process for use when writting output
!>  files.
!>
!>  @param[in] this A @ref gaussp_class instance.
!>  @returns A string describing the gaussian process type.
!-------------------------------------------------------------------------------
      FUNCTION gaussp_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: gaussp_get_signal_type
      TYPE (gaussp_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)
         CASE (gaussp_no_type)
            gaussp_get_signal_type = 'gp none'

         CASE (gaussp_sxrem_type)
            gaussp_get_signal_type = 'gp sxrem'
 
      END SELECT

      CALL profiler_set_stop_time('gaussp_get_signal_type',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref gaussp_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION gaussp_get_header(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   gaussp_get_header
      TYPE (gaussp_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gaussp_get_header(1) = 'model(2)'
      gaussp_get_header(2) = 'model(3)'
      gaussp_get_header(3) = 'model(4)'
      gaussp_get_header(4) = 'model_sig(1)'
      gaussp_get_header(5) = 'model_sig(2)'
      gaussp_get_header(6) = 'model_sig(3)'
      gaussp_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('gaussp_get_header', start_time)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out auxiliary signal information to an output file.
!>
!>  Writes out the s_name and coefficient of the gaussp signals.
!>
!>  @param[in] this  A @ref gaussp_class instance.
!>  @param[in] iou   A input/output representing the file to write to.
!>  @param[in] index A index of this signal.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou
      INTEGER, INTENT(in)                  :: index

!  local variables
      INTEGER                              :: i
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (iou,*)
      WRITE (iou,1000) index, gaussp_get_signal_type(this)
      WRITE (iou,1001)

      DO i = 1, SIZE(this%signals)
         WRITE (iou,1002) i, this%signals(i)%p%s_name, this%a(i)
      END DO

      CALL profiler_set_stop_time('gaussp_write_auxiliary',                 &
     &                            start_time)

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a gaussian process of using signals of type: ',a)
1001  FORMAT('term #',2x,'s_name',16x,'Coefficient')
1002  FORMAT(2x,i4,2x,a20,2x,es12.5)

      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
! ERIC GP TODO I need to define the following as private
!-------------------------------------------------------------------------------
!>  @brief Subroutine sets the r_array for the gaussp class object
!>
!>  Generates r_array for the gaussp class object. The subroutine
!>  gets the radial locations for the desired profile. It assume that the 
!>  pprofile is either a spline of line segment type. There are no checks
!>  
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[in]    a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_set_rarray(this, a_model)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)       :: a_model

!  Local Variables
      INTEGER                                 :: index
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

! The value of the index in call to model_get_param value differes depending 
! on wheather or not the profile is 1D or two 2D
!> @todo FIXME: I'm using the value of profile ID TO guess at the
!>  number of indices used to define the profile. This is not fool proof

      IF (this%profile_id > 0) THEN
        DO index = 1, this%n_points
           this%r_array(index) = model_get_param_value(a_model,            &
     &                         this%as_id, this%profile_id, index)
        ENDDO
      ELSE
        DO index = 1, this%n_points
           this%r_array(index) = model_get_param_value(a_model,            &
     &                         this%as_id, index, 0)
        ENDDO
      END IF

      CALL profiler_set_stop_time('gaussp_set_rarray',                     &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Subroutine sets model coefficents from the gaussp class object
!>
!>  Set the model to the most likely coefficents. It assume that the 
!>  pprofile is either a spline of line segment type. There are no checks
!>  
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[in]    a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_set_model(this, a_model)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp_class),INTENT(in) :: this
      TYPE (model_class), INTENT(inout)       :: a_model

!  Local Variables
      INTEGER                                 :: index
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!> @todo FIXME: THE LAST argument in the function call is an mpi communicator
!> for VMEC equilibria. I don't know what it should be here.

! The value of the index in the model_set_param value differes depending 
! on wheather or not the profile is 1D or two 2D
!> @todo FIXME: I'm using the value of profile ID TO guess at the
!>  number of indices used to define the profile. This is not fool proof

      IF (this%profile_id > 0) THEN
       DO index = 1, this%n_points
         CALL model_set_param(a_model, this%af_id, this%profile_id,      &
     &                             index, this%mfs_array(index), 0)
       ENDDO
      ELSE
       DO index = 1, this%n_points
        CALL model_set_param(a_model, this%af_id, index,                &
     &                             0, this%mfs_array(index), 0)
       ENDDO
      END IF
      CALL profiler_set_stop_time('gaussp_set_model',                   &
     &                            start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Subroutine to generate the sigmay matrix
!>
!>  Generates the convariance matrix for the measured signals. Initially
!>  we assume that the measured signals are not correlated so this matrix is
!>  a diagonal matix. Thus is can be defined by a 1 dimensional array.
!>  
!>
!>  @param[in]    this    A @ref gaussp_class instance.

!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_generate_sigmay(this)

      IMPLICIT NONE

!  Delcare Arguments
      TYPE (gaussp_class) :: this

!  local variables
      INTEGER                              :: index
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      DO index = 1, this%n_signals
        this%matrix_sigmay(index) =                                     &
     &                         this%signals(index)%p%observed_sigma**2
!     &                         this%signals(index)%p%observed_sigma
      END DO

      CALL profiler_set_stop_time('gaussp_generate_sigmay',                     &
     &                            start_time)

       END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Subroutine to generate the kss kernal matrix
!>
!>  Generates the kss matrix for the gaussian process. The kss matrix
!>  evaluates at kernal function at the points x_i and x_j and stores the 
!>  result in kss(i,j). Here x_i, x_j are in the set of points where we 
!>  want to sample the model function.
!>
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[in]    a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_generate_kss(this, a_model)
      USE model

      IMPLICIT NONE

!  Delcare Arguments
      TYPE (gaussp_class) :: this
      TYPE (model_class), INTENT(inout)       :: a_model
!  local variables
      INTEGER                              :: index, jndex
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

! kss is often symmetric. If Gauss is intensive I can speed this step up
! using this fact.
      DO index = 1, this%n_points
        DO jndex = 1 , this%n_points
          this%matrix_kss(index,jndex) = model_get_gpker(a_model,        &
     &                    this%r_array(index), this%r_array(jndex),      &
     &                    this%kernel_id)
        END DO
      END DO

      CALL profiler_set_stop_time('gaussp_generate_kss',start_time)  

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Subroutine to generate the kls kernal matrix
!>
!>  Generates the kss matrix for the gaussian process. The kss matrix
!>  evaluates a kernal function where the linear operate is applied x_i
!>  and elvautes x_j. The result is stored as  kls(i,j). 
!>  Here (Lx)_i is represnts the i-th signal, x_j is the set of points where we 
!>  want to sample the model function.
!>
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[in]    a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_generate_kls(this, a_model)
      USE model

      IMPLICIT NONE

!  Delcare Arguments
      TYPE (gaussp_class) :: this
      TYPE (model_class),TARGET, INTENT(inout)  :: a_model
!  local variables
      INTEGER                              :: index, jndex
      CHARACTER(len=1), ALLOCATABLE        :: context(:)
      INTEGER                              :: context_length
      TYPE (gpsxrem_context)               :: gpsxr_context
      REAL (rprec)                         :: start_time

      REAL(rprec), PARAMETER                 :: dxinput = 0.005
!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the gp soft x-ray context
      gpsxr_context%model => a_model
      gpsxr_context%kernel_id = this%kernel_id


! kls 
! ech temp increase step size to speed up integration
      DO jndex = 1 , this%n_points
         gpsxr_context%spoint = this%r_array(jndex)
         context_length = SIZE(TRANSFER(gpsxr_context, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(gpsxr_context, context)
         DO index = 1, this%n_signals
          
          this%matrix_kls(index,jndex) =  path_integrate(               &
     &    this%signals(index)%p%diagnostic%sxrem%chord_path,            &
     &                    gpsxr_kls_function, context, dxinput)         &
!     &                    gpsxr_kls_function, context)                  &
     &                    * this%signals(index)%p%diagnostic%sxrem%geo

        END DO
        DEALLOCATE(context)
      END DO

      CALL profiler_set_stop_time('gaussp_generate_kls',start_time)  

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Subroutine to generate the kll kernal matrix
!>
!>  Generates the kss matrix for the gaussian process. The kss matrix
!>  evaluates a kernal function where the linear operate is applied x_i
!>  and elvautes x_j. The result is stored as  kls(i,j). 
!>  Here (Lx)_i is represnts the i-th signal, (Lx)_j is j-th signal
!>
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[in]    a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp_generate_kll(this, a_model)
      USE model

      IMPLICIT NONE

!  Delcare Arguments
      TYPE (gaussp_class) :: this
      TYPE (model_class), TARGET, INTENT(inout) :: a_model
!  local variables
      INTEGER                              :: index, jndex
      CHARACTER(len=1), ALLOCATABLE        :: context(:)
      INTEGER                              :: context_length
      TYPE (gpsxrem_context)               :: gpsxr_context
      REAL (rprec)                         :: start_time

      REAL(rprec), PARAMETER                 :: dxinput = 0.005
!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the gp soft x-ray context
      gpsxr_context%model => a_model
      gpsxr_context%kernel_id = this%kernel_id


! kll
      DO jndex = 1 , this%n_signals
         gpsxr_context%geo = this%signals(jndex)%p%diagnostic%sxrem%geo
         gpsxr_context%chord_path =>                                    &
     &         this%signals(jndex)%p%diagnostic%sxrem%chord_path

         context_length = SIZE(TRANSFER(gpsxr_context, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(gpsxr_context, context)

!$OMP PARALLEL DO
!$OMP& SCHEDULE(DYNAMIC)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(index)

! ech temp increase step size to speed up integration
         DO index = jndex, this%n_signals
        
          this%matrix_kll(index,jndex) =  path_integrate(                &
     &        this%signals(index)%p%diagnostic%sxrem%chord_path,         &
!     &                    gpsxr_kll1_function, context)                  &
     &                    gpsxr_kll1_function, context, dxinput)          &
     &                  * this%signals(index)%p%diagnostic%sxrem%geo
          this%matrix_kll(jndex,index) = this%matrix_kll(index,jndex) 
        END DO
!$OMP END PARALLEL DO

        DEALLOCATE(context)
      END DO

      CALL profiler_set_stop_time('gaussp_generate_kll',start_time)  

      END SUBROUTINE






      END MODULE
