!*******************************************************************************
!>  @file signal.f
!>  @brief Contains module @ref signal.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref signal_class and the base class of
!>  @ref combination_class.
!>  @par
!>  signal Sub Classes:
!>  @ref diagnostic, @ref geometric, @ref combination_class
!>
!>  @par
!>  combination Super Class:
!>  @ref signal
!>
!>  @note
!>  Combination is defined in the signal module to avoid a circular module
!>  dependency.
!*******************************************************************************
      MODULE signal
      USE stel_kinds, only: rprec
      USE diagnostic
      USE geometric

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
!>  Maximum length of the @ref s_name string.
      INTEGER, PARAMETER :: signal_name_short_length = 30
!>  Maximum length of the @ref l_name string.
      INTEGER, PARAMETER :: signal_name_long_length = 80
!>  Maximum length of the @ref units string.
      INTEGER, PARAMETER :: signal_units_length = 30

!>  Type descriptor for the @ref signal_class no subclass.
      INTEGER, PARAMETER :: signal_no_type = -1
!>  Type descriptor for the @ref signal_class subclass @ref diagnostic
      INTEGER, PARAMETER :: signal_diagnostic_type = 0
!>  Type descriptor for the @ref signal_class subclass @ref geometric
      INTEGER, PARAMETER :: signal_geometric_type = 1
!>  Type descriptor for the @ref signal_class subclass @ref combination_class
!  Use combination_class since combination lacks it's own module.
      INTEGER, PARAMETER :: signal_combination_type = 2

!*******************************************************************************
!  combination module parameters
!*******************************************************************************
!>  Maximum length of the combination type descritpion string.
      INTEGER, PARAMETER :: combination_type_length= 3

!>  Type descriptor for combination type no type.
      INTEGER, PARAMETER :: combination_none = -1
!>  Type descriptor for combination type sum.
      INTEGER, PARAMETER :: combination_sum = 0
!>  Type descriptor for combination type max.
      INTEGER, PARAMETER :: combination_max = 1
!>  Type descriptor for combination type min.
      INTEGER, PARAMETER :: combination_min = 2

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
!>  @ref diagnostic, @ref geometric, @ref combination_class
!-------------------------------------------------------------------------------
      TYPE signal_class
!>  Type descirptor of the signal subclass type.
!>  @par Possible values are:
!>  * @ref signal_no_type
!>  * @ref signal_diagnostic_type,
!>  * @ref signal_geometric_type
!>  * @ref signal_combination_type
         INTEGER                            :: type = signal_no_type
!>  Short name of the signal
         CHARACTER (len=signal_name_short_length) :: s_name
!>  Long name of the signal
         CHARACTER (len=signal_name_long_length)  :: l_name
         CHARACTER (len=signal_units_length)      :: units
         REAL (rprec)                       :: observed
         REAL (rprec)                       :: observed_sigma
         REAL (rprec)                       :: weight
!> Cached value of the modeled signal.
         REAL (rprec), DIMENSION(4)         :: modeled
!> Cached value of the modeled sigma.
         REAL (rprec), DIMENSION(4)         :: modeled_sigma
!>  An instance of a @ref diagnostic object.
         TYPE (diagnostic_class), POINTER   :: diagnostic => null()
!>  An instance of a @ref geometric object.
         TYPE (geometric_class), POINTER    :: geometric => null()
!>  An instance of a @ref combination object.
         TYPE (combination_class), POINTER  :: combination => null()
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

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the constriction of @ref signal_class types using
!>  @ref signal_construct_diagnostic, @ref signal_construct_geometric,
!>  @ref signal_construct_combination or @ref signal_construct_diagnostic_netcdf
!-------------------------------------------------------------------------------
      INTERFACE signal_construct
         MODULE PROCEDURE signal_construct_diagnostic,                         &
     &                    signal_construct_geometric,                          &
     &                    signal_construct_combination,                        &
     &                    signal_construct_diagnostic_netcdf
      END INTERFACE

      CONTAINS 
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref diagnostic object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
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
      TYPE (diagnostic_class), POINTER, INTENT(in) :: diagnostic_object
      CHARACTER(len=*), INTENT(in)        :: s_name
      CHARACTER(len=*), INTENT(in)        :: l_name
      CHARACTER(len=*), INTENT(in)        :: units
      REAL(rprec), INTENT(in)             :: observed
      REAL(rprec), INTENT(in)             :: sigma
      REAL(rprec), INTENT(in)             :: weight

!  Start of executable code
      ALLOCATE(signal_construct_diagnostic)

      signal_construct_diagnostic%type = signal_diagnostic_type
      signal_construct_diagnostic%s_name = s_name
      signal_construct_diagnostic%l_name = l_name
      signal_construct_diagnostic%units = units
      signal_construct_diagnostic%observed = observed
      signal_construct_diagnostic%observed_sigma = sigma
      signal_construct_diagnostic%weight = weight
      signal_construct_diagnostic%diagnostic => diagnostic_object
      signal_construct_diagnostic%modeled = 0.0
      signal_construct_diagnostic%modeled_sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref geometric object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
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
      TYPE (geometric_class), POINTER, INTENT(in) :: geometric_object
      CHARACTER(len=*), INTENT(in)       :: s_name
      CHARACTER(len=*), INTENT(in)       :: l_name
      CHARACTER(len=*), INTENT(in)       :: units
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  Start of executable code
      ALLOCATE(signal_construct_geometric)

      signal_construct_geometric%type = signal_geometric_type
      signal_construct_geometric%s_name = s_name
      signal_construct_geometric%l_name = l_name
      signal_construct_geometric%units = units
      signal_construct_geometric%observed = observed
      signal_construct_geometric%observed_sigma = sigma
      signal_construct_geometric%weight = weight
      signal_construct_geometric%geometric => geometric_object
      signal_construct_geometric%modeled = 0.0
      signal_construct_geometric%modeled_sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref combination_class
!>  object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
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
      TYPE (combination_class), POINTER, INTENT(in)                            &
     &                             :: combination_object
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight

!  Start of executable code
      ALLOCATE(signal_construct_combination)

      signal_construct_combination%type = signal_combination_type
      signal_construct_combination%s_name = s_name
      signal_construct_combination%l_name = l_name
      signal_construct_combination%units = units
      signal_construct_combination%observed = observed
      signal_construct_combination%observed_sigma = sigma
      signal_construct_combination%weight = weight
      signal_construct_combination%combination => combination_object
      signal_construct_combination%modeled = 0.0
      signal_construct_combination%modeled_sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref diagnostic object.
!>
!>  Reads the s_name, l_name and units from a netcdf file. Used for constructing
!>  the @ref magnetic objects.
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
      TYPE (diagnostic_class), POINTER, INTENT(in) :: diagnostic_object
      INTEGER, INTENT(in)                      :: mdsig_iou
      REAL(rprec), INTENT(in)                  :: observed
      REAL(rprec), INTENT(in)                  :: sigma
      REAL(rprec), INTENT(in)                  :: weight

!  local variables
      CHARACTER (len=signal_name_short_length) :: s_name
      CHARACTER (len=signal_name_long_length)  :: l_name
      CHARACTER (len=signal_units_length)      :: units

!  Start of executable code
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_s_name', s_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_l_name', l_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_units', units)

      signal_construct_diagnostic_netcdf =>                                     &
     &   signal_construct(diagnostic_object, s_name, l_name, units,             &
     &                   observed, sigma, weight)

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
      TYPE (signal_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%diagnostic)) THEN
         CALL diagnostic_destruct(this%diagnostic)
         this%diagnostic => null()
      END IF

      IF (ASSOCIATED(this%geometric)) THEN
         CALL geometric_destruct(this%geometric)
         this%geometric => null()
      END IF

      IF (ASSOCIATED(this%combination)) THEN
         CALL combination_destruct(this%combination)
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
!>  @see combination_get_modeled_signal
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

!  Start of executable code
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

            CASE (signal_combination_type)
               signal_get_modeled_signal =                                     &
     &            combination_get_modeled_signal(this%combination,             &
     &                                           a_model, sigma)

         END SELECT

         this%modeled = signal_get_modeled_signal
         this%modeled_sigma = sigma

      END IF

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

!  Local Variables

!  Start of executable code
! g^2 = W*(O - M)^2/(sigma_o^2 + sigma_m^2) = e*e
      signal_get_g2 = signal_get_e(this, a_model, use_cache)**2.0

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

!  Local Variables
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma,                &
     &                                           use_cache)

! e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
      signal_get_e = (this%observed - modeled_signal(1))                       &
     &             * SQRT(this%weight/signal_get_sigma2(this))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the total sigma^2 of a signal.
!>
!>  Gets the total sigma^2 of an individual signal.
!>  sigma^2 = sigma_o^2 + sigma_m^2
!>
!>  @param[in]    this      A @ref diagnostic_class instance.
!>  @returns The e value.
!>  @note this assumes that modeled sigma value has already been cached by a
!>  call to @ref signal_get_modeled_signal
!-------------------------------------------------------------------------------
      FUNCTION signal_get_sigma2(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                    :: signal_get_sigma2
      TYPE (signal_class), INTENT(in) :: this

!  Start of executable code
!  sigma_o^2 + sigma_m^2)
      signal_get_sigma2 = this%observed_sigma**2.0                             &
     &                  + this%modeled_sigma(1)**2.0

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
!>  @see combination_get_signal_type
!>
!>  @param[in] this A @ref signal_class instance.
!>  @returns A string describing the signal type.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: signal_get_signal_type
      TYPE (signal_class), INTENT(in)  :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_signal_type =                                           &
     &         TRIM(diagnostic_get_signal_type(this%diagnostic))

         CASE (signal_geometric_type)
            signal_get_signal_type =                                           &
     &         TRIM(geometric_get_signal_type(this%geometric))

         CASE (signal_combination_type)
            signal_get_signal_type =                                           &
     &         TRIM(combination_get_signal_type(this%combination))

      END SELECT

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
!>  @see combination_get_header
!>
!>  @param[in] this A @ref signal_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   signal_get_header
      TYPE (signal_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_header = diagnostic_get_header(this%diagnostic)

         CASE (signal_geometric_type)
            signal_get_header = geometric_get_header(this%geometric)

         CASE (signal_combination_type)
            signal_get_header = combination_get_header(this%combination)

      END SELECT

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
      CHARACTER (len=signal_name_short_length) :: signal_make_short_name
      CHARACTER (len=*), INTENT(in)            :: name
      INTEGER, INTENT(in)                      :: index

!  Start of executable code
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

!  Start of executable code
      WRITE (iou,1000) index, signal_get_signal_type(this), this%s_name,       &
     &                 signal_get_g2(this, a_model, .true.),                   &
     &                 this%weight, this%observed, this%observed_sigma,        &
     &                 this%modeled, this%modeled_sigma

1000  FORMAT(i4,2x,a23,2x,a20,12(2x,es12.5))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out the signal header information to an output file.
!>
!>  Signals information is formated as a the index, signal type, s_name, g^2,
!>  weight, observed value, observed sigma, modeled value and modeled sigma.
!>
!>  @param[in] this    A @ref diagnostic_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @note This should only be called when the type of signal changes.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write_header(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,1000) signal_get_header(this)

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
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] index   A index of a signal.
!-------------------------------------------------------------------------------
      SUBROUTINE signals_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou
      INTEGER, INTENT(in)             :: index

!  Start of executable code
      SELECT CASE (this%type)

         CASE (signal_combination_type)
            CALL combination_write_auxiliary(this%combination,                 &
     &                                       iou, index)

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
!>  @param[in] this                     A @ref diagnostic_class instance.
!>  @param[in] a_model                  The equilibrium model.
!>  @param[in] result_ncid              A NetCDF id of the result file.
!>  @param[in] current_step             The surrent reconstruction step.
!>  @param[in] index                    A index of a signal.
!>  @param[in] signal_model_value_id    NetCDF variable id of the model value.
!>  @param[in] signal_sigma_value_id    NetCDF variable id of the total sigma.
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
      REAL (rprec), DIMENSION(4)         :: modeled_sigma
      REAL (rprec), DIMENSION(4)         :: modeled_signal

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .true.)

      status = nf_put_var1_double(result_ncid, signal_model_value_id,          &
     &                            (/ index, current_step /),                   &
     &                            modeled_signal(1))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var1_double(result_ncid, signal_sigma_value_id,          &
     &                            (/ index, current_step /),                   &
     &                            SQRT(signal_get_sigma2(this)))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

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

!  Start of executable code
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
      TYPE (combination_class), POINTER, INTENT(inout) :: this

!  Local Variables
      INTEGER                                          :: i

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
      TYPE (signal_class), INTENT(in), POINTER :: signal
      INTEGER, INTENT(in)                      :: index
      REAL (rprec), INTENT(in)                 :: a

!  Local Variables

!  Start of executable code
      this%a(index) = a
      this%signals(index)%p => signal

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

!  Start of executable code
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
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: combination_get_signal_type
      TYPE (combination_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)
         CASE (combination_sum)
            combination_get_signal_type = 'coosig sum'

         CASE (combination_max)
            combination_get_signal_type = 'coosig max'

         CASE (combination_min)
            combination_get_signal_type = 'coosig min'

      END SELECT

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
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   combination_get_header
      TYPE (combination_class), INTENT(in) :: this

!  Start of executable code
      combination_get_header(1) = 'model(2)'
      combination_get_header(2) = 'model(3)'
      combination_get_header(3) = 'model(4)'
      combination_get_header(4) = 'model_sig(1)'
      combination_get_header(5) = 'model_sig(2)'
      combination_get_header(6) = 'model_sig(3)'
      combination_get_header(7) = 'model_sig(4)'

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out auxiliary signal information to an output file.
!>
!>  Writes out the s_name and coefficient of the combined signals.
!>
!>  @param[in] this    A @ref combination_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] index   A index of this signal.
!-------------------------------------------------------------------------------
      SUBROUTINE combination_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou
      INTEGER, INTENT(in)                  :: index

!  local variables
      INTEGER                              :: i

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,1000) index, combination_get_signal_type(this)
      WRITE (iou,1001)

      DO i = 1, SIZE(this%signals)
         WRITE (iou,1002) i, this%signals(i)%p%s_name, this%a(i)
      END DO

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a combination of other signals, type: ',a)
1001  FORMAT('term #',2x,'s_name',16x,'Coefficient')
1002  FORMAT(2x,i4,2x,a20,2x,es12.5)

      END SUBROUTINE

      END MODULE
