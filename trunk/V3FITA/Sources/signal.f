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
      USE magnetic
      USE sxrem
      USE intpol
      USE thomson
      USE extcurz
      USE mse
      USE ece
      USE geometric
      USE prior_gaussian
      USE feedback
      USE data_parameters

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
!>  Type descriptor for the @ref signal_class no subclass.
      INTEGER, PARAMETER :: signal_no_type             = -1
!>  Type descriptor for the @ref diagnostic_class subclass @ref magnetic
      INTEGER, PARAMETER :: signal_magnetic_type       = 0
!>  Type descriptor for the @ref diagnostic_class subclass @ref sxrem
      INTEGER, PARAMETER :: signal_sxrem_type          = 1
!>  Type descriptor for the @ref diagnostic_class subclass @ref intpol
      INTEGER, PARAMETER :: signal_intpol_type         = 2
!>  Type descriptor for the @ref diagnostic_class subclass @ref thomson
      INTEGER, PARAMETER :: signal_thomson_type        = 3
!>  Type descriptor for the @ref diagnostic_class subclass @ref extcurz
      INTEGER, PARAMETER :: signal_extcurz_type        = 4
!>  Type descriptor for the @ref diagnostic_class subclass @ref mse
      INTEGER, PARAMETER :: signal_mse_type            = 5
!>  Type descriptor for the @ref diagnostic_class subclass @ref ece
      INTEGER, PARAMETER :: signal_ece_type            = 6
!>  Type descriptor for the @ref signal_class subclass @ref geometric
      INTEGER, PARAMETER :: signal_geometric_type      = 7
!>  Type descriptor for the @ref signal_class subclass @ref prior
      INTEGER, PARAMETER :: signal_prior_gaussian_type = 8
!>  Type descriptor for the @ref signal_class subclass @ref feedback
      INTEGER, PARAMETER :: signal_feedback_type       = 9
!>  Type descriptor for the @ref signal_class subclass @ref combination_class
!  Use combination_class since combination lacks it's own module.
      INTEGER, PARAMETER :: signal_combination_type    = 10

!*******************************************************************************
!  combination module parameters
!*******************************************************************************
!>  Maximum length of the combination type descritpion string.
      INTEGER, PARAMETER :: combination_type_length = 4

!>  Type descriptor for combination type no type.
      INTEGER, PARAMETER :: combination_none = -1
!>  Type descriptor for combination type sum.
      INTEGER, PARAMETER :: combination_sum  = 0
!>  Type descriptor for combination type max.
      INTEGER, PARAMETER :: combination_max  = 1
!>  Type descriptor for combination type min.
      INTEGER, PARAMETER :: combination_min  = 2
!>  Type descriptor for combination type weighted average
      INTEGER, PARAMETER :: combination_wavg = 3

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
!>  * @ref signal_magnetic_type,
!>  * @ref signal_sxrem_type
!>  * @ref signal_intpol_type,
!>  * @ref signal_thomson_type
!>  * @ref signal_extcurz_type
!>  * @ref signal_mse_type
!>  * @ref signal_ece_type
!>  * @ref signal_geometric_type
!>  * @ref signal_prior_gaussian_type
!>  * @ref signal_feedback_type
!>  * @ref signal_combination_type
         INTEGER                             :: type = signal_no_type
!>  Short name of the signal.
         CHARACTER (len=data_short_name_length) :: s_name
!>  Long name of the signal.
         CHARACTER (len=data_name_length)       :: l_name
!>  Physical units the signal measures.
         CHARACTER (len=data_short_name_length) :: units
!>  Eperimentally measured signal value.
         REAL (rprec)                        :: observed
!>  Eperimentally measured signal uncertainty.
         REAL (rprec)                        :: observed_sigma
!>  Weighting parameter of the signal.
         REAL (rprec)                        :: weight
!>  Cached value of the modeled signal.
         REAL (rprec), DIMENSION(4)          :: modeled
!>  Cached value of the modeled sigma.
         REAL (rprec), DIMENSION(4)          :: modeled_sigma
!>  Scale factor index.
         INTEGER                             :: scale_index
!>  An instance of a @ref magnetic object.
         TYPE (magnetic_class), POINTER      :: magnetic => null()
!>  An instance of a @ref sxrem object.
         TYPE (sxrem_class), POINTER         :: sxrem => null()
!>  An instance of a @ref intpol object.
         TYPE (intpol_class), POINTER        :: intpol => null()
!>  An instance of a @ref thomson object.
         TYPE (thomson_class), POINTER       :: thomson => null()
!>  An instance of a @ref extcurz object.
         TYPE (extcurz_class), POINTER       :: extcurz => null()
!>  An instance of a @ref mse object.
         TYPE (mse_class), POINTER           :: mse => null()
!>  An instance of a @ref mse object.
         TYPE (ece_class), POINTER           :: ece => null()
!>  An instance of a @ref geometric object.
         TYPE (geometric_class), POINTER     :: geometric => null()
!>  An instance of a @ref prior object.
        TYPE (prior_gaussian_class), POINTER :: gaussian => null()
!>  An instance of a @ref feedback object.
         TYPE (feedback_class), POINTER      :: feedback => null()
!>  An instance of a @ref combination object.
         TYPE (combination_class), POINTER   :: combination => null()
      END TYPE signal_class

!-------------------------------------------------------------------------------
!>  Pointer to a signal object. Used for creating arrays of signal pointers.
!>  This is needed because fortran does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE signal_pointer
!>  Pointer to a @ref signal_class. Used for building arrays of
!>  @ref signal_class objects.
         TYPE (signal_class), POINTER :: p => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the guassian process
!>  callback function.
!-------------------------------------------------------------------------------
      TYPE signal_gp_context
!>  First instance of an @ref signal_class object.
         TYPE (signal_class), POINTER :: signal1 => null()
!>  Second instance of an @ref signal_class object.
         TYPE (signal_class), POINTER :: signal2 => null()
!>  An instance of an @ref model_class object.
         TYPE (model_class), POINTER  :: model => null()
!>  Gaussian process kernel flags.
         INTEGER                      :: flags = model_state_all_off
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
!>  Array of @ref signal_class pointer presenting the signals to combine.
         TYPE (signal_pointer), DIMENSION(:), POINTER :: signals
!>  Array of coeffients to for the combined signals.
         REAL(rprec), DIMENSION(:), POINTER           :: a
!>  Index of the average weighting factor.
         INTEGER                                      :: wgt_index
      END TYPE combination_class

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
     &                    signal_construct_magnetic,                           &
     &                    signal_construct_sxrem,                              &
     &                    signal_construct_intpol,                             &
     &                    signal_construct_thomson,                            &
     &                    signal_construct_extcurz,                            &
     &                    signal_construct_mse,                                &
     &                    signal_construct_ece,                                &
     &                    signal_construct_geometric,                          &
     &                    signal_construct_prior_gaussian,                     &
     &                    signal_construct_feedback,                           &
     &                    signal_construct_combination,                        &
     &                    signal_construct_diagnostic_netcdf
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to get the guassian process kernel values.
!-------------------------------------------------------------------------------
      INTERFACE signal_get_gp
         MODULE PROCEDURE signal_get_gp_i,                                     &
     &                    signal_get_gp_s,                                     &
     &                    signal_get_gp_x
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
!>  @param[in] weight      Weight factor of the signal.
!>  @param[in] index       Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_new(signal_type, s_name, l_name, units,        &
     &                              observed, sigma, weight, index)

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
      INTEGER, INTENT(in)          :: index

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
      signal_construct_new%scale_index = index

      CALL profiler_set_stop_time('signal_construct_new', start_time)

1000  FORMAT('Warning: sigma for ',a,' is zero.')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref magnetic object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] magnetic_object An instance of a @ref magnetic subclass.
!>  @param[in] s_name          The short name of the signal.
!>  @param[in] l_name          The short name of the signal.
!>  @param[in] units           The units of the signal.
!>  @param[in] observed        Observed value of the signal.
!>  @param[in] sigma           Observed sigma of the signal.
!>  @param[in] weight          Weight factor of the signal.
!>  @param[in] index           Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_magnetic(magnetic_object, s_name,              &
     &                                   l_name, units, observed,              &
     &                                   sigma, weight, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER   :: signal_construct_magnetic
      TYPE (magnetic_class), POINTER :: magnetic_object
      CHARACTER(len=*), INTENT(in)   :: s_name
      CHARACTER(len=*), INTENT(in)   :: l_name
      CHARACTER(len=*), INTENT(in)   :: units
      REAL(rprec), INTENT(in)        :: observed
      REAL(rprec), INTENT(in)        :: sigma
      REAL(rprec), INTENT(in)        :: weight
      INTEGER, INTENT(in)            :: index

!  local variables
      REAL (rprec)                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_magnetic =>                                             &
     &   signal_construct(signal_magnetic_type, s_name, l_name, units,         &
     &                    observed, sigma, weight, index)
      signal_construct_magnetic%magnetic => magnetic_object

      CALL profiler_set_stop_time('signal_construct_magnetic',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref sxrem object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] sxrem_object An instance of a @ref sxrem subclass.
!>  @param[in] s_name       The short name of the signal.
!>  @param[in] l_name       The short name of the signal.
!>  @param[in] units        The units of the signal.
!>  @param[in] observed     Observed value of the signal.
!>  @param[in] sigma        Observed sigma of the signal.
!>  @param[in] weight       Weight factor of the signal.
!>  @param[in] index        Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_sxrem(sxrem_object, s_name, l_name,            &
     &                                units, observed, sigma, weight,          &
     &                                index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_sxrem
      TYPE (sxrem_class), POINTER  :: sxrem_object
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight
      INTEGER, INTENT(in)          :: index

!  local variables
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_sxrem =>                                                &
     &   signal_construct(signal_sxrem_type, s_name, l_name, units,            &
     &                    observed, sigma, weight, index)
      signal_construct_sxrem%sxrem => sxrem_object

      CALL profiler_set_stop_time('signal_construct_sxrem', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref intpol object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] intpol_object An instance of a @ref intpol subclass.
!>  @param[in] s_name        The short name of the signal.
!>  @param[in] l_name        The short name of the signal.
!>  @param[in] units         The units of the signal.
!>  @param[in] observed      Observed value of the signal.
!>  @param[in] sigma         Observed sigma of the signal.
!>  @param[in] weight        Weight factor of the signal.
!>  @param[in] index         Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_intpol(intpol_object, s_name, l_name,          &
     &                                 units, observed, sigma, weight,         &
     &                                 index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_intpol
      TYPE (intpol_class), POINTER :: intpol_object
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight
      INTEGER, INTENT(in)          :: index

!  local variables
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_intpol =>                                               &
     &   signal_construct(signal_intpol_type, s_name, l_name, units,           &
     &                    observed, sigma, weight, index)
      signal_construct_intpol%intpol => intpol_object

      CALL profiler_set_stop_time('signal_construct_intpol', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref thomson object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] thomson_object An instance of a @ref thomson subclass.
!>  @param[in] s_name         The short name of the signal.
!>  @param[in] l_name         The short name of the signal.
!>  @param[in] units          The units of the signal.
!>  @param[in] observed       Observed value of the signal.
!>  @param[in] sigma          Observed sigma of the signal.
!>  @param[in] weight         Weight factor of the signal.
!>  @param[in] index          Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_thomson(thomson_object, s_name, l_name,        &
     &                                  units, observed, sigma, weight,        &
     &                                  index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER  :: signal_construct_thomson
      TYPE (thomson_class), POINTER :: thomson_object
      CHARACTER(len=*), INTENT(in)  :: s_name
      CHARACTER(len=*), INTENT(in)  :: l_name
      CHARACTER(len=*), INTENT(in)  :: units
      REAL(rprec), INTENT(in)       :: observed
      REAL(rprec), INTENT(in)       :: sigma
      REAL(rprec), INTENT(in)       :: weight
      INTEGER, INTENT(in)           :: index

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_thomson =>                                              &
     &   signal_construct(signal_thomson_type, s_name, l_name, units,          &
     &                    observed, sigma, weight, index)
      signal_construct_thomson%thomson => thomson_object

      CALL profiler_set_stop_time('signal_construct_thomson',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref extcurz object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] extcurz_object An instance of a @ref extcurz subclass.
!>  @param[in] s_name         The short name of the signal.
!>  @param[in] l_name         The short name of the signal.
!>  @param[in] units          The units of the signal.
!>  @param[in] observed       Observed value of the signal.
!>  @param[in] sigma          Observed sigma of the signal.
!>  @param[in] weight         Weight factor of the signal.
!>  @param[in] index          Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_extcurz(extcurz_object, s_name, l_name,        &
     &                                  units, observed, sigma, weight,        &
     &                                  index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER  :: signal_construct_extcurz
      TYPE (extcurz_class), POINTER :: extcurz_object
      CHARACTER(len=*), INTENT(in)  :: s_name
      CHARACTER(len=*), INTENT(in)  :: l_name
      CHARACTER(len=*), INTENT(in)  :: units
      REAL(rprec), INTENT(in)       :: observed
      REAL(rprec), INTENT(in)       :: sigma
      REAL(rprec), INTENT(in)       :: weight
      INTEGER, INTENT(in)           :: index

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_extcurz =>                                              &
     &   signal_construct(signal_extcurz_type, s_name, l_name, units,          &
     &                    observed, sigma, weight, index)
      signal_construct_extcurz%extcurz => extcurz_object

      CALL profiler_set_stop_time('signal_construct_extcurz',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref mse object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] mse_object An instance of a @ref mse subclass.
!>  @param[in] s_name     The short name of the signal.
!>  @param[in] l_name     The short name of the signal.
!>  @param[in] units      The units of the signal.
!>  @param[in] observed   Observed value of the signal.
!>  @param[in] sigma      Observed sigma of the signal.
!>  @param[in] weight     Weight factor of the signal.
!>  @param[in] index      Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_mse(mse_object, s_name, l_name, units,         &
     &                              observed, sigma, weight, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_mse
      TYPE (mse_class), POINTER    :: mse_object
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight
      INTEGER, INTENT(in)          :: index

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_mse =>                                                  &
     &   signal_construct(signal_mse_type, s_name, l_name, units,              &
     &                    observed, sigma, weight, index)
      signal_construct_mse%mse => mse_object

      CALL profiler_set_stop_time('signal_construct_mse', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref ece object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] ece_object An instance of a @ref ece subclass.
!>  @param[in] s_name     The short name of the signal.
!>  @param[in] l_name     The short name of the signal.
!>  @param[in] units      The units of the signal.
!>  @param[in] observed   Observed value of the signal.
!>  @param[in] sigma      Observed sigma of the signal.
!>  @param[in] weight     Weight factor of the signal.
!>  @param[in] index      Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_ece(ece_object, s_name, l_name, units,         &
     &                              observed, sigma, weight, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_ece
      TYPE (ece_class), POINTER    :: ece_object
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight
      INTEGER, INTENT(in)          :: index

!  local variables
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_ece =>                                                  &
     &   signal_construct(signal_ece_type, s_name, l_name, units,              &
     &                    observed, sigma, weight, index)
      signal_construct_ece%ece => ece_object

      CALL profiler_set_stop_time('signal_construct_ece', start_time)

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
!>  @param[in] index            Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_geometric(geometric_object, s_name,            &
     &                                    l_name, units, observed,             &
     &                                    sigma, weight, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER    :: signal_construct_geometric
      TYPE (geometric_class), POINTER :: geometric_object
      CHARACTER(len=*), INTENT(in)    :: s_name
      CHARACTER(len=*), INTENT(in)    :: l_name
      CHARACTER(len=*), INTENT(in)    :: units
      REAL(rprec), INTENT(in)         :: observed
      REAL(rprec), INTENT(in)         :: sigma
      REAL(rprec), INTENT(in)         :: weight
      INTEGER, INTENT(in)             :: index

!  local variables
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_geometric =>                                            &
     &   signal_construct(signal_geometric_type, s_name, l_name, units,        &
     &                    observed, sigma, weight, index)
      signal_construct_geometric%geometric => geometric_object

      CALL profiler_set_stop_time('signal_construct_geometric',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref prior_gaussian
!>  object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] prior_object An instance of a @ref prior_gaussian subclass.
!>  @param[in] s_name       The short name of the signal.
!>  @param[in] l_name       The short name of the signal.
!>  @param[in] units        The units of the signal.
!>  @param[in] observed     Observed value of the signal.
!>  @param[in] sigma        Observed sigma of the signal.
!>  @param[in] weight       Weight factor of the signal.
!>  @param[in] index        Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_prior_gaussian(prior_gaussian_object,          &
     &                                         s_name, l_name, units,          &
     &                                         observed, sigma, weight,        &
     &                                         index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_prior_gaussian
      TYPE (prior_gaussian_class), POINTER  :: prior_gaussian_object
      CHARACTER(len=*), INTENT(in)          :: s_name
      CHARACTER(len=*), INTENT(in)          :: l_name
      CHARACTER(len=*), INTENT(in)          :: units
      REAL(rprec), INTENT(in)               :: observed
      REAL(rprec), INTENT(in)               :: sigma
      REAL(rprec), INTENT(in)               :: weight
      INTEGER, INTENT(in)                   :: index

!  local variables
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_prior_gaussian =>                                       &
     &   signal_construct(signal_prior_gaussian_type, s_name, l_name,          &
     &                    units, observed, sigma, weight, index)
      signal_construct_prior_gaussian%gaussian => prior_gaussian_object

      CALL profiler_set_stop_time('signal_construct_prior_gaussian',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref signal_class containing a @ref feedback object.
!>
!>  Allocates memory and initializes a @ref signal_class object.
!>
!>  @param[in] feedback_object An instance of a @ref feedback subclass.
!>  @param[in] s_name          The short name of the signal.
!>  @param[in] l_name          The short name of the signal.
!>  @param[in] units           The units of the signal.
!>  @param[in] observed        Observed value of the signal.
!>  @param[in] sigma           Observed sigma of the signal.
!>  @param[in] weight          Weight factor of the signal.
!>  @param[in] index           Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_feedback(feedback_object, s_name,              &
     &                                   l_name, units, observed, sigma,       &
     &                                   weight, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER   :: signal_construct_feedback
      TYPE (feedback_class), POINTER :: feedback_object
      CHARACTER(len=*), INTENT(in)   :: s_name
      CHARACTER(len=*), INTENT(in)   :: l_name
      CHARACTER(len=*), INTENT(in)   :: units
      REAL(rprec), INTENT(in)        :: observed
      REAL(rprec), INTENT(in)        :: sigma
      REAL(rprec), INTENT(in)        :: weight
      INTEGER, INTENT(in)            :: index

!  local variables
      REAL (rprec)                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_feedback =>                                             &
     &   signal_construct(signal_feedback_type, s_name, l_name, units,         &
     &                    observed, sigma, weight, index)
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
!>  @param[in] index              Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_combination(combination_object, s_name,        &
     &                                      l_name, units, observed,           &
     &                                      sigma, weight, index)

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
      INTEGER, INTENT(in)               :: index

!  local variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      signal_construct_combination =>                                          &
     &   signal_construct(signal_combination_type, s_name, l_name,             &
     &                    units, observed, sigma, weight, index)
      signal_construct_combination%combination => combination_object

      CALL profiler_set_stop_time('signal_construct_combination',              &
     &                            start_time)

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
!>  @param[in] index             Index of the model signal scale factor.
!>  @returns A pointer to a constructed @ref signal_class object.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_diagnostic_netcdf(magnetic_object,             &
     &                                            mdsig_iou, observed,         &
     &                                            sigma, weight, index)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_diagnostic_netcdf
      TYPE (magnetic_class), POINTER :: magnetic_object
      INTEGER, INTENT(in)                    :: mdsig_iou
      REAL(rprec), INTENT(in)                :: observed
      REAL(rprec), INTENT(in)                :: sigma
      REAL(rprec), INTENT(in)                :: weight
      INTEGER, INTENT(in)                    :: index

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
     &   signal_construct(magnetic_object, s_name, l_name, units,              &
     &                    observed, sigma, weight, index)

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
      IF (ASSOCIATED(this%magnetic)) THEN
         CALL magnetic_destruct(this%magnetic)
         this%magnetic => null()
      END IF

      IF (ASSOCIATED(this%sxrem)) THEN
         CALL sxrem_destruct(this%sxrem)
         this%sxrem => null()
      END IF

      IF (ASSOCIATED(this%intpol)) THEN
         CALL intpol_destruct(this%intpol)
         this%intpol => null()
      END IF

      IF (ASSOCIATED(this%thomson)) THEN
         CALL thomson_destruct(this%thomson)
         this%thomson => null()
      END IF

      IF (ASSOCIATED(this%extcurz)) THEN
         CALL extcurz_destruct(this%extcurz)
         this%extcurz => null()
      END IF

      IF (ASSOCIATED(this%mse)) THEN
         CALL mse_destruct(this%mse)
         this%mse => null()
      END IF

      IF (ASSOCIATED(this%ece)) THEN
         CALL ece_destruct(this%ece)
         this%ece => null()
      END IF

      IF (ASSOCIATED(this%geometric)) THEN
         CALL geometric_destruct(this%geometric)
         this%geometric => null()
      END IF

      IF (ASSOCIATED(this%gaussian)) THEN
         CALL prior_gaussian_destruct(this%gaussian)
         this%gaussian => null()
      END IF

      IF (ASSOCIATED(this%feedback)) THEN
         CALL feedback_destruct(this%feedback)
         this%feedback => null()
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
!>  recalulated. If the signal is to be recalculated, write these to the cache.
!>
!>  @see magnetic::magnetic_get_signal_type
!>  @see sxrem::sxrem_get_signal_type
!>  @see intpol::intpol_get_signal_type
!>  @see thomson::thomson_get_signal_type
!>  @see extcurz::extcurz_get_signal_type
!>  @see mse::mse_get_signal_type
!>  @see ece::ece_get_signal_type
!>  @see geometric::geometric_get_modeled_signal
!>  @see prior_guassian::prior_guassian_get_modeled_signal
!>  @see feedback::feedback_get_modeled_signal
!>  @see combination_get_modeled_signal
!>
!>  @param[inout] this      A @ref diagnostic_class instance.
!>  @param[inout] a_model   A @ref model instance.
!>  @param[out]   sigma     The modeled sigma.
!>  @param[in]    use_cache If false calculate the signal and sigma and write them
!>                          them to the cache. Otherwise read from the cache
!>                          without recalculating the signal.
!>  @param[in]  last_value  Last good value in case the signal did not change.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_modeled_signal(this, a_model, sigma,                 &
     &                                   use_cache, last_value)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: signal_get_modeled_signal
      TYPE (signal_class), INTENT(inout)      :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma
      LOGICAL, INTENT(in)                     :: use_cache
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value

!  local variables
      REAL (rprec)                            :: scale_factor
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (use_cache) THEN
         sigma = this%modeled_sigma
         signal_get_modeled_signal = this%modeled
      ELSE
         scale_factor = model_get_signal_factor(a_model,                       &
     &                                          this%scale_index)

         SELECT CASE(this%type)

            CASE (signal_magnetic_type)
               signal_get_modeled_signal =                                     &
     &            magnetic_get_modeled_signal(this%magnetic, a_model,          &
     &                                        sigma, last_value,               &
     &                                        scale_factor)

            CASE (signal_sxrem_type)
               signal_get_modeled_signal =                                     &
     &            sxrem_get_modeled_signal(this%sxrem, a_model, sigma,         &
     &                                     last_value, scale_factor)

            CASE (signal_intpol_type)
               signal_get_modeled_signal =                                     &
     &            intpol_get_modeled_signal(this%intpol, a_model, sigma,       &
     &                                      last_value, scale_factor)

            CASE (signal_thomson_type)
               signal_get_modeled_signal =                                     &
     &            thomson_get_modeled_signal(this%thomson, a_model,            &
     &                                       sigma, last_value,                &
     &                                       scale_factor)

            CASE (signal_extcurz_type)
               signal_get_modeled_signal =                                     &
     &            extcurz_get_modeled_signal(this%extcurz, a_model,            &
     &                                       sigma, last_value,                &
     &                                       scale_factor)

            CASE (signal_mse_type)
               signal_get_modeled_signal =                                     &
     &            mse_get_modeled_signal(this%mse, a_model, sigma,             &
     &                                   last_value, scale_factor)

            CASE (signal_ece_type)
               signal_get_modeled_signal =                                     &
     &            ece_get_modeled_signal(this%ece, a_model, sigma,             &
     &                                   last_value, scale_factor)

            CASE (signal_geometric_type)
               signal_get_modeled_signal =                                     &
     &            geometric_get_modeled_signal(this%geometric,                 &
     &                                         a_model, sigma,                 &
     &                                         last_value, scale_factor)

            CASE (signal_prior_gaussian_type)
!  There is no way to know if a prior changed. Since they are very low cost,
!  just recompute them anyway. So priors don't use the last_value argument.
               signal_get_modeled_signal =                                     &
     &            prior_gaussian_get_modeled_signal(this%gaussian,             &
     &                                              a_model, sigma)
               signal_get_modeled_signal(1) =                                  &
     &            signal_get_modeled_signal(1)*scale_factor

            CASE (signal_feedback_type)
               signal_get_modeled_signal =                                     &
     &            feedback_get_modeled_signal(this%feedback, a_model,          &
     &                                        sigma, last_value,               &
     &                                        scale_factor)

            CASE (signal_combination_type)
               signal_get_modeled_signal =                                     &
     &            combination_get_modeled_signal(this%combination,             &
     &                                           a_model, sigma,               &
     &                                           last_value,                   &
     &                                           scale_factor)

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
!>  @param[in] this      A @ref diagnostic_class instance.
!>  @param[in] a_model   A @ref model instance.
!>  @returns The observed value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_observed_signal(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: signal_get_observed_signal
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(in)     :: a_model

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
!>  @param[inout] this       A @ref diagnostic_class instance.
!>  @param[inout] a_model    A @ref model instance.
!>  @param[in]    use_cache  If the signals have been already calculated this
!>                           should be false. Otherwise set to true.
!>  @param[in]    last_value Last good value in case the signal did not change.
!>  @returns The g^2 value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_g2(this, a_model, use_cache, last_value)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: signal_get_g2
      TYPE (signal_class), INTENT(inout)     :: this
      TYPE (model_class), INTENT(inout)      :: a_model
      LOGICAL, INTENT(in)                    :: use_cache
      REAL (rprec), DIMENSION(4), INTENT(in) :: last_value

!  local Variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  g^2 = W*(O - M)^2/(sigma_o^2 + sigma_m^2) = e*e
      signal_get_g2 =                                                          &
     &   signal_get_e(this, a_model, use_cache, last_value)**2.0

      CALL profiler_set_stop_time('signal_get_g2', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the e contribution of a signal.
!>
!>  Gets the e contribution of an individual signal.
!>  e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
!>
!>  @param[inout] this       A @ref diagnostic_class instance.
!>  @param[inout] a_model    A @ref model instance.
!>  @param[in]    use_cache  If the signals have been already calculated this
!>                           should be false. Otherwise set to true.
!>  @param[in]    last_value Last good value in case the signal did not change.
!>  @returns The e value.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_e(this, a_model, use_cache, last_value)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: signal_get_e
      TYPE (signal_class), INTENT(inout)     :: this
      TYPE (model_class), INTENT(inout)      :: a_model
      LOGICAL, INTENT(in)                    :: use_cache
      REAL (rprec), DIMENSION(4), INTENT(in) :: last_value

!  local Variables
      REAL (rprec), DIMENSION(4)             :: modeled_signal
      REAL (rprec), DIMENSION(4)             :: modeled_sigma
      REAL (rprec)                           :: observed_signal
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      observed_signal = signal_get_observed_signal(this, a_model)

      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma,                &
     &                                           use_cache, last_value)

!  e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
!  Users are can be lax when specifiying the sigma for signals weighted to zero.
!  Avoid a possible divide by zero by error checking the weight first. If the
!  weight is zero, set the error vector explicitly to zero.
      IF (this%weight .eq. 0.0) THEN
         signal_get_e = 0.0
      ELSE
         signal_get_e = (modeled_signal(1) - observed_signal)                  &
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
!>  @see magnetic::magnetic_get_signal_type
!>  @see sxrem::sxrem_get_signal_type
!>  @see intpol::intpol_get_signal_type
!>  @see thomson::thomson_get_signal_type
!>  @see extcurz::extcurz_get_signal_type
!>  @see mse::mse_get_signal_type
!>  @see ece::ece_get_signal_type
!>  @see geometric::geometric_get_signal_type
!>  @see prior_guassian::prior_guassian_get_signal_type
!>  @see feedback::feedback_get_signal_type
!>  @see combination_get_signal_type
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

         CASE (signal_magnetic_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         magnetic_get_signal_type(this%magnetic))

         CASE (signal_sxrem_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         sxrem_get_signal_type(this%sxrem))

         CASE (signal_intpol_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         intpol_get_signal_type(this%intpol))

         CASE (signal_thomson_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         thomson_get_signal_type(this%thomson))

         CASE (signal_extcurz_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         extcurz_get_signal_type(this%extcurz))

         CASE (signal_mse_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         mse_get_signal_type(this%mse))

         CASE (signal_ece_type)
            signal_get_signal_type = TRIM('diagnostic ' //                     &
     &         ece_get_signal_type(this%ece))

         CASE (signal_geometric_type)
            signal_get_signal_type =                                           &
     &         TRIM(geometric_get_signal_type(this%geometric))

         CASE (signal_prior_gaussian_type)
            signal_get_signal_type =                                           &
     &         TRIM(prior_gaussian_get_signal_type(this%gaussian))

         CASE (signal_feedback_type)
            signal_get_signal_type =                                           &
     &         TRIM(feedback_get_signal_type(this%feedback))

         CASE (signal_combination_type)
            signal_get_signal_type =                                           &
     &         TRIM(combination_get_signal_type(this%combination))

         CASE DEFAULT
            signal_get_signal_type = ''

      END SELECT

      CALL profiler_set_stop_time('signal_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_header method of the subclass instance
!>  @ref signal_class was constructed with.
!>
!>  @see magnetic::magnetic_get_header
!>  @see mse::mse_get_header
!>  @see ece::ece_get_header
!>  @see geometric::geometric_get_header
!>  @see feedback::feedback_get_header
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

         CASE (signal_magnetic_type)
            signal_get_header = magnetic_get_header(this%magnetic)

         CASE (signal_mse_type)
            signal_get_header = mse_get_header(this%mse)

         CASE (signal_ece_type)
            signal_get_header = ece_get_header(this%ece)

         CASE (signal_geometric_type)
            signal_get_header = geometric_get_header(this%geometric)

         CASE (signal_feedback_type)
            signal_get_header = feedback_get_header(this%feedback)

         CASE (signal_combination_type)
            signal_get_header = combination_get_header(this%combination)

         CASE DEFAULT
            signal_get_header(1:3) = 'N/A'
            signal_get_header(4) = 'model_sig(1)'
            signal_get_header(5) = 'model_sig(2)'
            signal_get_header(6) = 'model_sig(3)'
            signal_get_header(7) = 'model_sig(4)'

      END SELECT

      CALL profiler_set_stop_time('signal_get_header', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a signal and a position.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp_i method of the subclass instance @ref signal_class was
!>  constructed with.
!>
!>  @see sxrem::sxrem_get_gp_i
!>  @see intpol::intpol_get_gp_i
!>  @see thomson::thomson_get_gp_i
!>  @see ece::ece_get_gp_i
!>  @see prior_gaussian::prior_gaussian_get_gp_i
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] a_model A @ref model_class instance.
!>  @param[in] i       Index of the position for the kernel.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns Kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_gp_i(this, a_model, i, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                    :: signal_get_gp_i
      TYPE (signal_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in)  :: a_model
      INTEGER, INTENT(in)             :: i
      INTEGER, INTENT(in)             :: flags

!  local variables
      REAL (rprec)                    :: scale_factor
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      scale_factor = model_get_signal_factor(a_model, this%scale_index)

      SELECT CASE(this%type)

         CASE (signal_sxrem_type)
            signal_get_gp_i = sxrem_get_gp(this%sxrem, a_model, i,             &
     &                                     flags)

         CASE (signal_intpol_type)
            signal_get_gp_i = intpol_get_gp(this%intpol, a_model, i,           &
     &                                      flags)

         CASE (signal_thomson_type)
            signal_get_gp_i = thomson_get_gp(this%thomson, a_model, i,         &
     &                                       flags)

         CASE (signal_ece_type)
            signal_get_gp_i = ece_get_gp(this%ece, a_model, i, flags)

         CASE (signal_prior_gaussian_type)
            signal_get_gp_i = prior_gaussian_get_gp(this%gaussian,             &
     &                                              a_model, i, flags)

         CASE DEFAULT
            signal_get_gp_i = 0.0

      END SELECT

      signal_get_gp_i = signal_get_gp_i*scale_factor

      CALL profiler_set_stop_time('signal_get_gp_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a signal and a signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp_s method of the subclass instance @ref signal_class was
!>  constructed with. The second signal is cast to a void pointer so it can be
!>  send to the subclass instance.
!>
!>  @see sxrem::sxrem_get_gp_s
!>  @see intpol::intpol_get_gp_s
!>  @see thomson::thomson_get_gp_s
!>  @see ece::ece_get_gp_s
!>  @see prior::prior_get_gp_s
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] a_model A @ref model_class instance.
!>  @param[in] signal  A @ref signal_class instance for the second signal.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_gp_s(this, a_model, signal, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                            :: signal_get_gp_s
      TYPE (signal_class), TARGET, INTENT(in) :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      TYPE (signal_class), TARGET, INTENT(in) :: signal
      INTEGER, INTENT(in)                     :: flags

!  local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (signal_gp_context)                :: gp_context
      REAL (rprec)                            :: start_time
      REAL (rprec)                            :: scale_factor

!  Start of executable code
      start_time = profiler_get_start_time()

      scale_factor = model_get_signal_factor(a_model, this%scale_index)

      gp_context%signal1 => this
      gp_context%signal2 => signal
      gp_context%model => a_model
      gp_context%flags = flags

!  Cast the other signal into a data context. This is the equivalent to casting
!  to a void pointer in C.
      context_length = SIZE(TRANSFER(gp_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(gp_context, context)

      SELECT CASE(this%type)

         CASE (signal_sxrem_type)
            signal_get_gp_s = sxrem_get_gp(this%sxrem, a_model, context,       &
     &                                     signal_kll_callback_p)

         CASE (signal_intpol_type)
            signal_get_gp_s = intpol_get_gp(this%intpol, a_model,              &
     &                                      context,                           &
     &                                      signal_kll_callback_p)

         CASE (signal_thomson_type)
            signal_get_gp_s = thomson_get_gp(this%thomson, context,            &
     &                                       signal_kll_callback_p)

         CASE (signal_ece_type)
            signal_get_gp_s = ece_get_gp(this%ece, a_model, context,           &
     &                                   signal_kll_callback_p)

         CASE (signal_prior_gaussian_type)
            signal_get_gp_s =                                                  &
     &         prior_gaussian_get_gp(this%gaussian, context,                   &
     &                               signal_kll_callback_i, flags)

         CASE DEFAULT
            signal_get_gp_s = 0.0

      END SELECT

      signal_get_gp_s = signal_get_gp_s*scale_factor

      DEALLOCATE(context)

      CALL profiler_set_stop_time('signal_get_gp_s', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a signal and cartesian position.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp_x method of the subclass instance @ref signal_class was
!>  constructed with.
!>
!>  @see sxrem::sxrem_get_gp_x
!>  @see intpol::intpol_get_gp_x
!>  @see thomson::thomson_get_gp_x
!>  @see ece::ece_get_gp_x
!>  @see prior::prior_get_gp_x
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] a_model A @ref model_class instance.
!>  @param[in] x_cart  The cartesian position of to get the kernel at.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns Kernel value for position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_gp_x(this, a_model, x_cart, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: signal_get_gp_x
      TYPE (signal_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(in)         :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: flags

!  local variables
      REAL (rprec)                           :: scale_factor
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      scale_factor = model_get_signal_factor(a_model, this%scale_index)

      SELECT CASE(this%type)

         CASE (signal_sxrem_type)
            signal_get_gp_x = sxrem_get_gp(this%sxrem, a_model, x_cart,        &
     &                                     flags)

         CASE (signal_intpol_type)
            signal_get_gp_x = intpol_get_gp(this%intpol, a_model,              &
     &                                      x_cart, flags)

         CASE (signal_thomson_type)
            signal_get_gp_x = thomson_get_gp(this%thomson, a_model,            &
     &                                       x_cart, flags)

         CASE (signal_ece_type)
            signal_get_gp_x = ece_get_gp(this%ece, a_model, x_cart,            &
     &                                   flags)

         CASE (signal_prior_gaussian_type)
            signal_get_gp_x = prior_gaussian_get_gp(this%gaussian,             &
     &                                              a_model, x_cart,           &
     &                                              flags)

         CASE DEFAULT
            signal_get_gp_x = 0.0

      END SELECT

      signal_get_gp_x = signal_get_gp_x*scale_factor

      CALL profiler_set_stop_time('signal_get_gp_x', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the integration element for the signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_dx method of the subclass instance @ref signal_class was
!>  constructed with.
!>
!>  @see diagnostic::sxrem_get_dx
!>  @see diagnostic::intpol_get_dx
!>
!>  @param[in] this    A @ref signal_class instance.
!>  @param[in] a_model A @ref model_class instance.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns dx value for the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_dx(this, a_model, xcart, dxcart, length, dx,         &
     &                       flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: signal_get_dx
      TYPE (signal_class), INTENT(in)        :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      INTEGER, INTENT(in)                    :: flags

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (signal_sxrem_type)
            signal_get_dx = sxrem_get_dx(this%sxrem, a_model, xcart,           &
     &                                   dxcart, length, dx, flags)

         CASE (signal_intpol_type)
            signal_get_dx = intpol_get_dx(this%intpol, a_model, xcart,         &
     &                                    dxcart, length, dx, flags)

        CASE DEFAULT
           signal_get_dx = 1.0

      END SELECT

      CALL profiler_set_stop_time('signal_get_dx', start_time)

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
!>  @param[inout] this    A @ref diagnostic_class instance.
!>  @param[in]    iou     A input/output representing the file to write to.
!>  @param[in]    index   A index of a signal.
!>  @param[inout] a_model The equilibrium model.
!>  @note This assumes that signals have already been calculated. The model
!>  isn't used by signal_get_g2 since it is reading from the cache but must be
!>  passed in anyway.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write(this, iou, index, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(inout)    :: this
      INTEGER, INTENT(in)                   :: iou
      INTEGER, INTENT(in)                   :: index
      TYPE (model_class), INTENT(inout)     :: a_model

!  local variables
      REAL (rprec), DIMENSION(4)            :: modeled_signal
      REAL (rprec), DIMENSION(4)            :: modeled_sigma
      REAL (rprec)                          :: start_time

!  local paramaters
      REAL (rprec), DIMENSION(4), PARAMETER :: dummy_value = 0.0

!  Start of executable code
      start_time = profiler_get_start_time()

!  When reading from the signal cache, the signal is not recomputed so any value
!  can be used for the last_value arument.
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .true.,        &
     &                                           dummy_value)

      WRITE (iou,1000) index, signal_get_signal_type(this), this%s_name,       &
     &                 signal_get_g2(this, a_model, .true.,                    &
     &                               modeled_signal),                          &
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

         CASE (signal_prior_gaussian_type)
            WRITE (iou,*)
            WRITE (iou,1000) index, signal_get_signal_type(this)
            CALL prior_gaussian_write_auxiliary(this%gaussian, iou,            &
     &                                          a_model)

         CASE (signal_feedback_type)
            CALL feedback_write_auxiliary(this%feedback, iou, index)

         CASE (signal_combination_type)
            CALL combination_write_auxiliary(this%combination,                 &
     &                                       iou, index)

      END SELECT

      CALL profiler_set_stop_time('signals_write_auxiliary', start_time)

1000  FORMAT('Signal',1x,i4,1x,'is a prior signal, type: ',a)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Call back function to generate the kll kernal matrix
!>
!>  This is callback function is needed so that the ith signal can call the j
!>  signal. Need to pass the @ref signal_class as a void pointer otherwise there
!>  is a circular dependancy.
!>
!>  @param[in] context A @ref signal_class instance.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_kll_callback_p(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: signal_kll_callback_p
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (signal_gp_context)               :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context = TRANSFER(context, gp_context)
      signal_kll_callback_p = signal_get_gp(gp_context%signal2,                &
     &                                      gp_context%model, xcart,           &
     &                                      gp_context%flags)                  &
     &                      * signal_get_dx(gp_context%signal1,                &
     &                                      gp_context%model, xcart,           &
     &                                      dxcart, length, dx,                &
     &                                      gp_context%flags)

      CALL profiler_set_stop_time('signal_kll_callback_p', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Call back function to generate the kll kernal matrix
!>
!>  This is callback function is needed so that the ith signal can call the j
!>  signal. Need to pass the @ref signal_class as a void pointer otherwise there
!>  is a circular dependancy.
!>
!>  @param[in] context A @ref signal_class instance.
!>  @param[in] index   Index of the first position.
!>  @returns The kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION signal_kll_callback_i(context, index)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: signal_kll_callback_i
      CHARACTER (len=1), INTENT(in) :: context(:)
      INTEGER, INTENT(in)           :: index

! local variables
      TYPE (signal_gp_context)      :: gp_context
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context = TRANSFER(context, gp_context)
      signal_kll_callback_i = signal_get_gp(gp_context%signal2,                &
     &                                      gp_context%model, index,           &
     &                                      gp_context%flags)

      CALL profiler_set_stop_time('signal_kll_callback_i', start_time)

      END FUNCTION

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the signal data for a step to the result netcdf file.
!>
!>  Writes out the observed, modeled and total sigma value of the signal.
!>
!>  @param[inout] this                  A @ref diagnostic_class instance.
!>  @param[inout] a_model               The equilibrium model.
!>  @param[in]    result_ncid           A NetCDF id of the result file.
!>  @param[in]    current_step          The surrent reconstruction step.
!>  @param[in]    index                 A index of a signal.
!>  @param[in]    signal_model_value_id NetCDF variable id of the model value.
!>  @param[in]    signal_sigma_value_id NetCDF variable id of the total sigma.
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
      TYPE (signal_class), INTENT(inout)    :: this
      TYPE (model_class), INTENT(inout)     :: a_model
      INTEGER, INTENT(in)                   :: result_ncid
      INTEGER, INTENT(in)                   :: current_step
      INTEGER, INTENT(in)                   :: index
      INTEGER, INTENT(in)                   :: signal_model_value_id
      INTEGER, INTENT(in)                   :: signal_sigma_value_id

!  local variables
      INTEGER                               :: status
      INTEGER                               :: varid
      REAL (rprec), DIMENSION(4)            :: modeled_sigma
      REAL (rprec), DIMENSION(4)            :: modeled_signal
      REAL (rprec)                          :: start_time

!  local paramaters
      REAL (rprec), DIMENSION(4), PARAMETER :: dummy_value = 0.0

!  Start of executable code
      start_time = profiler_get_start_time()

!  When reading from the signal cache, the signal is not recomputed so any value
!  can be used for the last_value arument.
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .true.,        &
     &                                           dummy_value)

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
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Syncronize a child signal state to the parent.
!>
!>  Syncs data between a child and the parent process. If MPI support is not
!>  compiled in this subroutine reduces to a no op.
!>
!>  @param[inout] this       A @ref signal_class instance.
!>  @param[in]    index      Reconstruction rank to sync.
!>  @param[in]    recon_comm MPI communicator for the reconstruction processes.
!-------------------------------------------------------------------------------
      SUBROUTINE signal_sync_child(this, index, recon_comm)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                :: index
      INTEGER, INTENT(in)                :: recon_comm

#if defined(MPI_OPT)
!  local variables
      INTEGER                            :: error
      INTEGER                            :: mpi_rank
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      CALL MPI_COMM_RANK(recon_comm, mpi_rank, error)

      IF (mpi_rank .eq. index) THEN
         CALL MPI_SSEND(this%modeled, 4, MPI_REAL8, 0, mpi_rank,               &
     &                  recon_comm, error)
         CALL MPI_SSEND(this%modeled_sigma, 4, MPI_REAL8, 0, mpi_rank,         &
     &                  recon_comm, error)
       ELSE IF (mpi_rank .eq. 0) THEN
         CALL MPI_RECV(this%modeled, 4, MPI_REAL8, index, index,               &
     &                 recon_comm, MPI_STATUS_IGNORE, error)
         CALL MPI_RECV(this%modeled_sigma, 4, MPI_REAL8, index, index,         &
     &                 recon_comm, MPI_STATUS_IGNORE, error)
      END IF

      CALL profiler_set_stop_time('signal_sync_child', start_time)
#endif

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
!>                              values are @ref combination_sum,
!>                              @ref combination_max, @ ref combination_min and
!>                              @ref combination_wavg
!>  @param[in] wgt_index        Index number of the average weighting.
!>  @returns A pointer to a constructed @ref combination_class object.
!-------------------------------------------------------------------------------
      FUNCTION combination_construct(n_signals, combination_type,              &
     &                               wgt_index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), POINTER :: combination_construct
      INTEGER, INTENT(in)               :: n_signals
      CHARACTER (len=combination_type_length), INTENT(in)                      &
     &                                  :: combination_type
      INTEGER, INTENT(IN)               :: wgt_index

!  Local Variables
      REAL (rprec)                      :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(combination_construct)

      ALLOCATE(combination_construct%signals(n_signals))
      ALLOCATE(combination_construct%a(n_signals))

      combination_construct%wgt_index = wgt_index

      SELECT CASE (combination_type)
         CASE ('sum')
            combination_construct%type = combination_sum

         CASE ('max')
            combination_construct%type = combination_max

         CASE ('min')
            combination_construct%type = combination_min

         CASE ('wavg')
            combination_construct%type = combination_wavg
            IF (n_signals .gt. 2) THEN
               WRITE (*,1000)
            END IF

      END SELECT

      CALL profiler_set_stop_time('combination_construct', start_time)

1000  FORMAT('More than two signals specified for weighted average.            &
     &       Only using the first two.')

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

!-------------------------------------------------------------------------------
!>  @brief Sets the weight factors from from the weight index.
!>
!>  Only assumed valid when the combination type is weighted average. This also
!>  assumes that the average of only two signals.
!>
!>  @param[inout] this  A @ref combination_class instance.
!>  @param[in]  a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE combination_set_weight(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(in)          :: a_model

!  local Variables
      REAL (rprec)                             :: start_time

      this%a(1) = a_model%coosig_wgts(this%wgt_index)
      this%a(2) = 1.0 - this%a(1)

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
!>  @param[in]    this         A @ref combination_class instance.
!>  @param[inout] a_model      A @ref model instance.
!>  @param[out]   sigma        The modeled sigma.
!>  @param[in]    last_value   Last good value in case the signal did not change.
!>  @param[in]    scale_factor Factor to scale the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION combination_get_modeled_signal(this, a_model, sigma,            &
     &                                        last_value, scale_factor)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)      :: combination_get_modeled_signal
      TYPE (combination_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor

!  Local Variables
      INTEGER                                 :: i
      REAL (rprec), DIMENSION(4)              :: sigma_local
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      combination_get_modeled_signal =                                         &
     &   this%a(1)*signal_get_modeled_signal(this%signals(1)%p, a_model,       &
     &                                       sigma, .true., last_value)

      SELECT CASE(this%type)
         CASE (combination_sum)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal =                                &
     &            combination_get_modeled_signal +                             &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.,                  &
     &                                                last_value)
               sigma = sigma + sigma_local
            END DO

         CASE (combination_max)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal = MAX(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.,                  &
     &                                                last_value))
               sigma = MAX(sigma, sigma_local)
            END DO

         CASE (combination_min)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal = MIN(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.,                  &
     &                                                sigma_local))
               sigma = MIN(sigma, sigma_local)
            END DO

!  ECH I assume that the combination wavg is composed of only two signals and
!  the there is only 1 coosig wgt. Second do I need to account for the
!  uncertainity in the wgt_ids when calculaing the sigma value.
         CASE (combination_wavg)
            CALL combination_set_weight(this, a_model)
            combination_get_modeled_signal =                                   &
     &         combination_get_modeled_signal +                                &
     &         this%a(2)*signal_get_modeled_signal(this%signals(2)%p,          &
     &                                             a_model, sigma_local,       &
     &                                             .true., last_value)

            sigma = sigma + this%a(2)*sigma_local

!  It doesn't always make sense to average the values in model signal 2-4.
            combination_get_modeled_signal(2:4) = 0.0

      END SELECT

      combination_get_modeled_signal(1) =                                      &
     &   combination_get_modeled_signal(1)*scale_factor

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

         CASE (combination_wavg)
            combination_get_signal_type = 'coosig wavg'

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
      IF (this%type .eq. combination_wavg) THEN
         WRITE (iou,1001) this%wgt_index
      END IF
      WRITE (iou,1002)

      DO i = 1, SIZE(this%signals)
         WRITE (iou,1003) i, this%signals(i)%p%s_name, this%a(i)
      END DO

      CALL profiler_set_stop_time('combination_write_auxiliary',               &
     &                            start_time)

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a combination of other signals, type: ',a)
1001  FORMAT('Weighted Average Index: ',i4)
1002  FORMAT('term #',2x,'s_name',16x,'Coefficient')
1003  FORMAT(2x,i4,2x,a20,2x,es12.5)

      END SUBROUTINE

      END MODULE
