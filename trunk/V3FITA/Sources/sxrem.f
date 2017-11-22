!*******************************************************************************
!>  @file sxrem.f
!>  @brief Contains module @ref sxrem.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref sxrem_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE sxrem

      USE stel_kinds, only: rprec
      USE integration_path
      USE model

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) sxrem base class
!  2) sxrem context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a soft x-ray emissivity signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE sxrem_class
!>  The complete path of the chord.
         TYPE (vertex), POINTER :: chord_path => null()
!>  Geometric chord factor.
         REAL (rprec)           :: geo
!>  Index of the soft x-ray emissivity profile for this chord.
         INTEGER                :: profile_number
      END TYPE sxrem_class

!-------------------------------------------------------------------------------
!> Structure to hold all memory needed to be sent to the callback function.
!-------------------------------------------------------------------------------
      TYPE sxrem_context
!>  The index of the emissivity profile model.
         INTEGER                     :: profile_number
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PRIVATE :: sxr_function

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_class object.
!>
!>  Allocates memory and initializes a @ref sxrem_class object.
!>
!>  @param[in] start_path     Starting point of a sxrem chord.
!>  @param[in] end_path       Ending point of a sxrem chord.
!>  @param[in] geo            Geometric factor of the chord.
!>  @param[in] profile_number Index of the soft x-ray emissivity profile.
!>  @returns A pointer to a constructed @ref sxrem_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct(start_path, end_path, geo,                      &
     &                         profile_number)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER            :: sxrem_construct
      REAL (rprec), DIMENSION(3), INTENT(in) :: start_path
      REAL (rprec), DIMENSION(3), INTENT(in) :: end_path
      REAL (rprec), INTENT(in)               :: geo
      INTEGER, INTENT(in)                    :: profile_number

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(sxrem_construct)

      CALL path_append_vertex(sxrem_construct%chord_path, start_path)
      CALL path_append_vertex(sxrem_construct%chord_path, end_path)

      sxrem_construct%geo = geo
      sxrem_construct%profile_number = profile_number

      CALL profiler_set_stop_time('sxrem_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref sxrem_class object.
!>
!>  Deallocates memory and uninitializes a @ref sxrem_class object.
!>
!>  @param[inout] this A @ref sxrem_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE sxrem_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%chord_path)) THEN
         CALL path_destruct(this%chord_path)
         this%chord_path => null()
      END IF

      this%profile_number = 0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the soft x-ray emissivity signal by integrating along the chord
!>  path. The emissivity is provided by the @ref sxr_function.
!>
!>  @param[in]  this       A @ref sxrem_class instance.
!>  @param[in]  a_model    A @ref model instance.
!>  @param[out] sigma      The modeled sigma.
!>  @param[in]  last_value Last good value in case the signal did not change.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_signal(this, a_model, sigma,                  &
     &                                  last_value)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_signal
      TYPE (sxrem_class), INTENT(in)          :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value

! local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (sxrem_context)                    :: sxr_context
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_get_modeled_signal = last_value
      sigma = 0.0

      IF (BTEST(a_model%state_flags, model_state_vmec_flag) .or.               &
     &    BTEST(a_model%state_flags, model_state_siesta_flag) .or.             &
     &    BTEST(a_model%state_flags, model_state_shift_flag) .or.              &
     &    BTEST(a_model%state_flags, model_state_sxrem_flag +                  &
     &                               (this%profile_number - 1))) THEN

!  The relevant data for the soft x-ray context.
         sxr_context%profile_number = this%profile_number
         sxr_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
         context_length = SIZE(TRANSFER(sxr_context, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(sxr_context, context)
         sxrem_get_modeled_signal(1) = path_integrate(this%chord_path,         &
     &                                                sxr_function,            &
     &                                                context)                 &
     &                               * this%geo

         DEALLOCATE(context)
      END IF

      CALL profiler_set_stop_time('sxrem_get_modeled_signal',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the sxrem type.
!>
!>  Returns a description of the sxrem type for use when writting output files.
!>
!>  @param[in] this A @ref sxrem_class instance.
!>  @returns A string describing the sxrem type.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: sxrem_get_signal_type
      TYPE (sxrem_class), INTENT(in)   :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (sxrem_get_signal_type, 1000) this%profile_number

      CALL profiler_set_stop_time('sxrem_get_signal_type', start_time)

1000  FORMAT('sxrch(',i2,')')

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Soft x-ray callback function.
!>
!>  Returns the value of the soft x-ray emissivity times the change in path
!>  length. This function is passed to @ref integration_path::path_integrate to
!>  act as a callback. The soft x-ray emissivity is provided by
!>  @ref model::model_get_sxrem.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref sxrem_context for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The sxrem(x)*dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION sxr_function(context, xcart, dxcart, length, dx)
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: sxr_function

! local variables
      TYPE (sxrem_context)                   :: sxr_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxr_context = TRANSFER(context, sxr_context)
      sxr_function = model_get_sxrem(sxr_context%model, xcart,                 &
     &                               sxr_context%profile_number)*dx

      CALL profiler_set_stop_time('sxr_function', start_time)

      END FUNCTION

      END MODULE
