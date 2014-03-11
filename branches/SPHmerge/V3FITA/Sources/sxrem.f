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
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a magnetic signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE sxrem_class
         TYPE (path) :: chord_path
!>  Index of the soft x-ray emissivity profile for this chord.
         INTEGER     :: profile_number
      END TYPE sxrem_class

!-------------------------------------------------------------------------------
!> Structure to hold all memory needed to be sent to the callback function.
!-------------------------------------------------------------------------------
      TYPE sxrem_context
         INTEGER                     :: profile_number
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
!>  @param[in] start_path Starting point of a sxrem chord.
!>  @param[in] end_path Ending point of a sxrem chord.
!>  @returns A pointer to a constructed @ref sxrem_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct(start_path, end_path, profile_number)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER :: sxrem_construct
      REAL(rprec), DIMENSION(3)   :: start_path
      REAL(rprec), DIMENSION(3)   :: end_path
      INTEGER                     :: profile_number

!  Start of executable code
      ALLOCATE(sxrem_construct)

      CALL path_append_vertex(sxrem_construct%chord_path, start_path(1),       &
     &                        start_path(2), start_path(3))
      CALL path_append_vertex(sxrem_construct%chord_path, end_path(1),         &
     &                        end_path(2), end_path(3))

      sxrem_construct%profile_number = profile_number

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
      TYPE (sxrem_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      CALL path_destroy(this%chord_path)
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
!>  @param[in]    this    A @ref sxrem_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_signal
      TYPE (sxrem_class), INTENT(inout)       :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

! local variables
      CHARACTER(len=1), ALLOCATABLE  :: context(:)
      INTEGER                        :: context_length
      TYPE (sxrem_context)           :: sxr_context

!  Start of executable code
      sxrem_get_modeled_signal = 0.0

!  The relevant data for the soft x-ray context.
      sxr_context%profile_number = this%profile_number
      sxr_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
      context_length = SIZE(TRANSFER(sxr_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(sxr_context, context)
      sxrem_get_modeled_signal(1) = path_integrate(this%chord_Path,            &
     &                                             sxr_function,               &
     &                                             context)

      sigma = 0.0
      DEALLOCATE(context)

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

!  Start of executable code
      WRITE (sxrem_get_signal_type, 1000) this%profile_number

1000  FORMAT('sxrch(',i2,')')

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Soft x-ray callback function.
!>
!>  Returns the value of the soft x-ray emissivity times the change in path
!>  length. This function is passed to @ref path_integrate to act as a callback.
!>  The soft x-ray emissivity is provided by @ref model::model_get_sxrem.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref the model.
!>  @param[in] xcart   A the integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] dx      A scalar change in path.
!>  @returns The sxrem(x)*dl at point x.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C. 
!-------------------------------------------------------------------------------
      FUNCTION sxr_function(context, xcart, dxcart, dx)
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=1), INTENT(in)          :: context(:)
      REAL(rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL(rprec), INTENT(in)               :: dx
      REAL(rprec)                           :: sxr_function

! local variables
      TYPE(sxrem_context)                   :: sxr_context

!  Start of executable code
      sxr_context = TRANSFER(context, sxr_context)
      sxr_function = model_get_sxrem(sxr_context%model, xcart,                 &
     &                               sxr_context%profile_number)*dx

      END FUNCTION

      END MODULE
