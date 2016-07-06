!*******************************************************************************
!>  @file ece_travis.f
!>  @brief Contains module @ref ece_travis.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref ece_travis_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE ece_travis

      USE stel_kinds, only: rprec
      USE integration_path
      USE model

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) ece_travis base class
!  2) ece_travis context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing an ECE signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE ece_travis_class
!>  TRAVIS Data File name
         CHARACTER (len=path_length) :: file
!>  The index of the ece signal in the TRAVIS data_file
         INTEGER                     :: index
      END TYPE ece_travis_class

!-------------------------------------------------------------------------------
!> Structure to hold all memory needed to be sent to the callback function.
!-------------------------------------------------------------------------------
      TYPE ece_travis_context
!>  The index of the emissivity profile model.
         REAL (rprec)                :: resonance
!>  Reference to a @ref model_class object.
         TYPE (model_class), POINTER :: model => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PRIVATE :: ece_function, is_in_range

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref ece_travis_class object.
!>
!>  Allocates memory and initializes a @ref ece_travis_class object.
!>
!>  @param[in] data_file TRAVIS data file where ece signal is defined.
!>  @param[in] index Index of ECE signal in TRAVIS data file.
!>  @returns A pointer to a constructed @ref ece_travis_class object.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_construct(data_file, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (ece_travis_class), POINTER        :: ece_travis_construct
      CHARACTER (len=path_length), INTENT(in) :: data_file
      INTEGER, INTENT(in)                     :: index
!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(ece_travis_construct)

      ece_travis_construct%file = data_file
      ece_travis_construct%index = index

      CALL profiler_set_stop_time('ece_travis_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref ece_travis_class object.
!>
!>  Deallocates memory and uninitializes a @ref ece_class object.
!>
!>  @param[inout] this A @ref ece_travis_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE ece_travis_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (ece_travis_class), POINTER :: this

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the ECE signal by searching along the chord for the resonance.
!>  The search criteria is provided by the @ref ece_function.
!>
!>  @param[in]    this    A @ref ece_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: ece_travis_get_modeled_signal
      TYPE (ece_travis_class), INTENT(inout)         :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

! local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (ece_travis_context)                      :: temp_context
      LOGICAL                                 :: found
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ece_travis_get_modeled_signal = 0.0

!  The relevant data for the ece context.
!      temp_context%resonance = this%resonance
!      temp_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
!      context_length = SIZE(TRANSFER(temp_context, context))
!      ALLOCATE(context(context_length))
!      context = TRANSFER(temp_context, context)
!         ece_get_modeled_signal(2:4) = path_search(this%chord_path,            &
!     &                                             ece_function,               &
!     &                                             context, found)

!      IF (found) THEN
!         ece_get_modeled_signal(1) =                                           &
!     &      model_get_te(a_model, ece_get_modeled_signal(2:4))
!      END IF

      sigma = 0.0
!      DEALLOCATE(context)

      CALL profiler_set_stop_time('ece_travis_get_modeled_signal',             &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the ece type.
!>
!>  Returns a description of the ece type for use when writting output files.
!>
!>  @param[in] this A @ref ece_class instance.
!>  @returns A string describing the ece type.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: ece_travis_get_signal_type
      TYPE (ece_travis_class), INTENT(in)     :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ece_travis_get_signal_type = 'ece_travis'

      CALL profiler_set_stop_time('ece_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  Returns a description of the array indices for use when writting output
!>  files.
!>
!>  @param[in] this A @ref ece_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION ece_travis_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &                                 ece_travis_get_header
      TYPE (ece_travis_class), INTENT(in)                   :: this

!  local variables
      REAL (rprec)                                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ece_travis_get_header(1) = 'x (m)'
      ece_travis_get_header(2) = 'y (m)'
      ece_travis_get_header(3) = 'z (m)'

      ece_travis_get_header(4) = 'model_sig(1)'
      ece_travis_get_header(5) = 'model_sig(2)'
      ece_travis_get_header(6) = 'model_sig(3)'
      ece_travis_get_header(7) = 'model_sig(4)'

      CALL profiler_set_stop_time('ece_travis_get_header', start_time)

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief ECE callback function.
!>
!>  Returns true if the magnetic field resonance if found with in a specified
!>  range.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref ece_context for the model.
!>  @param[in] xcart   The search point.
!>  @returns True if the resonace was found.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C.
!-------------------------------------------------------------------------------
      FUNCTION ece_function(context, xcart1, xcart2)
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart1
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart2
      LOGICAL                                :: ece_function

! local variables
      TYPE (ece_travis_context)                     :: temp_context
      REAL (rprec), DIMENSION(3)             :: bcart1
      REAL (rprec), DIMENSION(3)             :: bcart2
      REAL (rprec)                           :: bmod1
      REAL (rprec)                           :: bmod2
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ece_function = .false.

      temp_context = TRANSFER(context, temp_context)

      bcart1 = equilibrium_get_B_vec(temp_context%model%equilibrium,           &
     &                               xcart1, .false.)
      bcart2 = equilibrium_get_B_vec(temp_context%model%equilibrium,           &
     &                               xcart2, .false.)

      bmod1 = SQRT(DOT_PRODUCT(bcart1, bcart1))
      bmod2 = SQRT(DOT_PRODUCT(bcart2, bcart2))

      IF (ABS(bmod1 - bmod2) .lt.                                              &
     &    temp_context%model%resonace_range) THEN
         ece_function = is_in_range(temp_context%resonance,                    &
     &                              bmod1, bmod2)
      END IF

      CALL profiler_set_stop_time('ece_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check if values is in range.
!>
!>  Returns true the value if x_low <= xp <= x_high
!>
!>  @param[in] xp Value to see if it is in range.
!>  @param[in] x1 First bounding point to check.
!>  @param[in] x2 Second bounding point to check.
!>  @returns True if the resonace was found.
!-------------------------------------------------------------------------------
      PURE FUNCTION is_in_range(xp, x1, x2)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                 :: is_in_range
      REAL(rprec), INTENT(in) :: xp
      REAL(rprec), INTENT(in) :: x1
      REAL(rprec), INTENT(in) :: x2

!  Start of executable code
      IF (x1 .gt. x2) THEN
         is_in_range = xp .ge. x2 .and. xp .le. x1
      ELSE
         is_in_range = xp .ge. x1 .and. xp .le. x2
      END IF

      END FUNCTION

      END MODULE
