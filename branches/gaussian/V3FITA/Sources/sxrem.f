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
!>  Reference to a @ref model_class object.
         TYPE (model_class), POINTER :: model => null()
      END TYPE

!-------------------------------------------------------------------------------
!> Structure to hold all memory needed to be sent to the callback function.
!-------------------------------------------------------------------------------
      TYPE gpsxrem_context
!>  The index of the kernel function
         INTEGER                     :: kernel_id = -1
!>  Reference to a @ref model_class object.
         TYPE (model_class), POINTER :: model => null()
!>  Cartesian position of the 1st position when integrating the 2nd chord
         REAL (rprec), DIMENSION(3)  :: xpoint = 0.0
!>  Radial position of the model position when integrating chord in KLS eval
         REAL (rprec)                :: spoint = 0.0
!>  Geometric chord factor of the 2nd chord
         REAL (rprec)                :: geo = 0.0
!>  The complete path of the 2nd chord
         TYPE (vertex), POINTER :: chord_path => null()
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
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (sxrem_context)                    :: sxr_context
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_get_modeled_signal = 0.0

!  The relevant data for the soft x-ray context.
      sxr_context%profile_number = this%profile_number
      sxr_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
      context_length = SIZE(TRANSFER(sxr_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(sxr_context, context)
      sxrem_get_modeled_signal(1) = path_integrate(this%chord_path,            &
     &                                             sxr_function,               &
     &                                             context)                    &
     &                            * this%geo

      sigma = 0.0
      DEALLOCATE(context)

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
!  CALL BACK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Soft x-ray callback function for the gaussian process integration
!>
!>  Returns the value of the gaussian process kernel functon times the change in
!>  the path length. This function is passed to @ref path_integrate to act as a
!>  callback. The kernel function is provided by @ref model::model_get_gpker.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref gpsxrem_context for the model and 2nd point
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The gpker(x,s1)*dl at point x. s1 i
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C. 
!-------------------------------------------------------------------------------
      FUNCTION gpsxr_kls_function(context, xcart, dxcart, length, dx)
      USE coordinate_utilities

      IMPLICIT NONE

!     Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: gpsxr_kls_function

! local variable
      TYPE(gpsxrem_context)                  :: gpsxr_context
      REAL(rprec)                            :: start_time

      gpsxr_context = TRANSFER(context, gpsxr_context)
      gpsxr_kls_function = model_get_gpker(gpsxr_context%model,         &
     &                    xcart, gpsxr_context%spoint,                  &
     &                    gpsxr_context%kernel_id)*dx

      CALL profiler_set_stop_time('gpsxr_kls_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Soft x-ray callback function for the gaussian process integration
!>
!>  Returns the value of the gaussian process kernel functon integrated along
!>  a second chord times the change in path lenth times the 2nd chords geometric
!>  factor. This function is passed to @ref path_integrate to act as a
!>  callback. The kernel function is provided by @ref model::model_get_gpker.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref gpsxrem_context for the model and 2nd integration 
!>  path.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The integral gpker(x,x2)*dl*geo at point x along the path x2.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C. 
!-------------------------------------------------------------------------------
      FUNCTION gpsxr_kll1_function(context, xcart, dxcart, length, dx)
      USE coordinate_utilities

      IMPLICIT NONE

!     Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: gpsxr_kll1_function

! local variable
      TYPE(gpsxrem_context)                  :: gpsxr_context
      INTEGER                                :: context2_length
      CHARACTER(len=1), ALLOCATABLE          :: context2(:)
      REAL(rprec)                            :: start_time

      REAL(rprec), PARAMETER                 :: dxinput = 0.005

!  Start of executable code
      start_time = profiler_get_start_time()
      gpsxr_kll1_function = 0.0

! Add the current position to gp soft x-ray context. And then recast into a new
! context
      gpsxr_context = TRANSFER(context, gpsxr_context)
      gpsxr_context%xpoint = xcart

      context2_length = SIZE(TRANSFER(gpsxr_context, context2))
      ALLOCATE(context2(context2_length))
      context2 = TRANSFER(gpsxr_context, context2)
! ech temp increase step size to speed up integration
      gpsxr_kll1_function = path_integrate(gpsxr_context%chord_path,    &
     &                                             gpsxr_kll2_function, &
!     &                                             context2)            &
     &                                             context2,dxinput)     &
     &                            * gpsxr_context%geo*dx

      DEALLOCATE(context2)

      CALL profiler_set_stop_time('gpsxr_kll1_function', start_time)

      END FUNCTION
!-------------------------------------------------------------------------------
!>  @brief Soft x-ray callback function for the gaussian process integration
!>
!>  Returns the value of the gaussian process kernel functon times the change in
!>  the path length. This function is passed to @ref path_integrate to act as a
!>  callback. The kernel function is provided by @ref model::model_get_gpker.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref gpsxrem_context for the model and 2nd point
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The gpker(x,x2)*dl at point x. s2 is the 2nd point
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C. 
!-------------------------------------------------------------------------------
      FUNCTION gpsxr_kll2_function(context, xcart, dxcart, length, dx)
      USE coordinate_utilities

      IMPLICIT NONE

!     Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: gpsxr_kll2_function

! local variable
      TYPE(gpsxrem_context)                  :: gpsxr_context
      REAL(rprec)                            :: start_time

      gpsxr_context = TRANSFER(context, gpsxr_context)
      gpsxr_kll2_function = model_get_gpker(gpsxr_context%model,        &
     &                    xcart, gpsxr_context%xpoint,                  &
     &                    gpsxr_context%kernel_id)*dx

      CALL profiler_set_stop_time('gpsxr_kll2_function', start_time)

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
!>  @param[in] context A @ref sxrem_context for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The sxrem(x)*dl at point x.
!>
!>  @note context in this case if not that same as @ref context. It refers to
!>  generic passed to the function. This is the fortran eqivalent to a void
!>  pointer in C. 
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
