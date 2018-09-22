!*******************************************************************************
!>  @file intpol.f
!>  @brief Contains module @ref intpol.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Implements interferometry/polarimetry diagnostic. Defines the base class of
!>  the type @ref intpol_class.
!>  @par Super Class:
!>  @ref diagnostic
!*******************************************************************************

      MODULE intpol

      USE stel_kinds, only: rprec, dp
      USE integration_path
      USE coordinate_utilities
      USE profiler
      USE model

      IMPLICIT NONE
!*******************************************************************************
!  intpol module parameters
!*******************************************************************************
!>  Type descriptor for intpol no type.
      INTEGER, PARAMETER :: intpol_no_type             = -1
!>  Type descriptor for intpol type represented a interferometer diagnostic.
      INTEGER, PARAMETER :: intpol_interferometry_type = 0
!>  Type descriptor for intpol type represented a polarmetry diagnostic.
      INTEGER, PARAMETER :: intpol_polarimetry_type    = 1

!>  Constant term for the polarimety.
      REAL(rprec), PARAMETER :: intpol_polar_constant = 2.62E-13_dp

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) intpol base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a intpol signal.
!>  @par Super Class:
!>  @ref diagnostic
!-------------------------------------------------------------------------------
      TYPE intpol_class
!>  Type descirptor of the intpol type.
!>  @par Possible values are:
!>  * @ref intpol_no_type
!>  * @ref intpol_interferometry_type
!>  * @ref intpol_polarimetry_type
         INTEGER                :: type
!>  Wavelength if the beam.
         REAL (rprec)           :: wavelength
!>  Controls if the results is computed in degrees or radians.
         LOGICAL                :: in_degrees
!>  The complete path of the chord.
         TYPE (vertex), POINTER :: chord_path => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the guassian process
!>  callback function.
!-------------------------------------------------------------------------------
      TYPE intpol_gp_context_i
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
!>  Position index.
         INTEGER                     :: i
      END TYPE

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the guassian process
!>  callback function for signal.
!-------------------------------------------------------------------------------
      TYPE intpol_gp_context_x
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
!>  First position.
         REAL (rprec), DIMENSION(3)  :: xcart
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref intpol_class types using
!>  @ref intpol_construct_int or @ref intpol_construct_pol.
!-------------------------------------------------------------------------------
      INTERFACE intpol_construct
         MODULE PROCEDURE intpol_construct_int,                                &
     &                    intpol_construct_pol
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to get the guassian process kernel values.
!-------------------------------------------------------------------------------
      INTERFACE intpol_get_gp
         MODULE PROCEDURE intpol_get_gp_i,                                     &
     &                    intpol_get_gp_s,                                     &
     &                    intpol_get_gp_x
      END INTERFACE

      PRIVATE :: int_function, pol_function, gp_function_i,                    &
     &           gp_function_x

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref intpol_class object representing a interferometer
!>  diagnostic.
!>
!>  Allocates memory and initializes a @ref intpol_class object.
!>
!>  @param[in] chord_paths The nodes of a multi point segment integration path.
!>  @returns A pointer to a constructed @ref intpol_class object.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_int(chord_paths)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: intpol_construct_int
      REAL (rprec), DIMENSION(:,:) :: chord_paths

!  local variables
      INTEGER                      :: i
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(intpol_construct_int)

      DO i = 1, SIZE(chord_paths, 1)
         CALL path_append_vertex(intpol_construct_int%chord_path,              &
     &                           chord_paths(i,:))
      END DO

      intpol_construct_int%type = intpol_interferometry_type
      intpol_construct_int%wavelength = 0.0
      intpol_construct_int%in_degrees = .false.

      CALL profiler_set_stop_time('intpol_construct_int', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref intpol_class object representing a polarimetry
!>  diagnostic.
!>
!>  Allocates memory and initializes a @ref intpol_class object.
!>
!>  @param[in] wavelength  Wavelength of the polarimetry beam.
!>  @param[in] in_degrees  Specifies if the sigals is degrees or radians. True
!>                         specifies degree while false specifies radians.
!>  @param[in] chord_paths The nodes of a multi point segment integration path.
!>  @returns A pointer to a constructed @ref intpol_class object.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct_pol(wavelength, in_degrees, chord_paths)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: intpol_construct_pol
      REAL (rprec), INTENT(in)     :: wavelength
      LOGICAL, INTENT(in)          :: in_degrees
      REAL (rprec), DIMENSION(:,:) :: chord_paths

!  local variables
      INTEGER                      :: i
      REAL (rprec)                 :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(intpol_construct_pol)

      DO i = 1, SIZE(chord_paths, 1)
         CALL path_append_vertex(intpol_construct_pol%chord_path,              &
     &                           chord_paths(i,:))
      END DO

      intpol_construct_pol%type = intpol_polarimetry_type
      intpol_construct_pol%wavelength = wavelength
      intpol_construct_pol%in_degrees = in_degrees

      CALL profiler_set_stop_time('intpol_construct_pol', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref intpol_class object.
!>
!>  Deallocates memory and uninitializes a @ref intpol_class object.
!>
!>  @param[inout] this A @ref intpol_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE intpol_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), POINTER :: this

!  Start of executable code
      this%type = intpol_no_type
      this%wavelength = 0.0
      this%in_degrees = .false.

      IF (ASSOCIATED(this%chord_path)) THEN
         CALL path_destruct(this%chord_path)
         this%chord_path => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  Calculates the modeled signal by integrating along the chord path. The
!>  integration function is provided by @ref int_function for interferometer
!>  diagnostics and @ ref pol_function for polarimetry diagnostics.
!>
!>  @param[in]  this          A @ref intpol_class instance.
!>  @param[in]  a_model       A @ref model instance.
!>  @param[out] sigma         The modeled sigma.
!>  @param[in]  last_value    Last good value in case the signal did not change.
!>  @param[in]  scale_factor  Factor to scale the modeled signal.
!>  @param[in]  offset_factor Factor to offset the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_modeled_signal(this, a_model, sigma,                 &
     &                                   last_value, scale_factor,             &
     &                                   offset_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: intpol_get_modeled_signal
      TYPE (intpol_class), INTENT(in)         :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor
      REAL (rprec), INTENT(in)                :: offset_factor

!  local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      intpol_get_modeled_signal = last_value
      sigma = 0.0

      IF (BTEST(a_model%state_flags, model_state_vmec_flag)   .or.             &
     &    BTEST(a_model%state_flags, model_state_siesta_flag) .or.             &
     &    BTEST(a_model%state_flags, model_state_ne_flag)     .or.             &
     &    BTEST(a_model%state_flags, model_state_shift_flag)  .or.             &
     &    BTEST(a_model%state_flags, model_state_signal_flag)) THEN

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
         context_length = SIZE(TRANSFER(a_model, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(a_model, context)

         SELECT CASE (this%type)

            CASE (intpol_interferometry_type)
               intpol_get_modeled_signal(1) =                                  &
     &            path_integrate(a_model%int_params, this%chord_path,          &
     &                           int_function, context)

            CASE (intpol_polarimetry_type)
               intpol_get_modeled_signal(1) = intpol_polar_constant *          &
     &            (this%wavelength**2.0_dp) *                                  &
     &            path_integrate(a_model%int_params, this%chord_path,          &
     &                           pol_function, context)
               IF (this%in_degrees) THEN
                   intpol_get_modeled_signal(1) =                              &
     &                intpol_get_modeled_signal(1)/degree
               ENDIF

         END SELECT

         intpol_get_modeled_signal(1) =                                        &
     &      (intpol_get_modeled_signal(1) + offset_factor)*scale_factor

         DEALLOCATE(context)
      END IF

      CALL profiler_set_stop_time('intpol_get_modeled_signal',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the intpol type.
!>
!>  Returns a description of the intpol type for use when writting output files.
!>
!>  @param[in] this A @ref intpol_class instance.
!>  @returns A string describing the intpol type.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: intpol_get_signal_type
      TYPE (intpol_class), INTENT(in)  :: this

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      intpol_get_signal_type = 'ipch'

      CALL profiler_set_stop_time('intpol_get_signal_type', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for an inteferometry/polarimetry
!>  signal and a position.
!>
!>  Calculates the guassian process kernel between the signal and the position.
!>  Density kernels are provided by @ref model::model_get_gp_ne. Only valid for
!>  interferometry types.
!>
!>  @param[in] this         A @ref intpol_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] i            Index of the position for the kernel.
!>  @param[in] flags        State flags to send to the kernel.
!>  @returns Kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_gp_i(this, a_model, i, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: intpol_get_gp_i
      TYPE (intpol_class), INTENT(in)        :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      INTEGER, INTENT(in)                    :: i
      INTEGER, INTENT(in)                    :: flags

!  local variables
      CHARACTER(len=1), ALLOCATABLE          :: context(:)
      INTEGER                                :: context_length
      TYPE (intpol_gp_context_i)             :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the guassian process context.
      gp_context%model => a_model
      gp_context%i = i

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
      context_length = SIZE(TRANSFER(gp_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(gp_context, context)

      SELECT CASE (this%type)

         CASE (intpol_interferometry_type)
            intpol_get_gp_i =                                                  &
     &         path_integrate(a_model%int_params, this%chord_path,             &
     &                        gp_function_i, context)

         CASE (intpol_polarimetry_type)
            intpol_get_gp_i = intpol_polar_constant                            &
     &                      * (this%wavelength**2.0_dp)                        &
     &                      * path_integrate(a_model%int_params,               &
     &                                       this%chord_path,                  &
     &                                       gp_function_pol_i,                &
     &                                       context)

            IF (this%in_degrees) THEN
               intpol_get_gp_i = intpol_get_gp_i/degree
            END IF

         CASE DEFAULT
            intpol_get_gp_i = 0.0

      END SELECT

      DEALLOCATE(context)

      CALL profiler_set_stop_time('intpol_get_gp_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for an inteferometry signal and a
!>  signal.
!>
!>  Calculates the guassian process kernel between the signal and a signal.
!>  Calls back to the @ref signal module to call the other signal. Only valid
!>  for interferometry types. This is the first signal.
!>
!>  @param[in] this         A @ref intpol_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] gp_context   A @ref signal_class instance.
!>  @param     gp_callback  Function pointer to call back to signal module.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_gp_s(this, a_model, context, callback)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                    :: intpol_get_gp_s
      TYPE (intpol_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in)  :: a_model
      CHARACTER (len=1), INTENT(in)   :: context(:)
      INTERFACE
         FUNCTION callback(context, xcart, dxcart, length, dx)
         USE stel_kinds
         REAL (rprec)                           :: callback
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
         REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
         REAL (rprec), INTENT(in)               :: length
         REAL (rprec), INTENT(in)               :: dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec)                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (intpol_interferometry_type)
            intpol_get_gp_s = path_integrate(a_model%int_params,               &
     &                                       this%chord_path, callback,        &
     &                                       context)

         CASE (intpol_polarimetry_type)
            intpol_get_gp_s = intpol_polar_constant                            &
     &                      * (this%wavelength**2.0_dp)                        &
     &                      * path_integrate(a_model%int_params,               &
     &                                       this%chord_path, callback,        &
     &                                       context)

            IF (this%in_degrees) THEN
               intpol_get_gp_s = intpol_get_gp_s/degree
            ENDIF

         CASE DEFAULT
            intpol_get_gp_s = 0.0

      END SELECT

      CALL profiler_set_stop_time('intpol_get_gp_s', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a inteferometry signal and a
!>  cartesian position.
!>
!>  Calculates the guassian process kernel between the signal and the position.
!>  Only valid for interferometry types. This is the second signal.
!>
!>  @param[in] this         A @ref intpol_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] x_cart       The cartesian position of to get the kernel at.
!>  @param[in] flags        State flags to send to the kernel.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_gp_x(this, a_model, x_cart, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: intpol_get_gp_x
      TYPE (intpol_class), INTENT(in)        :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: flags

!  local variables
      REAL (rprec)                           :: start_time
      CHARACTER(len=1), ALLOCATABLE          :: context(:)
      INTEGER                                :: context_length
      TYPE (intpol_gp_context_x)             :: gp_context

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context%model => a_model
      gp_context%xcart = x_cart


!  Cast model into a data context. This is the equivalent to casting to a void
!  pointer in C.
      context_length = SIZE(TRANSFER(gp_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(gp_context, context)

      SELECT CASE (this%type)

         CASE (intpol_interferometry_type)
            intpol_get_gp_x = path_integrate(a_model%int_params,               &
     &                                       this%chord_path,                  &
     &                                       gp_function_x, context)

         CASE (intpol_polarimetry_type)
            intpol_get_gp_x = intpol_polar_constant                            &
     &                      * (this%wavelength**2.0_dp)                        &
     &                      * path_integrate(a_model%int_params,               &
     &                                       this%chord_path,                  &
     &                                       gp_function_pol_x,                &
     &                                       context)

            IF (this%in_degrees) THEN
               intpol_get_gp_x = intpol_get_gp_x/degree
            END IF

         CASE DEFAULT
            intpol_get_gp_x = 0.0

      END SELECT

      DEALLOCATE(context)

      CALL profiler_set_stop_time('intpol_get_gp_x', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the integration element for the signal.
!>
!>  Integrated signals may contain extra quanties that are integrated. In the on
!>  context of gaussian process, the correlation of two signals measureing
!>  different quanties but using the same kernel can be computed. This ensures
!>  the one of the signals has the correct value.
!>
!>  @param[in] this         A @ref intpol_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_dx(this, a_model, xcart, dxcart, length, dx,         &
     &                       flags)
      REAL (rprec)                           :: intpol_get_dx
      TYPE (intpol_class), INTENT(in)        :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      INTEGER, INTENT(in)                    :: flags

! local variables
      REAL (rprec), DIMENSION(3)             :: bcart
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (intpol_interferometry_type)
            intpol_get_dx = dx

         CASE (intpol_polarimetry_type)
            bcart = equilibrium_get_B_vec(a_model%equilibrium, xcart,          &
     &                                    .false.)
            intpol_get_dx = DOT_PRODUCT(bcart, dxcart)

         CASE DEFAULT
            intpol_get_dx = 0.0

      END SELECT

      CALL profiler_set_stop_time('intpol_get_dx', start_time)

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Interferometer callback function.
!>
!>  Returns the value of the denisty times the change in path length. This
!>  function is passed to @ref integration_path::path_integrate to act as a
!>  callback. The denisty is provided by @ref model::model_get_ne.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref model.
!>  @param[in] xcart   A integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration. Not used in this function.
!>  @param[in] dx      A scalar change in path.
!>  @returns The ne(x)*dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION int_function(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: int_function

!  local variables
      TYPE (model_class)                     :: a_model
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      a_model = TRANSFER(context, a_model)
      int_function = model_get_ne(a_model, xcart)*dx

      CALL profiler_set_stop_time('int_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Polarmetry callback function.
!>
!>  Returns the value of the denisty times the dot product of the magnetic
!>  field vector and the path direction. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The denisty is
!>  provided by @ref model::model_get_ne. The magnetic field vector is proved by
!>  @ref equilibrium::equilibrium_get_B_vec
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref model.
!>  @param[in] xcart   A integration point.
!>  @param[in] dxcart  A vector change in path.
!>  @param[in] length  Length along the integration. Not used in this function.
!>  @param[in] dx      A scalar change in path. Not used in this function.
!>  @returns The ne(x)*B.dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION pol_function(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx
      REAL (rprec)                           :: pol_function

!  local variables
      TYPE (model_class)                     :: a_model
      REAL (rprec), DIMENSION(3)             :: bcart
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      a_model = TRANSFER(context, a_model)

      bcart = equilibrium_get_B_vec(a_model%equilibrium, xcart, .false.)
      pol_function =                                                           &
     &   model_get_ne(a_model, xcart)*DOT_PRODUCT(bcart, dxcart)

      CALL profiler_set_stop_time('pol_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Electron density gaussian process callback function for signal point
!>  kernel evaluation.
!>
!>  Returns the value of the density guassian process kernel times the change in
!>  path length. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The density
!>  kernel is provided by @ref model::model_get_gp_ne. This is the second
!>  signal.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref intpol_gp_context_i for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The k(x, i)*dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION gp_function_i(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: gp_function_i
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (intpol_gp_context_i)             :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context = TRANSFER(context, gp_context)
      gp_function_i = model_get_gp_ne(gp_context%model, xcart,                 &
     &                                gp_context%i)*dx

      CALL profiler_set_stop_time('gp_function_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Electron density gaussian process callback function for signal point
!>  kernel evaluation for a polarimetry signal.
!>
!>  Returns the value of the density guassian process kernel times the change in
!>  path length. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The density
!>  kernel is provided by @ref model::model_get_gp_ne. This is the second
!>  signal.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref intpol_gp_context_i for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The k(x, i)*dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION gp_function_pol_i(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: gp_function_pol_i
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (intpol_gp_context_i)             :: gp_context
      REAL (rprec), DIMENSION(3)             :: bcart
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context = TRANSFER(context, gp_context)
      bcart = equilibrium_get_B_vec(gp_context%model%equilibrium, xcart,       &
     &                              .false.)
      gp_function_pol_i = model_get_gp_ne(gp_context%model, xcart,             &
     &                                    gp_context%i)                        *
     &                  * DOT_PRODUCT(bcart, dxcart)

      CALL profiler_set_stop_time('gp_function_pol_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Electron density gaussian process callback function for kernel
!>  evaluation of two signals.
!>
!>  Calculates the guassian process kernel between the signal and a signal.
!>  Calls back to the @ref signal module to call the other signal. This is the
!>  second signal.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref intpol_gp_context_x for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The k(x, y)*dl at point y.
!-------------------------------------------------------------------------------
      FUNCTION gp_function_x(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: gp_function_x
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (intpol_gp_context_x)             :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  This is the second signal so put xcart in the second position.
      gp_context = TRANSFER(context, gp_context)
      gp_function_x = model_get_gp_ne(gp_context%model,                        &
     &                                xcart, gp_context%xcart)*dx

      CALL profiler_set_stop_time('gp_function_x', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Electron density gaussian process callback function for kernel
!>  evaluation of two signals.
!>
!>  Calculates the guassian process kernel between the signal and a signal.
!>  Calls back to the @ref signal module to call the other signal. This is the
!>  second signal for a polarimetry signal.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref intpol_gp_context_x for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The k(x, y)*dl at point y.
!-------------------------------------------------------------------------------
      FUNCTION gp_function_pol_x(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: gp_function_pol_x
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (intpol_gp_context_x)             :: gp_context
      REAL (rprec), DIMENSION(3)             :: bcart
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  This is the second signal so put xcart in the second position.
      gp_context = TRANSFER(context, gp_context)
      bcart = equilibrium_get_B_vec(gp_context%model%equilibrium, xcart,       &
     &                              .false.)
      gp_function_pol_x = model_get_gp_ne(gp_context%model,                    &
     &                                    xcart, gp_context%xcart)             &
     &                  * DOT_PRODUCT(bcart, dxcart)

      CALL profiler_set_stop_time('gp_function_x', start_time)

      END FUNCTION

      END MODULE
