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
!  sxrem module parameters
!*******************************************************************************
!>  Type descriptor for sxrem no type.
      INTEGER, PARAMETER :: sxrem_no_type    = -1
!>  Type descriptor for sxrem type represented a standard emission diagnostic.
      INTEGER, PARAMETER :: sxrem_emiss_type = 0
!>  Type descriptor for sxrem type represented a te diagnostic.
      INTEGER, PARAMETER :: sxrem_ti_type    = 1

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
!>  Type descirptor of the intpol type.
!>  @par Possible values are:
!>  * @ref sxrem_no_type
!>  * @ref sxrem_emiss_type
!>  * @ref sxrem_ti_type
         INTEGER                :: type
!>  The complete path of the chord.
         TYPE (vertex), POINTER :: chord_path => null()
!>  Geometric chord factor.
         REAL (rprec)           :: geo
!>  Index of the soft x-ray emissivity profile for this chord.
         INTEGER                :: profile_number
      END TYPE sxrem_class

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the callback function.
!-------------------------------------------------------------------------------
      TYPE sxrem_context
!>  The index of the emissivity profile model.
         INTEGER                     :: profile_number
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the guassian process
!>  callback function for point.
!-------------------------------------------------------------------------------
      TYPE sxrem_gp_context_i
!>  The index of the emissivity profile model.
         INTEGER                     :: profile_number
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
!>  Position index.
         INTEGER                     :: i
!>  Gaussian process kernel flags.
         INTEGER                     :: flags = model_state_all_off
      END TYPE

!-------------------------------------------------------------------------------
!>  Structure to hold all memory needed to be sent to the guassian process
!>  callback function for signal.
!-------------------------------------------------------------------------------
      TYPE sxrem_gp_context_x
!>  The index of the emissivity profile model.
         INTEGER                     :: profile_number
!>  Reference to a @ref model::model_class object.
         TYPE (model_class), POINTER :: model => null()
!>  First position.
         REAL (rprec), DIMENSION(3)  :: xcart
!>  Gaussian process kernel flags.
         INTEGER                     :: flags
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref sxrem_class types using
!>  @ref intpol_construct_int or @ref intpol_construct_pol.
!-------------------------------------------------------------------------------
      INTERFACE sxrem_construct
         MODULE PROCEDURE sxrem_construct_emiss,                                &
     &                    sxrem_construct_ti
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to get the guassian process kernel values.
!-------------------------------------------------------------------------------
      INTERFACE sxrem_get_gp
         MODULE PROCEDURE sxrem_get_gp_i,                                      &
     &                    sxrem_get_gp_s,                                      &
     &                    sxrem_get_gp_x
      END INTERFACE

      PRIVATE :: sxr_function, gp_function_i, gp_function_x

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_class object for emission.
!>
!>  Allocates memory and initializes a @ref sxrem_class object.
!>
!>  @param[in] start_path     Starting point of a sxrem chord.
!>  @param[in] end_path       Ending point of a sxrem chord.
!>  @param[in] geo            Geometric factor of the chord.
!>  @param[in] profile_number Index of the soft x-ray emissivity profile.
!>  @returns A pointer to a constructed @ref sxrem_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct_emiss(start_path, end_path, geo,                &
     &                               profile_number)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER            :: sxrem_construct_emiss
      REAL (rprec), DIMENSION(3), INTENT(in) :: start_path
      REAL (rprec), DIMENSION(3), INTENT(in) :: end_path
      REAL (rprec), INTENT(in)               :: geo
      INTEGER, INTENT(in)                    :: profile_number

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(sxrem_construct_emiss)

      sxrem_construct_emiss%type = sxrem_emiss_type

      CALL path_append_vertex(sxrem_construct_emiss%chord_path,                &
     &                        start_path)
      CALL path_append_vertex(sxrem_construct_emiss%chord_path,                &
     &                        end_path)

      sxrem_construct_emiss%geo = geo
      sxrem_construct_emiss%profile_number = profile_number

      CALL profiler_set_stop_time('sxrem_construct_emiss', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref sxrem_class object for ti.
!>
!>  Allocates memory and initializes a @ref sxrem_class object.
!>
!>  @param[in] start_path     Starting point of a sxrem chord.
!>  @param[in] end_path       Ending point of a sxrem chord.
!>  @param[in] profile_number Index of the soft x-ray emissivity profile.
!>  @returns A pointer to a constructed @ref sxrem_class object.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct_ti(start_path, end_path, profile_number)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (sxrem_class), POINTER            :: sxrem_construct_ti
      REAL (rprec), DIMENSION(3), INTENT(in) :: start_path
      REAL (rprec), DIMENSION(3), INTENT(in) :: end_path
      INTEGER, INTENT(in)                    :: profile_number

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(sxrem_construct_ti)

      sxrem_construct_ti%type = sxrem_ti_type

      CALL path_append_vertex(sxrem_construct_ti%chord_path,                   &
     &                        start_path)
      CALL path_append_vertex(sxrem_construct_ti%chord_path, end_path)

      sxrem_construct_ti%geo = 1.0
      sxrem_construct_ti%profile_number = profile_number

      CALL profiler_set_stop_time('sxrem_construct_ti', start_time)

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
!>  path. The signal is computed based on the type.
!>
!>  @param[in]  this         A @ref sxrem_class instance.
!>  @param[in]  a_model      A @ref model instance.
!>  @param[out] sigma        The modeled sigma.
!>  @param[in]  last_value   Last good value in case the signal did not change.
!>  @param[in]  scale_factor Factor to scale the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_signal(this, a_model, sigma,                  &
     &                                  last_value, scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_signal
      TYPE (sxrem_class), INTENT(in)          :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor

! local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (sxrem_context)                    :: sxr_context
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE (this%type)

         CASE (sxrem_emiss_type)
            sxrem_get_modeled_signal =                                         &
     &         sxrem_get_modeled_emiss_signal(this, a_model, sigma,            &
     &                                        last_value, scale_factor)

         CASE (sxrem_ti_type)
            sxrem_get_modeled_signal =                                         &
     &         sxrem_get_modeled_ti_signal(this, a_model, sigma,               &
     &                                     last_value, scale_factor)

         CASE DEFAULT
            sxrem_get_modeled_signal = 0.0
            sigma = 0.0

      END SELECT

      CALL profiler_set_stop_time('sxrem_get_modeled_signal',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled emiss signal.
!>
!>  Calculates the soft x-ray emissivity signal by integrating along the chord
!>  path. The emissivity is provided by the @ref sxr_function.
!>
!>  @param[in]  this         A @ref sxrem_class instance.
!>  @param[in]  a_model      A @ref model instance.
!>  @param[out] sigma        The modeled sigma.
!>  @param[in]  last_value   Last good value in case the signal did not change.
!>  @param[in]  scale_factor Factor to scale the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_emiss_signal(this, a_model, sigma,            &
     &                                        last_value, scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_emiss_signal
      TYPE (sxrem_class), INTENT(in)          :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor

! local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (sxrem_context)                    :: sxr_context
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_get_modeled_emiss_signal = last_value
      sigma = 0.0

      IF (BTEST(a_model%state_flags, model_state_vmec_flag)     .or.           &
     &    BTEST(a_model%state_flags, model_state_siesta_flag)   .or.           &
     &    BTEST(a_model%state_flags, model_state_shift_flag)    .or.           &
     &    BTEST(a_model%state_flags, model_state_sxrem_flag +                  &
     &                               (this%profile_number - 1)) .or.           &
     &    BTEST(a_model%state_flags, model_state_signal_flag)) THEN

!  The relevant data for the soft x-ray context.
         sxr_context%profile_number = this%profile_number
         sxr_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
         context_length = SIZE(TRANSFER(sxr_context, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(sxr_context, context)

         sxrem_get_modeled_emiss_signal(1) =                                   &
     &      path_integrate(a_model%int_params, this%chord_path,                &
     &                     sxr_function, context)*                             &
     &         this%geo*scale_factor

         DEALLOCATE(context)
      END IF

      CALL profiler_set_stop_time('sxrem_get_modeled_emiss_signal',            &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled ion temperature signal.
!>
!>  Calculates the soft x-ray emissivity/ion temperature signal by integrating
!>  along the chord path. The emissivity is provided by the @ref ti_function.
!>
!>  @param[in]  this         A @ref sxrem_class instance.
!>  @param[in]  a_model      A @ref model instance.
!>  @param[out] sigma        The modeled sigma.
!>  @param[in]  last_value   Last good value in case the signal did not change.
!>  @param[in]  scale_factor Factor to scale the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_modeled_ti_signal(this, a_model, sigma,               &
     &                                     last_value, scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: sxrem_get_modeled_ti_signal
      TYPE (sxrem_class), INTENT(in)          :: this
      TYPE (model_class), TARGET, INTENT(in)  :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor

! local variables
      CHARACTER(len=1), ALLOCATABLE           :: context(:)
      INTEGER                                 :: context_length
      TYPE (sxrem_context)                    :: sxr_context
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxrem_get_modeled_ti_signal = last_value
      sigma = 0.0

      IF (BTEST(a_model%state_flags, model_state_vmec_flag)     .or.           &
     &    BTEST(a_model%state_flags, model_state_siesta_flag)   .or.           &
     &    BTEST(a_model%state_flags, model_state_shift_flag)    .or.           &
     &    BTEST(a_model%state_flags, model_state_ti_flag)       .or.           &
     &    BTEST(a_model%state_flags, model_state_sxrem_flag +                  &
     &                               (this%profile_number - 1)) .or.           &
     &    BTEST(a_model%state_flags, model_state_signal_flag)) THEN

!  The relevant data for the soft x-ray context.
         sxr_context%profile_number = this%profile_number
         sxr_context%model => a_model

!  Cast model into a data to a context. This is the equivalent to casting to a
!  void pointer in C.
         context_length = SIZE(TRANSFER(sxr_context, context))
         ALLOCATE(context(context_length))
         context = TRANSFER(sxr_context, context)

         sxrem_get_modeled_ti_signal(1) =                                      &
     &      path_integrate(a_model%int_params, this%chord_path,                &
     &                     ti_function, context)*scale_factor

         DEALLOCATE(context)
      END IF

      CALL profiler_set_stop_time('sxrem_get_modeled_ti_signal',               &
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

      IF (this%type .eq. sxrem_emiss_type) THEN
         WRITE (sxrem_get_signal_type, 1000) this%profile_number
      ELSE
         WRITE (sxrem_get_signal_type, 1001) this%profile_number
      END IF

      CALL profiler_set_stop_time('sxrem_get_signal_type', start_time)

1000  FORMAT('sxrch(',i2,')')
1001  FORMAT('sxrch_ti(',i2,')')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a sxrem signal and a position.
!>
!>  Calculates the guassian process kernel between the signal and the position.
!>  Soft x-ray emission kernels are provided by
!>  @ref model::model_get_gp_sxrem.
!>
!>  @param[in] this         A @ref sxrem_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] i            Index of the position for the kernel.
!>  @param[in] flags        State flags to send to the kernel.
!>  @returns Kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_gp_i(this, a_model, i, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: sxrem_get_gp_i
      TYPE (sxrem_class), INTENT(in)         :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      INTEGER, INTENT(in)                    :: i

!  local variables
      REAL (rprec)                           :: start_time
      CHARACTER(len=1), ALLOCATABLE          :: context(:)
      INTEGER                                :: context_length
      TYPE (sxrem_gp_context_i)              :: gp_context
      INTEGER, INTENT(in)                    :: flags

!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the guassian process context.
      gp_context%profile_number = this%profile_number
      gp_context%model => a_model
      gp_context%i = i
      gp_context%flags = flags

!  Cast model into a data context. This is the equivalent to casting to a void
!  pointer in C.
      context_length = SIZE(TRANSFER(gp_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(gp_context, context)

      SELECT CASE (this%type)

         CASE (sxrem_emiss_type)
            sxrem_get_gp_i = path_integrate(a_model%int_params,                &
     &                                      this%chord_path,                   &
     &                                      gp_function_i, context)            &
     &                     * this%geo

         CASE (sxrem_ti_type)
            sxrem_get_gp_i = path_integrate(a_model%int_params,                &
     &                                      this%chord_path,                   &
     &                                      gp_function_i, context)

         CASE DEFAULT
            sxrem_get_gp_i = 0.0

      END SELECT

      DEALLOCATE(context)

      CALL profiler_set_stop_time('sxrem_get_gp_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a sxrem signal and a signal.
!>
!>  Calculates the guassian process kernel between the signal and a signal.
!>  Calls back to the @ref signal module to call the other signal. This is the
!>  first signal.
!>
!>  @param[in] this         A @ref sxrem_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] gp_context   A @ref signal_class instance.
!>  @param     gp_callback  Function pointer to call back to signal module.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_gp_s(this, a_model, context, callback)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: sxrem_get_gp_s
      TYPE (sxrem_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model
      CHARACTER (len=1), INTENT(in)  :: context(:)
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
      REAL (rprec)                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the guassian process context.
      sxrem_get_gp_s = path_integrate(a_model%int_params,                      &
     &                                this%chord_path, callback,               &
     &                                context)

      IF (this%type .eq. sxrem_emiss_type) THEN
         sxrem_get_gp_s = sxrem_get_gp_s*this%geo
      END IF

      CALL profiler_set_stop_time('sxrem_get_gp_s', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a sxrem signal and a cartesian
!>  position.
!>
!>  Calculates the guassian process kernel between the signal and the position.
!>  This is the second signal.
!>
!>  @param[in] this         A @ref sxrem_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] x_cart       The cartesian position of to get the kernel at.
!>  @param[in] flags        State flags to send to the kernel.
!>  @returns Kernel value for the signal and the signal.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_gp_x(this, a_model, x_cart, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: sxrem_get_gp_x
      TYPE (sxrem_class), INTENT(in)         :: this
      TYPE (model_class), TARGET, INTENT(in) :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      INTEGER, INTENT(in)                    :: flags

!  local variables
      REAL (rprec)                           :: start_time
      CHARACTER(len=1), ALLOCATABLE          :: context(:)
      INTEGER                                :: context_length
      TYPE (sxrem_gp_context_x)              :: gp_context

!  Start of executable code
      start_time = profiler_get_start_time()

!  The relevant data for the guassian process context.
      gp_context%profile_number = this%profile_number
      gp_context%model => a_model
      gp_context%xcart = x_cart
      gp_context%flags = flags

!  Cast model into a data context. This is the equivalent to casting to a void
!  pointer in C.
      context_length = SIZE(TRANSFER(gp_context, context))
      ALLOCATE(context(context_length))
      context = TRANSFER(gp_context, context)

      sxrem_get_gp_x =  path_integrate(a_model%int_params,                     &
     &                                 this%chord_path, gp_function_x,         &
     &                                 context)

      IF (this%type .eq. sxrem_emiss_type) THEN
         sxrem_get_gp_x = sxrem_get_gp_x*this%geo
      END IF

      DEALLOCATE(context)

      CALL profiler_set_stop_time('sxrem_get_gp_x', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the integration element for the signal.
!>
!>  Integrated signals may contain extra quanties that are integrated. In the on
!>  context of gaussian process, the correlation of two signals measureing
!>  different quanties but using the same kernel can be computed. This ensures
!>  the one of the signals has the correct value. For XICS signals, dx can
!>  change depending on if the emission or the ion temperature is the current
!>  guassian process.
!>
!>  @param[in] this         A @ref intpol_class instance.
!>  @param[in] a_model      A @ref model instance.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @param[in] flags   State flags to send to the kernel.
!>  @returns dx value for the signal.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_get_dx(this, a_model, xcart, dxcart, length, dx,          &
     &                      flags)
      REAL (rprec)                           :: sxrem_get_dx
      TYPE (sxrem_class), INTENT(in)         :: this
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

         CASE (sxrem_emiss_type)
            sxrem_get_dx = dx

         CASE (sxrem_ti_type)
            IF (BTEST(flags, model_state_sxrem_flag +                          &
     &                       (this%profile_number - 1))) THEN
               sxrem_get_dx = model_get_ti(a_model, xcart)*dx
            ELSE IF (BTEST(flags, model_state_ti_flag)) THEN
               sxrem_get_dx = model_get_sxrem(a_model, xcart,                  &
     &                                        this%profile_number)*dx
            ELSE
               sxrem_get_dx = 0.0
            END IF

         CASE DEFAULT
            sxrem_get_dx = 0.0

      END SELECT

      CALL profiler_set_stop_time('sxrem_get_dx', start_time)

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

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: sxr_function
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

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

!-------------------------------------------------------------------------------
!>  @brief Ion temperature callback function.
!>
!>  Returns the value of the ion temperature times the soft x-ray emissivity
!>  times the change in path length. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The soft x-ray
!>  emissivity is provided by @ref model::model_get_sxrem and the ion
!>  temperature is provided by @ref model::model_get_ti.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref sxrem_context for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The te(x)*dl at point x.
!-------------------------------------------------------------------------------
      FUNCTION ti_function(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: ti_function
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (sxrem_context)                   :: sxr_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      sxr_context = TRANSFER(context, sxr_context)
      ti_function = model_get_sxrem(sxr_context%model, xcart,                  &
     &                              sxr_context%profile_number)                &
     &            * model_get_ti(sxr_context%model, xcart)*dx

      CALL profiler_set_stop_time('ti_function', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Soft x-ray gaussian process callback function for signal point
!>  kernel evaluation.
!>
!>  Returns the value of the soft x-ray guassian process kernel times the change
!>  in path length. This function is passed to
!>  @ref integration_path::path_integrate to act as a callback. The soft x-ray
!>  kernel is provided by @ref model::model_get_gp_sxrem.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref sxrem_gp_context_i for the model.
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
      TYPE (sxrem_gp_context_i)              :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      gp_context = TRANSFER(context, gp_context)
      IF (BTEST(gp_context%flags, model_state_sxrem_flag +                     &
     &                            (gp_context%profile_number - 1))) THEN
         gp_function_i = model_get_gp_sxrem(gp_context%model, xcart,           &
     &                                      gp_context%i,                      &
     &                                      gp_context%profile_number)
      ELSE IF (BTEST(gp_context%flags, model_state_ti_flag)) THEN
         gp_function_i = model_get_gp_ti(gp_context%model, xcart,              &
     &                                   gp_context%i)
      ELSE
         gp_function_i = 0.0
      END IF

      gp_function_i = gp_function_i*dx

      CALL profiler_set_stop_time('gp_function_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Soft x-ray gaussian process callback function for kernel evaluation
!>  of two signals.
!>
!>  Calculates the guassian process kernel between the signal and a signal.
!>  Calls back to the @ref signal module to call the other signal. This is the
!>  second signal so xcart should go in the first index.
!>
!>  @see integration_path
!>
!>  @param[in] context A @ref sxrem_gp_context_x for the model.
!>  @param[in] xcart   The integration point.
!>  @param[in] dxcart  A vector change in path. Not used in this function.
!>  @param[in] length  Length along the integration.
!>  @param[in] dx      A scalar change in path.
!>  @returns The k(x, y)*dl at point y.
!-------------------------------------------------------------------------------
      FUNCTION gp_function_x(context, xcart, dxcart, length, dx)
      USE v3_utilities, ONLY: assert

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: gp_function_x
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

! local variables
      TYPE (sxrem_gp_context_x)              :: gp_context
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  This is the second signal so put xcart in the second position.
      gp_context = TRANSFER(context, gp_context)
      IF (BTEST(gp_context%flags, model_state_sxrem_flag +                     &
     &                            (gp_context%profile_number - 1))) THEN
         gp_function_x = model_get_gp_sxrem(gp_context%model,                  &
     &                                      xcart, gp_context%xcart,           &
     &                                      gp_context%profile_number)
      ELSE IF (BTEST(gp_context%flags, model_state_ti_flag)) THEN
         gp_function_x = model_get_gp_ti(gp_context%model, xcart,              &
     &                                   gp_context%xcart)
      ELSE
         gp_function_x = 0.0
      END IF

      gp_function_x = gp_function_x*dx

      CALL profiler_set_stop_time('gp_function_x', start_time)

      END FUNCTION

      END MODULE
