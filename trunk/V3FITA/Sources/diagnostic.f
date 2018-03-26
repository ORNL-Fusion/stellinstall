!*******************************************************************************
!>  @file diagnostic.f
!>  @brief Contains module @ref diagnostic
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref diagnostic_class.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref magnetic, @ref sxrem, @ref intpol, @ref thomson, @ref extcurz
!*******************************************************************************

      MODULE diagnostic

      USE magnetic
      USE sxrem
      USE intpol
      USE thomson
      USE extcurz
      USE mse
      USE ece

      IMPLICIT NONE

!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
!>  Type descriptor for the @ref diagnostic_class no sub class.
      INTEGER, PARAMETER :: diagnostic_no_type       = -1
!>  Type descriptor for the @ref diagnostic_class subclass @ref magnetic
      INTEGER, PARAMETER :: diagnostic_magnetic_type = 0
!>  Type descriptor for the @ref diagnostic_class subclass @ref sxrem
      INTEGER, PARAMETER :: diagnostic_sxrem_type    = 1
!>  Type descriptor for the @ref diagnostic_class subclass @ref intpol
      INTEGER, PARAMETER :: diagnostic_intpol_type   = 2
!>  Type descriptor for the @ref diagnostic_class subclass @ref thomson
      INTEGER, PARAMETER :: diagnostic_thomson_type  = 3
!>  Type descriptor for the @ref diagnostic_class subclass @ref extcurz
      INTEGER, PARAMETER :: diagnostic_extcurz_type  = 4
!>  Type descriptor for the @ref diagnostic_class subclass @ref mse
      INTEGER, PARAMETER :: diagnostic_mse_type      = 5
!>  Type descriptor for the @ref diagnostic_class subclass @ref ece
      INTEGER, PARAMETER :: diagnostic_ece_type      = 6

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) diagnostic base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a diagnostic signal.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref magnetic, @ref sxrem, @ref intpol, @ref thomson, @ref extcurz
!-------------------------------------------------------------------------------
      TYPE diagnostic_class
!>  Type descirptor of the diagnostic subclass type.
!>  @par Possible values are:
!>  * @ref diagnostic_no_type
!>  * @ref diagnostic_magnetic_type,
!>  * @ref diagnostic_sxrem_type
!>  * @ref diagnostic_intpol_type,
!>  * @ref diagnostic_thomson_type
!>  * @ref diagnostic_extcurz_type
!>  * @ref diagnostic_mse_type
!>  * @ref diagnostic_ece_type
         INTEGER                       :: type = diagnostic_no_type
!>  An instance of a @ref magnetic object. 
         TYPE (magnetic_class), POINTER :: magnetic => null()
!>  An instance of a @ref sxrem object.
         TYPE (sxrem_class), POINTER    :: sxrem => null()
!>  An instance of a @ref intpol object.
         TYPE (intpol_class), POINTER   :: intpol => null()
!>  An instance of a @ref thomson object.
         TYPE (thomson_class), POINTER  :: thomson => null()
!>  An instance of a @ref extcurz object.
         TYPE (extcurz_class), POINTER  :: extcurz => null()
!>  An instance of a @ref mse object.
         TYPE (mse_class), POINTER      :: mse => null()
!>  An instance of a @ref mse object.
         TYPE (ece_class), POINTER      :: ece => null()
      END TYPE diagnostic_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref diagnostic_class types using
!>  @ref diagnostic_construct_magnetic, @ref diagnostic_construct_sxrem,
!>  @ref diagnostic_construct_intpol, @ref diagnostic_construct_thomson,
!>  @ref diagnostic_construct_extcurz, @ref diagnostic_construct_mse or
!>  @ref diagnostic_construct_ece.
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_construct
         MODULE PROCEDURE diagnostic_construct_magnetic,                       &
     &                    diagnostic_construct_sxrem,                          &
     &                    diagnostic_construct_intpol,                         &
     &                    diagnostic_construct_thomson,                        &
     &                    diagnostic_construct_extcurz,                        &
     &                    diagnostic_construct_mse,                            &
     &                    diagnostic_construct_ece
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to get the guassian process kernel values.
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_get_gp
         MODULE PROCEDURE diagnostic_get_gp_i,                                 &
     &                    diagnostic_get_gp_s,                                 &
     &                    diagnostic_get_gp_x
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref magnetic object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] magnetic_object An instance of a @ref magnetic subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_magnetic(magnetic_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_magnetic
      TYPE (magnetic_class), POINTER   :: magnetic_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_magnetic)

      diagnostic_construct_magnetic%type = diagnostic_magnetic_type
      diagnostic_construct_magnetic%magnetic => magnetic_object

      CALL profiler_set_stop_time('diagnostic_construct_magnetic',             &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref sxrem object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] sxrem_object An instance of a @ref sxrem subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_sxrem(sxrem_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_sxrem
      TYPE (sxrem_class), POINTER      :: sxrem_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_sxrem)

      diagnostic_construct_sxrem%type = diagnostic_sxrem_type
      diagnostic_construct_sxrem%sxrem => sxrem_object

      CALL profiler_set_stop_time('diagnostic_construct_sxrem',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref intpol object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] intpol_object An instance of a @ref intpol subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_intpol(intpol_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_intpol
      TYPE (intpol_class), POINTER     :: intpol_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_intpol)

      diagnostic_construct_intpol%type = diagnostic_intpol_type
      diagnostic_construct_intpol%intpol => intpol_object

      CALL profiler_set_stop_time('diagnostic_construct_intpol',               &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref thomson object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] thomson_object An instance of a @ref thomson subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_thomson(thomson_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_thomson
      TYPE (thomson_class), POINTER    :: thomson_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_thomson)

      diagnostic_construct_thomson%type = diagnostic_thomson_type
      diagnostic_construct_thomson%thomson => thomson_object

      CALL profiler_set_stop_time('diagnostic_construct_thomson',              &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref extcurz object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] extcurz_object An instance of a @ref extcurz subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_extcurz(extcurz_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_extcurz
      TYPE (extcurz_class), POINTER    :: extcurz_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_extcurz)

      diagnostic_construct_extcurz%type = diagnostic_extcurz_type
      diagnostic_construct_extcurz%extcurz => extcurz_object

      CALL profiler_set_stop_time('diagnostic_construct_extcurz',              &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref mse object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] mse_object An instance of a @ref mse subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_mse(mse_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_mse
      TYPE (mse_class), POINTER        :: mse_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_mse)

      diagnostic_construct_mse%type = diagnostic_mse_type
      diagnostic_construct_mse%mse => mse_object

      CALL profiler_set_stop_time('diagnostic_construct_mse',                  &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref diagnostic_class containing a @ref ece object.
!>
!>  Allocates memory and initializes a @ref diagnostic_class object.
!>
!>  @param[in] ece_object An instance of a @ref ece subclass.
!>  @returns A pointer to a constructed @ref diagnostic_class object.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_ece(ece_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_ece
      TYPE (ece_class), POINTER        :: ece_object

!  local variables
      REAL (rprec)                     :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(diagnostic_construct_ece)

      diagnostic_construct_ece%type = diagnostic_ece_type
      diagnostic_construct_ece%ece => ece_object

      CALL profiler_set_stop_time('diagnostic_construct_ece',                  &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref diagnostic_class object.
!>
!>  Deallocates memory and uninitializes a @ref diagnostic_class object.
!>
!>  @param[inout] this A @ref diagnostic_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER :: this

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

      this%type = diagnostic_no_type

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
!>  @ref diagnostic_class was constructed with.
!>
!>  @see magnetic::magnetic_get_modeled_signal
!>  @see sxrem::sxrem_get_modeled_signal
!>  @see intpol::intpol_get_modeled_signal
!>  @see thomson::thomson_get_modeled_signal
!>  @see extcurz::extcurz_get_modeled_signal
!>  @see mse::mse_get_modeled_signal
!>  @see ece::ece_get_modeled_signal
!>
!>  @param[in]  this         A @ref diagnostic_class instance.
!>  @param[in]  a_model      A @ref model instance.
!>  @param[out] sigma        The modeled sigma.
!>  @param[in]  last_value   Last good value in case the signal did not change.
!>  @param[in]  scale_factor Factor to scale the modeled signal.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_modeled_signal(this, a_model, sigma,             &
     &                                       last_value, scale_factor)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: diagnostic_get_modeled_signal
      TYPE (diagnostic_class), INTENT(in)     :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma
      REAL (rprec), DIMENSION(4), INTENT(in)  :: last_value
      REAL (rprec), INTENT(in)                :: scale_factor

!  local variables
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)
         CASE (diagnostic_magnetic_type)
            diagnostic_get_modeled_signal =                                    &
     &         magnetic_get_modeled_signal(this%magnetic, a_model,             &
     &                                     sigma, last_value,                  &
     &                                     scale_factor)

         CASE (diagnostic_sxrem_type)
            diagnostic_get_modeled_signal =                                    &
     &         sxrem_get_modeled_signal(this%sxrem, a_model, sigma,            &
     &                                  last_value, scale_factor)

         CASE (diagnostic_intpol_type)
            diagnostic_get_modeled_signal =                                    &
     &         intpol_get_modeled_signal(this%intpol, a_model, sigma,          &
     &                                   last_value, scale_factor)

         CASE (diagnostic_thomson_type)
            diagnostic_get_modeled_signal =                                    &
     &         thomson_get_modeled_signal(this%thomson, a_model, sigma,        &
     &                                    last_value, scale_factor)

         CASE (diagnostic_extcurz_type)
            diagnostic_get_modeled_signal =                                    &
     &         extcurz_get_modeled_signal(this%extcurz, a_model, sigma,        &
     &                                    last_value, scale_factor)

         CASE (diagnostic_mse_type)
            diagnostic_get_modeled_signal =                                    &
     &         mse_get_modeled_signal(this%mse, a_model, sigma,                &
     &                                last_value, scale_factor)

         CASE (diagnostic_ece_type)
            diagnostic_get_modeled_signal =                                    &
     &         ece_get_modeled_signal(this%ece, a_model, sigma,                &
     &                                last_value, scale_factor)

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_modeled_signal',             &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the diagnostic type.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_signal_type method of the subclass instance
!>  @ref diagnostic_class was constructed with.
!>
!>  @see magnetic::magnetic_get_signal_type
!>  @see sxrem::sxrem_get_signal_type
!>  @see intpol::intpol_get_signal_type
!>  @see thomson::thomson_get_signal_type
!>  @see extcurz::extcurz_get_signal_type
!>  @see mse::mse_get_signal_type
!>  @see ece::ece_get_signal_type
!>
!>  @param[in] this A @ref diagnostic_class instance.
!>  @returns A string describing the diagnostic type.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: diagnostic_get_signal_type
      TYPE (diagnostic_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (diagnostic_magnetic_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         magnetic_get_signal_type(this%magnetic))

         CASE (diagnostic_sxrem_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         sxrem_get_signal_type(this%sxrem))

         CASE (diagnostic_intpol_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         intpol_get_signal_type(this%intpol))

         CASE (diagnostic_thomson_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         thomson_get_signal_type(this%thomson))

         CASE (diagnostic_extcurz_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         extcurz_get_signal_type(this%extcurz))

         CASE (diagnostic_mse_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         mse_get_signal_type(this%mse))

         CASE (diagnostic_ece_type)
            diagnostic_get_signal_type = TRIM('diagnostic ' //                 &
     &         ece_get_signal_type(this%ece))

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_signal_type',                &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  This method is meant to be overwritten by a subclass method if need be. As a
!>  result, returns the get_header method of the subclass instance
!>  @ref diagnostic_class was constructed with. Otherwise uses the default.
!>
!>  @see magnetic::magnetic_get_header
!>
!>  @param[in] this A @ref diagnostic_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) ::                        &
     &   diagnostic_get_header
      TYPE (diagnostic_class), INTENT(in)            :: this

!  local variables
      REAL (rprec)                                   :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (diagnostic_magnetic_type)
            diagnostic_get_header = magnetic_get_header(this%magnetic)

         CASE (diagnostic_mse_type)
            diagnostic_get_header = mse_get_header(this%mse)

         CASE (diagnostic_ece_type)
            diagnostic_get_header = ece_get_header(this%ece)

         CASE DEFAULT
            diagnostic_get_header(1:3) = 'N/A'
            diagnostic_get_header(4) = 'model_sig(1)'
            diagnostic_get_header(5) = 'model_sig(2)'
            diagnostic_get_header(6) = 'model_sig(3)'
            diagnostic_get_header(7) = 'model_sig(4)'

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_header', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a diagnostic and a position.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp method of the subclass instance @ref diagnostic_class was
!>  constructed with. Only diagnostics that use te, ne, or sxrem profile have a
!>  valid implementation.
!>
!>  @see sxrem::sxrem_get_gp_i
!>  @see intpol::intpol_get_gp_i
!>  @see thomson::thomson_get_gp_i
!>  @see ece::ece_get_gp_i
!>
!>  @param[in] this         A @ref diagnostic_class instance.
!>  @param[in] a_model      A @ref model_class instance.
!>  @param[in] i            Index of the position for the kernel.
!>  @param[in] scale_factor Factor to scale the modeled signal.
!>  @returns Kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_gp_i(this, a_model, i, scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                        :: diagnostic_get_gp_i
      TYPE (diagnostic_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in)      :: a_model
      INTEGER, INTENT(in)                 :: i
      REAL (rprec), INTENT(in)            :: scale_factor

!  local variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (diagnostic_sxrem_type)
            diagnostic_get_gp_i = sxrem_get_gp(this%sxrem, a_model, i,         &
     &                                         scale_factor)

         CASE (diagnostic_intpol_type)
            diagnostic_get_gp_i = intpol_get_gp(this%intpol, a_model, i,       &
     &                                          scale_factor)

         CASE (diagnostic_thomson_type)
            diagnostic_get_gp_i = thomson_get_gp(this%thomson, a_model,        &
     &                                           i, scale_factor)

         CASE (diagnostic_ece_type)
            diagnostic_get_gp_i = ece_get_gp(this%ece, a_model, i,             &
     &                                       scale_factor)

         CASE DEFAULT
            diagnostic_get_gp_i = 0.0

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_gp_i', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a diagnostic and a signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp method of the subclass instance @ref diagnostic_class was
!>  constructed with. Only diagnostics that use te, ne, or sxrem profile have a
!>  valid implementation. The second signal was cast to a void pointer so it can
!>  be send to the subclass instance.
!>
!>  @see sxrem::sxrem_get_gp_s
!>  @see intpol::intpol_get_gp_s
!>  @see thomson::thomson_get_gp_s
!>  @see ece::ece_get_gp_s
!>
!>  @param[in] this         A @ref diagnostic_class instance.
!>  @param[in] a_model      A @ref model_class instance.
!>  @param[in] context      A @ref signal_class instance.
!>  @param     callback     Function pointer to call back to signal module.
!>  @param[in] scale_factor Factor to scale the modeled signal.
!>  @returns Kernel value for the position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_gp_s(this, a_model, context, callback,           &
     &                             scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                        :: diagnostic_get_gp_s
      TYPE (diagnostic_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in)      :: a_model
      CHARACTER (len=1), INTENT(in)       :: context(:)
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
      REAL (rprec), INTENT(in)            :: scale_factor

!  local variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (diagnostic_sxrem_type)
            diagnostic_get_gp_s = sxrem_get_gp(this%sxrem, a_model,            &
     &                                         context, callback,              &
     &                                         scale_factor)

         CASE (diagnostic_intpol_type)
            diagnostic_get_gp_s = intpol_get_gp(this%intpol, a_model,          &
     &                                          context, callback,             &
     &                                          scale_factor)

         CASE (diagnostic_thomson_type)
            diagnostic_get_gp_s = thomson_get_gp(this%thomson, context,        &
     &                                           callback, scale_factor)

         CASE (diagnostic_ece_type)
            diagnostic_get_gp_s = ece_get_gp(this%ece, a_model,                &
     &                                       context, callback,                &
     &                                       scale_factor)

         CASE DEFAULT
            diagnostic_get_gp_s = 0.0

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_gp_s', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the guassian process kernel for a diagnostic and a cartesian
!>  position.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_gp method of the subclass instance @ref diagnostic_class was
!>  constructed with. Only diagnostics that use te, ne, or sxrem profile have a
!>  valid implementation.
!>
!>  @see sxrem::sxrem_get_gp_x
!>  @see intpol::intpol_get_gp_x
!>  @see thomson::thomson_get_gp_x
!>  @see ece::ece_get_gp_x
!>
!>  @param[in] this         A @ref diagnostic_class instance.
!>  @param[in] a_model      A @ref model_class instance.
!>  @param[in] x_cart       The cartesian position of to get the kernel at.
!>  @param[in] scale_factor Factor to scale the modeled signal.
!>  @returns Kernel value for position and the signal.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_gp_x(this, a_model, x_cart, scale_factor)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: diagnostic_get_gp_x
      TYPE (diagnostic_class), INTENT(in)    :: this
      TYPE (model_class), INTENT(in)         :: a_model
      REAL (rprec), DIMENSION(3), INTENT(in) :: x_cart
      REAL (rprec), INTENT(in)               :: scale_factor

!  local variables
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      SELECT CASE(this%type)

         CASE (diagnostic_sxrem_type)
            diagnostic_get_gp_x = sxrem_get_gp(this%sxrem, a_model,            &
     &                                         x_cart, scale_factor)

         CASE (diagnostic_intpol_type)
            diagnostic_get_gp_x = intpol_get_gp(this%intpol, a_model,          &
     &                                          x_cart, scale_factor)

         CASE (diagnostic_thomson_type)
            diagnostic_get_gp_x = thomson_get_gp(this%thomson, a_model,        &
     &                                           x_cart, scale_factor)

         CASE (diagnostic_ece_type)
            diagnostic_get_gp_x = ece_get_gp(this%ece, a_model,                &
     &                                       x_cart, scale_factor)

         CASE DEFAULT
            diagnostic_get_gp_x = 0.0

      END SELECT

      CALL profiler_set_stop_time('diagnostic_get_gp_x', start_time)

      END FUNCTION

      END MODULE
