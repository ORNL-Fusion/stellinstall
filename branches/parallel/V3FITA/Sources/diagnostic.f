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

      IMPLICIT NONE

!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
!>  Type descriptor for the @ref diagnostic_class no sub class.
      INTEGER, PARAMETER :: diagnostic_no_type = -1
!>  Type descriptor for the @ref diagnostic_class subclass @ref magnetic
      INTEGER, PARAMETER :: diagnostic_magnetic_type = 0
!>  Type descriptor for the @ref diagnostic_class subclass @ref sxrem
      INTEGER, PARAMETER :: diagnostic_sxrem_type = 1
!>  Type descriptor for the @ref diagnostic_class subclass @ref intpol
      INTEGER, PARAMETER :: diagnostic_intpol_type = 2
!>  Type descriptor for the @ref diagnostic_class subclass @ref thomson
      INTEGER, PARAMETER :: diagnostic_thomson_type = 3
!>  Type descriptor for the @ref diagnostic_class subclass @ref extcurz
      INTEGER, PARAMETER :: diagnostic_extcurz_type = 4

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
      END TYPE diagnostic_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref diagnostic_class types using
!>  @ref diagnostic_construct_magnetic, @ref diagnostic_construct_sxrem,
!>  @ref diagnostic_construct_intpol, @ref diagnostic_construct_thomson or
!>  @ref diagnostic_construct_extcurz
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_construct
         MODULE PROCEDURE diagnostic_construct_magnetic,                       &
     &                    diagnostic_construct_sxrem,                          &
     &                    diagnostic_construct_intpol,                         &
     &                    diagnostic_construct_thomson,                        &
     &                    diagnostic_construct_extcurz
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
      TYPE (diagnostic_class), POINTER  :: diagnostic_construct_magnetic
      TYPE (magnetic_class), POINTER    :: magnetic_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_magnetic)

      diagnostic_construct_magnetic%type = diagnostic_magnetic_type
      diagnostic_construct_magnetic%magnetic => magnetic_object

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

!  Start of executable code
      ALLOCATE(diagnostic_construct_sxrem)

      diagnostic_construct_sxrem%type = diagnostic_sxrem_type
      diagnostic_construct_sxrem%sxrem => sxrem_object

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

!  Start of executable code
      ALLOCATE(diagnostic_construct_intpol)

      diagnostic_construct_intpol%type = diagnostic_intpol_type
      diagnostic_construct_intpol%intpol => intpol_object

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

!  Start of executable code
      ALLOCATE(diagnostic_construct_thomson)

      diagnostic_construct_thomson%type = diagnostic_thomson_type
      diagnostic_construct_thomson%thomson => thomson_object

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

!  Start of executable code
      ALLOCATE(diagnostic_construct_extcurz)

      diagnostic_construct_extcurz%type = diagnostic_extcurz_type
      diagnostic_construct_extcurz%extcurz => extcurz_object

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
!>
!>  @param[in]    this    A @ref diagnostic_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: diagnostic_get_modeled_signal
      TYPE (diagnostic_class), INTENT(inout)  :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
      SELECT CASE(this%type)
         CASE (diagnostic_magnetic_type)
            diagnostic_get_modeled_signal =                                    &
     &         magnetic_get_modeled_signal(this%magnetic, a_model,             &
     &                                     sigma)

         CASE (diagnostic_sxrem_type)
            diagnostic_get_modeled_signal =                                    &
     &         sxrem_get_modeled_signal(this%sxrem, a_model, sigma)

         CASE (diagnostic_intpol_type)
            diagnostic_get_modeled_signal =                                    &
     &         intpol_get_modeled_signal(this%intpol, a_model, sigma)

         CASE (diagnostic_thomson_type)
            diagnostic_get_modeled_signal =                                    &
     &         thomson_get_modeled_signal(this%thomson, a_model, sigma)

         CASE (diagnostic_extcurz_type)
            diagnostic_get_modeled_signal =                                    &
     &         extcurz_get_modeled_signal(this%extcurz, a_model, sigma)

      END SELECT

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

!  Start of executable code
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

      END SELECT

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
      TYPE (diagnostic_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (diagnostic_magnetic_type)
            diagnostic_get_header = magnetic_get_header(this%magnetic)

         CASE DEFAULT
            diagnostic_get_header(1:3) = 'N/A'
            diagnostic_get_header(4) = 'model_sig(1)'
            diagnostic_get_header(5) = 'model_sig(2)'
            diagnostic_get_header(6) = 'model_sig(3)'
            diagnostic_get_header(7) = 'model_sig(4)'

      END SELECT

      END FUNCTION

      END MODULE
