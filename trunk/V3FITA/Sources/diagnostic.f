!*******************************************************************************
!  File diagnostic.f
!  Contains module diagnostic
!  Defines the base class of the type signal:
!
!  Super Class: signal
!  Sub Classes: magnetic, sxrem, intpol, thomson, extcurz
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
      INTEGER, PARAMETER :: diagnostic_no_type = -1
      INTEGER, PARAMETER :: diagnostic_magnetic_type = 0
      INTEGER, PARAMETER :: diagnostic_sxrem_type = 1
      INTEGER, PARAMETER :: diagnostic_intpol_type = 2
      INTEGER, PARAMETER :: diagnostic_thomson_type = 3
      INTEGER, PARAMETER :: diagnostic_extcurz_type = 4

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) diagnostic base class
!
!*******************************************************************************

      TYPE diagnostic_class
         INTEGER                       :: type = diagnostic_no_type
         TYPE(magnetic_class), POINTER :: magnetic => null()
         TYPE(sxrem_class), POINTER    :: sxrem => null()
         TYPE(intpol_class), POINTER   :: intpol => null()
         TYPE(thomson_class), POINTER  :: thomson => null()
         TYPE(extcurz_class), POINTER  :: extcurz => null()
      END TYPE diagnostic_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
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
!  Construct a diagnostic containing a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_magnetic(magnetic_object)

      IMPLICIT NONE
      TYPE (diagnostic_class), POINTER  :: diagnostic_construct_magnetic
      TYPE (magnetic_class), POINTER, INTENT(in) :: magnetic_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_magnetic)

      diagnostic_construct_magnetic%type = diagnostic_magnetic_type
      diagnostic_construct_magnetic%magnetic => magnetic_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a diagnostic containing a sxrem.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_sxrem(sxrem_object)

      IMPLICIT NONE
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_sxrem
      TYPE (sxrem_class), POINTER, INTENT(in)   :: sxrem_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_sxrem)

      diagnostic_construct_sxrem%type = diagnostic_sxrem_type
      diagnostic_construct_sxrem%sxrem => sxrem_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a diagnostic containing a intpol.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_intpol(intpol_object)

      IMPLICIT NONE
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_intpol
      TYPE (intpol_class), POINTER, INTENT(in) :: intpol_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_intpol)

      diagnostic_construct_intpol%type = diagnostic_intpol_type
      diagnostic_construct_intpol%intpol => intpol_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a diagnostic containing a thomson.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_thomson(thomson_object)

      IMPLICIT NONE
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_thomson
      TYPE (thomson_class), POINTER, INTENT(in) :: thomson_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_thomson)

      diagnostic_construct_thomson%type = diagnostic_thomson_type
      diagnostic_construct_thomson%thomson => thomson_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a diagnostic containing a extcurz.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_construct_extcurz(extcurz_object)

      IMPLICIT NONE
      TYPE (diagnostic_class), POINTER :: diagnostic_construct_extcurz
      TYPE (extcurz_class), POINTER, INTENT(in) :: extcurz_object

!  Start of executable code
      ALLOCATE(diagnostic_construct_extcurz)

      diagnostic_construct_extcurz%type = diagnostic_extcurz_type
      diagnostic_construct_extcurz%extcurz => extcurz_object

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE diagnostic_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_class), POINTER, INTENT(inout) :: this

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
!  Get the modeled signal
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
!  Get the diagnostic type
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
!  Get the diagnostic header
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(3) ::                        &
     &   diagnostic_get_header
      TYPE (diagnostic_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (diagnostic_magnetic_type)
            diagnostic_get_header = magnetic_get_header(this%magnetic)

         CASE DEFAULT
            diagnostic_get_header = 'N/A'

      END SELECT

      END FUNCTION

      END MODULE
