!*******************************************************************************
!  File magnetic.f
!  Contains module thomson
!  Defines the base class of the type thomson:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE thomson

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  diagnostic module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: thomson_no_type = -1
      INTEGER, PARAMETER :: thomson_temperature_type = 0
      INTEGER, PARAMETER :: thomson_density_type = 1
      INTEGER, PARAMETER :: thomson_pressure_type = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) thomson base class
!
!*******************************************************************************

      TYPE thomson_class
         INTEGER                    :: type
         REAL (rprec), DIMENSION(3) :: xcart
      END TYPE thomson_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a thomson.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct(type, xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct
      INTEGER, INTENT(in)                    :: type
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      ALLOCATE(thomson_construct)

      thomson_construct%type = type
      thomson_construct%xcart = xcart

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE thomson_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%type = thomson_no_type
      this%xcart = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_modeled_signal(this, a_model)
      USE model
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)        :: thomson_get_modeled_signal
      TYPE (thomson_class), INTENT(in)  :: this
      TYPE (model_class), INTENT(inout) :: a_model

!  Local variables

!  Start of executable code

      SELECT CASE(this%type)

         CASE (thomson_temperature_type)
            thomson_get_modeled_signal(1) =                                    &
     &         model_get_te(a_model, this%xcart)

         CASE (thomson_density_type)
            thomson_get_modeled_signal(1) =                                    &
     &         model_get_te(a_model, this%xcart)

         CASE (thomson_pressure_type)
            thomson_get_modeled_signal(1) =                                    &
     &         equilibrium_get_p(a_model%equilibrium, this%xcart)

      END SELECT

      END FUNCTION

      END MODULE
