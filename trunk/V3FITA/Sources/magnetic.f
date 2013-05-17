!*******************************************************************************
!  File magnetic.f
!  Contains module magnetic
!  Defines the base class of the type magnetic:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE magnetic

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************

      TYPE magnetic_class
      END TYPE magnetic_class

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION magnetic_construct()

      IMPLICIT NONE
      TYPE (magnetic_class) :: magnetic_construct

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE magnetic_destruct(this)

      IMPLICIT NONE

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2),       :: diagnostic_get_modeled_signal
      TYPE (magnetic_class), INTENT(in) :: this
      TYPE (model), INTENT(in)          :: a_model

!  Start of executable code
      magnetic_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
