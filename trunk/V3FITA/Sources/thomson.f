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
!  DERIVED-TYPE DECLARATIONS
!  1) thomson base class
!
!*******************************************************************************

      TYPE thomson_class
      END TYPE thomson_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class) :: thomson_construct

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE thomson_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), INTENT(in) :: this

!  Start of executable code

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2)       :: thomson_get_modeled_signal
      TYPE (thomson_class), INTENT(in) :: this
      TYPE (model), INTENT(in)         :: a_model

!  Start of executable code
      thomson_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
