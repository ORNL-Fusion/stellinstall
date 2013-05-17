!*******************************************************************************
!  File magnetic.f
!  Contains module sxrem
!  Defines the base class of the type sxrem:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE sxrem

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) sxrem base class
!
!*******************************************************************************

      TYPE sxrem_class
      END TYPE sxrem_class

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION sxrem_construct()

      IMPLICIT NONE
      TYPE (sxrem_class) :: sxrem_construct

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE sxrem_destruct(this)

      IMPLICIT NONE

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION sxrem_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2),    :: sxrem_get_modeled_signal
      TYPE (sxrem_class), INTENT(in) :: this
      TYPE (model), INTENT(in)            :: a_model

!  Start of executable code
      sxrem_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
