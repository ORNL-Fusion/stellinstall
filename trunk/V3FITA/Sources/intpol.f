!*******************************************************************************
!  File magnetic.f
!  Contains module intpol
!  Defines the base class of the type intpol:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE intpol

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) intpol base class
!
!*******************************************************************************

      TYPE intpol_class
      END TYPE intpol_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION intpol_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class) :: intpol_construct

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE intpol_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (intpol_class), INTENT(inout) :: this

!  Start of executable code

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION intpol_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2)      :: intpol_get_modeled_signal
      TYPE (intpol_class), INTENT(in) :: this
      TYPE (model), INTENT(in)        :: a_model

!  Start of executable code
      intpol_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
