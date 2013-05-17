!*******************************************************************************
!  File magnetic.f
!  Contains module extcurz
!  Defines the base class of the type extcurz:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE extcurz

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) extcurz base class
!
!*******************************************************************************

      TYPE extcurz_class
      END TYPE extcurz_class

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_construct()

      IMPLICIT NONE
      TYPE (extcurz_class) :: extcurz_construct

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE extcurz_destruct(this)

      IMPLICIT NONE

!  Start of executable code

      END FUNCTION

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION extcurz_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2),    :: extcurz_get_modeled_signal
      TYPE (extcurz_class), INTENT(in) :: this
      TYPE (model), INTENT(in)            :: a_model

!  Start of executable code
      extcurz_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
