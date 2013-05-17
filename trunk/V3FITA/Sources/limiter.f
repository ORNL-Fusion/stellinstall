!*******************************************************************************
!  File limiter.f
!  Contains module limiter
!  Defines the base class of the type limiter:
!
!  Super Class: geometric
!  Sub Classes: None
!*******************************************************************************

      MODULE limiter

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) limiter base class
!
!*******************************************************************************

      TYPE limiter_class
      END TYPE limiter_class

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION limiter_construct()

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class) :: limiter_construct

!  Start of executable code
      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE limiter_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), INTENT(inout) :: this

!  Start of executable code
      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2) :: limiter_get_modeled_signal
      TYPE (limiter_class), INTENT(in) :: this
      TYPE (model), INTENT(in)         :: a_model

!  Start of executable code
      limiter_get_modeled_signal = 0.0

      END FUNCTION

      END MODULE
