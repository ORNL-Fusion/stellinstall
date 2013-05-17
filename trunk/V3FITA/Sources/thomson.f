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
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct(type, xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class)                   :: thomson_construct
      INTEGER, INTENT(in)                    :: type
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct%type = type
      thomson_construct%xcart = xcart

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE thomson_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), INTENT(inout) :: this

!  Start of executable code
      this%type = thomson_no_type
      this%xcart = 0.0

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_modeled_signal(this, a_model)
      USE model_T
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)       :: thomson_get_modeled_signal
      TYPE (thomson_class), INTENT(in) :: this
      TYPE (model), INTENT(inout)      :: a_model

!  Local variables
      REAL(rprec), DIMENSION(3) :: r_cyl

!  Start of executable code
      r_cyl = cart_to_cyl(this%xcart)

      SELECT CASE(this%type)

         CASE (thomson_temperature_type)
            thomson_get_modeled_signal(1) = model_get_te_xcyl(a_model,         &
     &                                                        r_cyl)

         CASE (thomson_density_type)
            thomson_get_modeled_signal(1) = model_get_ne_xcyl(a_model,         &
     &                                                        r_cyl)

         CASE (thomson_pressure_type)
            thomson_get_modeled_signal(1) = model_get_p_xcyl(a_model,          &
     &                                                       r_cyl)

      END SELECT

      END FUNCTION

      END MODULE
