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
!  Construct a thomson measuring temperature.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_te(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_te
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_te =>                                                  &
     &   thomson_construct(thomson_temperature_type, xcart)

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a thomson.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_ne(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_ne
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_ne => thomson_construct(thomson_density_type,          &
     &                                          xcart)

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a thomson.
!-------------------------------------------------------------------------------
      FUNCTION thomson_construct_p(xcart)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (thomson_class), POINTER          :: thomson_construct_p
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart

!  Start of executable code
      thomson_construct_p => thomson_construct(thomson_pressure_type,          &
     &                                         xcart)

      END FUNCTION

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
      FUNCTION thomson_get_modeled_signal(this, a_model, sigma)
      USE model
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: thomson_get_modeled_signal
      TYPE (thomson_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Local variables

!  Start of executable code
      thomson_get_modeled_signal = 0.0

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

      sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the thomson type
!-------------------------------------------------------------------------------
      FUNCTION thomson_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: thomson_get_signal_type
      TYPE (thomson_class), INTENT(in) :: this

!  Start of executable code
      thomson_get_signal_type = 'thscte'

      END FUNCTION

      END MODULE
