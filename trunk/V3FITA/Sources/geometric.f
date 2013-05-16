!*******************************************************************************
!  File geometric.f
!  Contains module geometric
!  Defines the base class of the type signal:
!
!  Super Class: geometric
!  Sub Classes: limiter, boundary
!*******************************************************************************

      MODULE geometric

      USE stel_kinds, only: rprec
      USE limiter
      USE boundary

      IMPLICIT NONE

!*******************************************************************************
!  geometric module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: geometric_no_type = -1
      INTEGER, PARAMETER :: geometric_limiter_type = 0
      INTEGER, PARAMETER :: geometric_boundary_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) geometric base class
!
!*******************************************************************************

      TYPE geometric_class
         INTEGER                       :: type = geometric_no_type
         TYPE(limiter_class),  POINTER :: limiter => null()
         TYPE(boundary_class), POINTER :: boundary => null()
      END TYPE geometric_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE geometric_construct
      MODULE PROCEDURE geometric_construct_limiter,                            &
     &                 geometric_construct_boundary
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a geometric containing a limiter.
!-------------------------------------------------------------------------------
      FUNCTION geometric_construct_limiter(limiter_object)

      IMPLICIT NONE
      TYPE (geometric_class) :: geometric_construct_limiter
      TYPE (limiter_class), INTENT(in) :: limiter_class

!  Start of executable code
      geometric_construct_limiter%type = geometric_limiter_type
      geometric_construct_limiter%limiter = limiter_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a geometric containing a boundary.
!-------------------------------------------------------------------------------
      FUNCTION geometric_construct_boundary(boundary_object)

      IMPLICIT NONE
      TYPE (geometric_class) :: geometric_construct_boundary
      TYPE (boundary_class), INTENT(in) :: boundary_object

!  Start of executable code
      geometric_construct_boundary%type = geometric_boundary_type
      geometric_construct_boundary%boundary = boundary_object

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE geometric_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (geometric_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%limiter)) THEN
         limiter_destruct(this%limiter)
         this%limiter => null()
      END IF

      IF (ASSOCIATED(this%boundary)) THEN
         boundary_destruct(this%boundary)
         this%boundary => null()
      END IF

      this%type = geometric_no_type

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION geometric_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2),     :: geometric_get_modeled_signal
      TYPE (geometric_class), INTENT(in) :: this
      TYPE (model), INTENT(in)            :: a_model

!  Start of executable code
      SELECT CASE(this%type)
         CASE (geometric_limiter_type)
            geometric_get_modeled_signal =                                     &
     &         limiter_get_modeled_signal(this%limiter, a_model)

         CASE (diagnostic_boundary_type)
            geometric_get_modeled_signal =                                    &
     &         boundary_get_modeled_signal(this%boundary, a_model)

      END SELECT

      END FUNCTION

      END MODULE
