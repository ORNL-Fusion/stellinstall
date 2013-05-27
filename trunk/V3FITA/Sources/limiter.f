!*******************************************************************************
!  File limiter.f
!  Contains module limiter
!  Defines the base class of the type limiter:
!
!  Super Class: geometric
!  Sub Classes: None
!*******************************************************************************

      MODULE limiter

      USE stel_kinds, only: rprec, dp
      USE limiter_iso_T

      IMPLICIT NONE
!*******************************************************************************
!  geometric module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: limiter_no_type = -1
      INTEGER, PARAMETER :: limiter_iso_type = 0
      INTEGER, PARAMETER :: limiter_polygon_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) limiter base class
!
!*******************************************************************************

      TYPE limiter_class
         INTEGER                     :: type
         LOGICAL                     :: on_edge
         TYPE (limiter_iso), POINTER :: iso => null()
      END TYPE limiter_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE limiter_construct
         MODULE PROCEDURE limiter_construct_iso,                               &
     &                    limiter_construct_polygon
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a limiter containing an iso function.
!-------------------------------------------------------------------------------
      FUNCTION limiter_construct_iso(iso, on_edge)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER           :: limiter_construct_iso
      TYPE (limiter_iso), POINTER, INTENT(in) :: iso
      LOGICAL, INTENT(in)                     :: on_edge

!  Start of executable code
      ALLOCATE(limiter_construct_iso)

      limiter_construct_iso%type = limiter_iso_type
      limiter_construct_iso%on_edge = on_edge
      limiter_construct_iso%iso => iso

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a limiter containing a polygon.
!-------------------------------------------------------------------------------
      FUNCTION limiter_construct_polygon(on_edge)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER :: limiter_construct_polygon
      LOGICAL, POINTER, INTENT(in)  :: on_edge

!  Start of executable code
      ALLOCATE(limiter_construct_polygon)

      limiter_construct_polygon%type = limiter_iso_type
      limiter_construct_polygon%on_edge = on_edge

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE limiter_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%iso)) THEN
         CALL limiter_iso_destroy(this%iso)
         this%iso => null()
      END IF

      this%type = limiter_no_type
      this%on_edge = .false.

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_signal(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_signal
      TYPE (limiter_class), INTENT(in)  :: this
      TYPE (model_class), INTENT(inout) :: a_model

!  Start of executable code
      SELECT CASE (this%type)

         CASE (limiter_iso_type)
            limiter_get_modeled_signal =                                       &
     &         limiter_get_modeled_iso_signal(this, a_model)

         CASE (limiter_polygon_type)
            limiter_get_modeled_signal(1) = 0.0

      END SELECT

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the iso function signal
!-------------------------------------------------------------------------------
      FUNCTION limiter_get_modeled_iso_signal(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_iso_signal
      TYPE (limiter_class), INTENT(in) :: this
      TYPE (model_class), INTENT(inout)      :: a_model

!  local variables
      INTEGER                                :: nu, j, k
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: r_at_seq1, z_at_seq1
      REAL(rprec), DIMENSION(3)              :: rphiz_at_max, rphiz_edge
      REAL(rprec)                            :: fval_max, fval

!  Start of executable code

      END FUNCTION

      END MODULE
