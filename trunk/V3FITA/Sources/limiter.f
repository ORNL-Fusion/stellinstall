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
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_signal
      TYPE (limiter_class), INTENT(in) :: this
      TYPE (model), INTENT(inout)      :: a_model

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
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: limiter_get_modeled_iso_signal
      TYPE (limiter_class), INTENT(in) :: this
      TYPE (model), INTENT(inout)      :: a_model

!  local variables
      INTEGER                                :: nu, j, k
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: r_at_seq1, z_at_seq1
      REAL(rprec), DIMENSION(3)              :: rphiz_at_max, rphiz_edge
      REAL(rprec)                            :: fval_max, fval

!  Start of executable code
      fval_max = -1.0E10_dp

      nu = MAX(this%iso%numin,                                                 &
     &   a_model%eqstate%fixp%ns_array(a_model%eqstate%varp%ns_index),         &
     &   2*a_model%eqstate%fixp%mpol + 6)

      ALLOCATE(r_at_seq1(nu), z_at_seq1(nu))

!  loop around the s = 1 surface
!  Evaluate limiter function for each point
!  take maximum, store R, Phi, and Z at location of maximum
      rphiz_at_max = 0.0

      DO k = 1, SIZE(this%iso%vgrid, 1)
         CALL model_get_seq1_rz(a_model, this%iso%vgrid(k),                    &
     &                          r_at_seq1, z_at_seq1)

         DO j = 1, nu
            rphiz_edge(1) = r_at_seq1(j)
            rphiz_edge(2) = this%iso%vgrid(k)
            rphiz_edge(3) = z_at_seq1(j)

            CALL limiter_iso_f_cyl(this%iso, rphiz_edge, fval)

            IF (fval .gt. fval_max) THEN
               fval_max = fval
               rphiz_at_max = rphiz_edge
            END IF
         END DO
      END DO

      IF (this%on_edge) THEN
         limiter_get_modeled_iso_signal(1) = fval_max
      ELSE
         limiter_get_modeled_iso_signal(1) = MAX(fval_max, 0.0)
      ENDIF
      limiter_get_modeled_iso_signal(2:4) = rphiz_at_max

      DEALLOCATE(r_at_seq1, z_at_seq1)

      END FUNCTION

      END MODULE
