!*******************************************************************************
!  File prior.f
!  Contains module prior
!  Defines the base class of the type prior:
!
!  Super Class: None
!  Sub Classes: prior_range, prior_guassian
!*******************************************************************************

      MODULE prior

      USE prior_range
      USE prior_gaussian

      IMPLICIT NONE

!*******************************************************************************
!  prior module parameters
!*******************************************************************************

      INTEGER, PARAMETER :: prior_no_type = -1
      INTEGER, PARAMETER :: prior_range_type = 0
      INTEGER, PARAMETER :: prior_gaussian_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior base class
!
!*******************************************************************************

      TYPE prior_class
         INTEGER                             :: type = prior_no_type
         TYPE(prior_range_class), POINTER    :: range => null()
         TYPE(prior_gaussian_class), POINTER :: gaussian => null()
      END TYPE prior_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE prior_construct
         MODULE PROCEDURE prior_construct_range,                               &
     &                    prior_construct_gaussian
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a prior containing a range.
!-------------------------------------------------------------------------------
      FUNCTION prior_construct_range(range_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), POINTER          :: prior_construct_range
      TYPE (prior_range_class), POINTER, INTENT(in) :: range_object

!  Start of executable code
      ALLOCATE(prior_construct_range)

      prior_construct_range%type = prior_range_type
      prior_construct_range%range => range_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a prior containing a range.
!-------------------------------------------------------------------------------
      FUNCTION prior_construct_gaussian(gaussian_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), POINTER                      ::                      &
     &   prior_construct_gaussian
      TYPE (prior_gaussian_class), POINTER, INTENT(in) ::                      &
     &   gaussian_object

!  Start of executable code
      ALLOCATE(prior_construct_gaussian)

      prior_construct_gaussian%type = prior_gaussian_type
      prior_construct_gaussian%gaussian => gaussian_object

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE prior_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%range)) THEN
         CALL prior_range_destruct(this%range)
         this%range => null()
      END IF

      IF (ASSOCIATED(this%gaussian)) THEN
         CALL prior_gaussian_destruct(this%gaussian)
         this%gaussian => null()
      END IF

      this%type = prior_no_type

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
      FUNCTION prior_is_in_range(this, value)
      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: prior_is_in_range
      TYPE(prior_class), INTENT(in) :: this
      REAL(rprec), INTENT(in)       :: value

!  Start of executable code
      SELECT CASE(this%type)
         CASE (prior_range_type)
            prior_is_in_range = prior_range_is_in_range(this%range,            &
     &                                                  value)

         CASE DEFAULT
            prior_is_in_range = .true.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  Getters SUBROUTINES
!*******************************************************************************
      FUNCTION prior_get_g2(this, value)
      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec)                   :: prior_get_g2
      TYPE(prior_class), INTENT(in) :: this
      REAL(rprec), INTENT(in)       :: value

!  Start of executable code
      SELECT CASE(this%type)
         CASE (prior_range_type)
            prior_get_g2 = prior_range_get_g2(this%range, value)

         CASE (prior_gaussian_type)
            prior_get_g2 = prior_gaussian_get_g2(this%gaussian, value)

         CASE DEFAULT
            prior_get_g2 = 0.0

      END SELECT

      END FUNCTION

      END MODULE
