!*******************************************************************************
!  File prior.f
!  Contains module prior
!  Defines the base class of the type prior:
!
!  Super Class: None
!  Sub Classes: prior_range, prior_guassian
!*******************************************************************************

      MODULE prior

      USE stel_kinds, only : rprec
      USE prior_range
      USE prior_guassian

      IMPLICIT NONE

!*******************************************************************************
!  prior module parameters
!*******************************************************************************

      INTEGER, PARAMETER :: prior_label_length = 80

      INTEGER, PARAMETER :: prior_range_type = 0
      INTEGER, PARAMETER :: prior_guassian_type = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior base class
!
!*******************************************************************************

      TYPE prior_class
         TYPE(prior_range_class), POINTER    :: range => null()
         TYPE(prior_guassian_class), POINTER :: gaussian => null()
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
      TYPE (prior_class)                   :: prior_construct_range
      TYPE (prior_range_class), INTENT(in) :: range_object

!  Start of executable code
      prior_construct_range%range = range_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a prior containing a range.
!-------------------------------------------------------------------------------
      FUNCTION prior_construct_gaussian(gaussian_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class)                      :: prior_construct_range
      TYPE (prior_gaussian_class), INTENT(in) :: gaussian_object

!  Start of executable code
      prior_construct_range%range = range_object

      END FUNCTION

      END MODULE
