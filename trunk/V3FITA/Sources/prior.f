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
      TYPE (prior_class)                   :: prior_construct_range
      TYPE (prior_range_class), INTENT(in) :: range_object

!  Start of executable code
      prior_construct_range%type = prior_range_type
      prior_construct_range%range = range_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a prior containing a range.
!-------------------------------------------------------------------------------
      FUNCTION prior_construct_gaussian(gaussian_object)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class)                      ::                               &
     &   prior_construct_gaussian
      TYPE (prior_gaussian_class), INTENT(in) :: gaussian_object

!  Start of executable code
      prior_construct_gaussian%type = prior_gaussian_type
      prior_construct_gaussian%gaussian = gaussian_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a prior_range.
!-------------------------------------------------------------------------------
      FUNCTION prior_range_construct(lower_range_label,                        &
     &                               upper_range_label,                        &
     &                               lower_value,                              &
     &                               upper_value,                              &
     &                               lower_param_indices,                      &
     &                               upper_param_indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class) :: prior_range_construct
      CHARACTER (len=prior_label_length), INTENT(in) ::                        &
     &   lower_range_label
      CHARACTER (len=prior_label_length), INTENT(in) ::                        &
     &   upper_range_label
      REAL(rprec), INTENT(in)           :: lower_value
      REAL(rprec), INTENT(in)           :: upper_value
      INTEGER, DIMENSION(2), INTENT(in) :: lower_param_indices
      INTEGER, DIMENSION(2), INTENT(in) :: upper_param_indices

!  Local Variables
      TYPE (prior_range_class), ALLOCATABLE :: range_object
      INTEGER                               :: status

!  Start of executable code
      ALLOCATE(range_object, STAT=status)

      SELECT CASE (TRIM(ADJUSTL(lower_range_label)))
         CASE ('value')
            range_object%lower_value = lower_value

         CASE ('parameter')
            range_object%lower_param_indices = lower_param_indices

         CASE ('infinity', '-infinity')
!  Do nothing

         CASE DEFAULT
            CALL err_fatal('prior_range_construct unrecongnized ' //           &
     &                     'lower_range_label: ' // lower_range_label)
      END SELECT

      range_object%lower_range_label = TRIM(ADJUSTL(lower_range_label))

      SELECT CASE (TRIM(ADJUSTL(upper_range_label)))
         CASE ('value')
            range_object%upper_value = upper_value

         CASE ('parameter')
            range_object%upper_param_indices = upper_param_indices

         CASE ('infinity', '-infinity')
!  Do nothing

         CASE DEFAULT
            CALL err_fatal('prior_range_construct unrecongnized ' //           &
     &                     'upper_range_label: ' // upper_range_label)
      END SELECT

      range_object%upper_range_label = TRIM(ADJUSTL(upper_range_label))

      prior_range_construct = prior_construct(range_object)

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a prior_gaussian.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_construct(mu, sigma)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class)      :: prior_gaussian_construct
      REAL(rprec), INTENT(in) :: mu
      REAL(rprec), INTENT(in) :: sigma

!  Local Variables
      TYPE (prior_gaussian_class), ALLOCATABLE :: gaussian_object
      INTEGER                                  :: status

!  Start of executable code
      ALLOCATE(gaussian_object, STAT=status)

      gaussian_object%mu = mu
      gaussian_object%sigma = sigma

      prior_gaussian_construct = prior_construct(gaussian_object)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE prior_destruct(this)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%range)) THEN
         DEALLOCATE(this%range)
         this%range => null()
      END IF

      IF (ASSOCIATED(this%gaussian)) THEN
         DEALLOCATE(this%gaussian)
         this%gaussian => null()
      END IF

      this%type = prior_no_type

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
