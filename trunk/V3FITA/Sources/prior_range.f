!*******************************************************************************
!  File prior.f
!  Contains module prior
!  Defines the base class of the type prior_range:
!
!  Super Class: prior
!  Sub Classes: None
!
!  This class implements prios of the type that
!
!         { 0 if B < lower_range
!  P(B) = { 1 if B >= lower_range and B <= upper_range
!         { 0 if B > upper_range
!
!*******************************************************************************

      MODULE prior_range

      USE stel_kinds, only : rprec
      USE prior

      IMPLICIT NONE

!*******************************************************************************
!  prior_range module parameters
!*******************************************************************************
      REAL(rprec), PARAMETER :: prior_range_large = 1.0E100

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior_range base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  prior_range_class
!
!  lower_range_label and upper_range_label hold the range type. Possible values
!  are:
!      'value'
!      'infinity'
!      '-infinity'
!      'parameter'
!
!  When the range either range label is a value, then upper_value or lower_value
!  contain the bound values. When either range label is a parameter,
!  lower_param_indices or upper_param_indices contains the index of the
!  reconstruction parameter.
!-------------------------------------------------------------------------------
      TYPE prior_range_class
         CHARACTER (len=prior_label_length) :: lower_range_label
         CHARACTER (len=prior_label_length) :: upper_range_label
         REAL(rprec)                        :: lower_value
         REAL(rprec)                        :: upper_value
         INTEGER, DIMENSION(2)              :: lower_param_indices
         INTEGER, DIMENSION(2)              :: upper_param_indices
      END TYPE prior_range_class

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a prior_range.
!-------------------------------------------------------------------------------
      FUNCTION prior_range_construct(lower_range_label,                        &
     &                               upper_range_label                         &
     &                               lower_value,                              &
     &                               upper_value,                              &
     &                               lower_param_indices,                      &
     &                               upper_param_indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class) :: prior_range_construct
      CHARACTER (len=prior_label_length), INTENT(in) :: lower_range_label
      CHARACTER (len=prior_label_length), INTENT(in) :: upper_range_label
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

      END MODULE
