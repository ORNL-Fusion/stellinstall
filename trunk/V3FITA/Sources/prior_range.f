!*******************************************************************************
!  File prior_range.f
!  Contains module prior_range
!  Defines the base class of the type prior_range:
!
!  Super Class: prior
!  Sub Classes: None
!
!  This class implements priors of the type that
!
!         { 0 if B < lower_range
!  P(B) = { 1 if B >= lower_range and B <= upper_range
!         { 0 if B > upper_range
!
!*******************************************************************************

      MODULE prior_range

      USE stel_kinds, only : rprec, dp
      USE v3_utilities

      IMPLICIT NONE

!*******************************************************************************
!  prior_range module parameters
!*******************************************************************************
      REAL(rprec), PARAMETER :: prior_range_large = 1.0E100_dp
      INTEGER, PARAMETER     :: prior_label_length = 80

      INTEGER, PARAMETER :: range_no_type = -1
      INTEGER, PARAMETER :: range_value_type = 0
      INTEGER, PARAMETER :: range_parameter_type = 1
      INTEGER, PARAMETER :: range_infinity_type = 2

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
!      'parameter_****'
!
!  When the range either range label is a value, then upper_value or lower_value
!  contain the bound values. When either range label is a parameter,
!  lower_param_indices or upper_param_indices contains the index of the
!  reconstruction parameter.
!-------------------------------------------------------------------------------
      TYPE prior_range_class
         INTEGER                            :: lower_range_type
         INTEGER                            :: upper_range_type
         CHARACTER (len=prior_label_length) :: lower_range_label
         CHARACTER (len=prior_label_length) :: upper_range_label
         REAL(rprec)                        :: lower_value
         REAL(rprec)                        :: upper_value
         INTEGER, DIMENSION(2)              :: lower_param_indices
         INTEGER, DIMENSION(2)              :: upper_param_indices
      END TYPE prior_range_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION prior_range_construct(lower_range_label,                        &
     &                               upper_range_label,                        &
     &                               lower_value,                              &
     &                               upper_value,                              &
     &                               lower_param_indices,                      &
     &                               upper_param_indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_range_class), POINTER :: prior_range_construct
      CHARACTER (len=prior_label_length), INTENT(in) ::                        &
     &   lower_range_label
      CHARACTER (len=prior_label_length), INTENT(in) ::                        &
     &   upper_range_label
      REAL(rprec), INTENT(in)           :: lower_value
      REAL(rprec), INTENT(in)           :: upper_value
      INTEGER, DIMENSION(2), INTENT(in) :: lower_param_indices
      INTEGER, DIMENSION(2), INTENT(in) :: upper_param_indices

!  Local Variables
      INTEGER                           :: status, i

!  Start of executable code
      ALLOCATE(prior_range_construct)

! Check to see where the first '_' character is located. This is used to split
! strings of the format 'parameter_param'
      i = INDEX(TRIM(ADJUSTL(lower_range_label)), '_')
      IF (i .eq. 0) THEN
         i = LEN_TRIM(lower_range_label) + 1
      ENDIF

      SELECT CASE (lower_range_label(1:i - 1))

         CASE ('value')
            prior_range_construct%lower_range_type = range_value_type
            prior_range_construct%lower_value = lower_value

         CASE ('parameter')
            prior_range_construct%lower_range_type                             &
     &         = range_parameter_type
            prior_range_construct%lower_param_indices                          &
     &         = lower_param_indices
            prior_range_construct%lower_range_label                            &
     &         = lower_range_label(i + 1:LEN_TRIM(lower_range_label))

         CASE ('infinity')
            prior_range_construct%lower_range_type = range_infinity_type

         CASE DEFAULT
            CALL err_fatal('prior_range_construct unrecongnized ' //           &
     &                     'lower_range_label: ' // lower_range_label)
      END SELECT

      prior_range_construct%lower_range_label                                  &
     &   = TRIM(ADJUSTL(lower_range_label))

! Check to see where the first '_' character is located. This is used to split
! strings of the format 'parameter_param'
      i = INDEX(TRIM(ADJUSTL(upper_range_label)), '_')
      IF (i .eq. 0) THEN
         i = LEN_TRIM(lower_range_label) + 1
      ENDIF

      SELECT CASE (upper_range_label(1:i - 1))

         CASE ('value')
            prior_range_construct%upper_range_type = range_value_type
            prior_range_construct%upper_value = upper_value

         CASE ('parameter')
            prior_range_construct%upper_range_type                             &
     &         = range_parameter_type
            prior_range_construct%upper_param_indices                          &
     &         = upper_param_indices
            prior_range_construct%upper_range_label                            &
     &         = upper_range_label(i + 1:LEN_TRIM(upper_range_label))

         CASE ('infinity')
            prior_range_construct%upper_range_type = range_infinity_type

         CASE DEFAULT
            CALL err_fatal('prior_range_construct unrecongnized ' //           &
     &                     'upper_range_label: ' // upper_range_label)
      END SELECT

      prior_range_construct%upper_range_label                                  &
     &   = TRIM(ADJUSTL(upper_range_label))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE prior_range_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_range_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%lower_range_type = range_no_type
      this%lower_range_label = ' '
      this%lower_value = 0.0
      this%lower_param_indices = 0

      this%upper_range_type = range_no_type
      this%upper_range_label = ' '
      this%upper_value = 0.0
      this%upper_param_indices = 0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
      FUNCTION prior_range_is_in_range(this, value)
      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                             :: prior_range_is_in_range
      TYPE(prior_range_class), INTENT(in) :: this
      REAL(rprec), INTENT(in)             :: value

!  Start of executable code
      SELECT CASE (this%lower_range_type)

         CASE (range_value_type)
            prior_range_is_in_range = value .ge. this%lower_value

         CASE (range_parameter_type)
            prior_range_is_in_range = .true. ! Set to true until parameters are implemented.

         CASE (range_infinity_type)
            prior_range_is_in_range = .true.

      END SELECT

      SELECT CASE (this%upper_range_type)
         CASE (range_value_type)
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                value .le. this%upper_value

         CASE (range_parameter_type)
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                .true. ! Set to true until parameters are implemented.

         CASE (range_infinity_type)
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                .true.

      END SELECT

      END FUNCTION

!*******************************************************************************
!  Getters SUBROUTINES
!*******************************************************************************
      FUNCTION prior_range_get_g2(this, value)
      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec)                         :: prior_range_get_g2
      TYPE(prior_range_class), INTENT(in) :: this
      REAL(rprec), INTENT(in)             :: value

!  Start of executable code
      IF (prior_range_is_in_range(this, value)) THEN
         prior_range_get_g2 = 0.0
      ELSE
         prior_range_get_g2 = prior_range_large
      END IF

      END FUNCTION

      END MODULE
