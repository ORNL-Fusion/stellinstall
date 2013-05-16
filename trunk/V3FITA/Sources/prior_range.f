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

      CONTAINS
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
      SELECT CASE (TRIM(ADJUSTL(this%lower_range_label)))
         CASE ('value')
            prior_range_is_in_range = value .ge. this%lower_value

         CASE ('parameter')
            prior_range_is_in_range = .true.

         CASE ('infinity', '-infinity')
            prior_range_is_in_range = .true.

         CASE DEFAULT
            CALL err_fatal('prior_range_is_in_range unrecongnized ' //         &
     &                     'lower_range_label: ' //                            &
     &                     this%lower_range_label)
       END SELECT

      SELECT CASE (TRIM(ADJUSTL(this%upper_range_label)))
         CASE ('value')
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                value .le. this%upper_value

         CASE ('parameter')
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                .true.

         CASE ('infinity', '-infinity')
            prior_range_is_in_range = prior_range_is_in_range .and.            &
     &                                .true.

         CASE DEFAULT
            CALL err_fatal('prior_range_is_in_range unrecongnized ' //         &
     &                     'upper_range_label: ' //                            &
     &                     this%upper_range_label)
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
