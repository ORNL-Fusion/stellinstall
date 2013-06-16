!*******************************************************************************
!  File parameter.f
!  Contains module param
!  Defines the base class of the type parameter:
!
!  Super Class: None
!  Sub Classes: None
!*******************************************************************************

      MODULE param

      USE stel_kinds, only: rprec
      USE data_parameters
      USE model

      IMPLICIT NONE

!*******************************************************************************
!  param module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: param_range_no_type = -1
      INTEGER, PARAMETER :: param_range_value_type = 0
      INTEGER, PARAMETER :: param_range_infinity_type = 1
      INTEGER, PARAMETER :: param_range_parameter_type = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) parameter base class
!
!*******************************************************************************
      TYPE param_class
         TYPE (model_class), POINTER :: eq_model
         INTEGER                     :: param_id = data_no_id
         INTEGER, DIMENSION(2)       :: indicies = 0
         REAL (rprec)                :: vrnc = 0.0

!  Parameter ranges
         INTEGER, DIMENSION(2) :: range_type = param_range_no_type
         INTEGER, DIMENSION(2)       :: range_id = data_no_id
         INTEGER, DIMENSION(2,2)     :: range_indicies = 0
         REAL (rprec), DIMENSION(2)  :: range_value = 0.0
      END TYPE

!-------------------------------------------------------------------------------
!  Derived type to hold a pointer to a parameter. This is needed because fortran
!  does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE param_pointer
         TYPE (param_class), POINTER :: p => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE param_increment
         MODULE PROCEDURE param_increment_single,                              &
     &                    param_increment_size
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION param_construct(eq_model, param_name, indicies,                 &
     &                         vrnc, range_type, range_indicies,               &
     &                         range_value)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), POINTER                 :: param_construct
      TYPE (model_class), POINTER, INTENT(in)     :: eq_model
      CHARACTER (len=*), INTENT(in)               :: param_name
      INTEGER, DIMENSION(2), INTENT(in)           :: indicies
      REAL (rprec), INTENT(in)                    :: vrnc
      CHARACTER (len=*), DIMENSION(2), INTENT(in) :: range_type
      INTEGER, DIMENSION(2,2), INTENT(in)         :: range_indicies
      REAL (rprec), DIMENSION(2), INTENT(in)      :: range_value

!  local variables
      INTEGER                                     :: id

!  Start of executable code
      id = model_get_param_id(eq_model, TRIM(param_name))
      IF (.not.model_is_recon_param(eq_model,                                  &
     &                              param_construct%param_id)) THEN
         CALL err_fatal('param_construct: invalid ' //                         &
     &                  'reconstruction parameter')
      END IF

      ALLOCATE(param_construct)

      param_construct%eq_model => eq_model
      param_construct%param_id = id
      param_construct%indicies = indicies
      param_construct%vrnc = vrnc

      SELECT CASE (TRIM(range_type(1)))

         CASE ('value')
            param_construct%range_type(1) = param_range_value_type
            param_construct%range_value(1) = range_value(1)
        
         CASE ('infinity')
            param_construct%range_type(1) = param_range_infinity_type

         CASE DEFAULT
            param_construct%range_type(1) = param_range_parameter_type
            param_construct%range_id(1) =                                      &
     &         model_get_param_id(eq_model, TRIM(range_type(1)))
            param_construct%range_indicies(1,:) = range_indicies(1,:)

      END SELECT

      SELECT CASE (TRIM(range_type(2)))

         CASE ('value')
            param_construct%range_type(2) = param_range_value_type
            param_construct%range_value(2) = range_value(2)

         CASE ('infinity')
            param_construct%range_type(2) = param_range_infinity_type

         CASE DEFAULT
            param_construct%range_type(2) = param_range_parameter_type
            param_construct%range_id(2) =                                      &
     &         model_get_param_id(eq_model, TRIM(range_type(2)))
            param_construct%range_indicies(2,:) = range_indicies(2,:)

      END SELECT

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE param_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), POINTER, INTENT(inout) :: this

!  Start of executable code
!  Do not call the model destructor here. Parameter only point to a model but
!  don't own it.
      this%eq_model => null()
      this%param_id = data_no_id
      this%indicies = 0
      this%vrnc = 0.0
      this%range_type = param_range_no_type
      this%range_id = data_no_id
      this%range_indicies = 0
      this%range_value = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Set reconstruction parameter
!
!  This function returns the actual value, the range was set to.
!-------------------------------------------------------------------------------
      FUNCTION param_set_value(this, value)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: param_set_value
      TYPE (param_class), INTENT(inout) :: this
      REAL (rprec), INTENT(in)          :: value

!  Local variables
      REAL (rprec)                      :: upper_value
      REAL (rprec)                      :: lower_value

!  Start of executable code
      param_set_value = value

!  Check if the value is within the lower range. If the value is outside,
!  truncate the value to the lower range.
      IF (this%range_type(1) .ne. param_range_infinity_type) THEN
         lower_value = param_get_lower_range_value(this)

         IF (param_set_value .lt. lower_value) THEN
            CALL model_set_param(this%eq_model, this%param_id,                 &
     &              this%indicies(1), this%indicies(2), lower_value)
            param_set_value = lower_value
            RETURN
         END IF
      END IF

!  The value is in the lower range. Check upper range. If the value is outside,
!  truncate the value to the upper range.
      IF (this%range_type(2) .ne. param_range_infinity_type) THEN
         upper_value = param_get_lower_range_value(this)

         IF (param_set_value .gt. upper_value) THEN
            CALL model_set_param(this%eq_model, this%param_id,                 &
     &             this%indicies(1), this%indicies(2), upper_value)
            param_set_value = upper_value
            RETURN
         END IF
      END IF

!  The value is in range. Set the value to this.
      CALL model_set_param(this%eq_model, this%param_id,                       &
     &                     this%indicies(1), this%indicies(2),                 &
     &                     param_set_value)

      END FUNCTION

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the value of the parameter
!-------------------------------------------------------------------------------
      FUNCTION param_get_value(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: param_get_value
      TYPE (param_class), INTENT(inout) :: this

!  Start of executable code
      param_get_value = model_get_param_value(this%eq_model,                   &
     &                                        this%param_id,                   &
     &                                        this%indicies(1),                &
     &                                        this%indicies(2))

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the value of the lower range
!-------------------------------------------------------------------------------
      FUNCTION param_get_lower_range_value(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: param_get_lower_range_value
      TYPE (param_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%range_type(1))

         CASE DEFAULT
            param_get_lower_range_value = this%range_value(1)

         CASE (param_range_parameter_type)
            param_get_lower_range_value =                                      &
     &         model_get_param_value(this%eq_model, this%range_id(1),          &
     &                               this%range_indicies(1,1),                 &
     &                               this%range_indicies(1,2))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the value of the upper range
!-------------------------------------------------------------------------------
      FUNCTION param_get_upper_range_value(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: param_get_upper_range_value
      TYPE (param_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%range_type(2))

         CASE DEFAULT
            param_get_upper_range_value = this%range_value(2)

         CASE (param_range_parameter_type)
            param_get_upper_range_value =                                      &
     &         model_get_param_value(this%eq_model, this%range_id(2),          &
     &                               this%range_indicies(2,1),                 &
     &                               this%range_indicies(2,2))

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Query value is in the lower range.
!-------------------------------------------------------------------------------
      FUNCTION param_is_in_lower_range(this, value)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                        :: param_is_in_lower_range
      TYPE (param_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)       :: value

!  Start of executable code
      IF (this%range_type(1) .ne. param_range_infinity_type) THEN
         param_is_in_lower_range =                                             &
     &      value .ge. param_get_lower_range_value(this)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Query value is in the lower range.
!-------------------------------------------------------------------------------
      FUNCTION param_is_in_upper_range(this, value)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                        :: param_is_in_upper_range
      TYPE (param_class), INTENT(in) :: this
      REAL (rprec), INTENT(in)       :: value

!  Start of executable code
      IF (this%range_type(2) .ne. param_range_infinity_type) THEN
         param_is_in_upper_range =                                             &
     &      value .le. param_get_upper_range_value(this)
      END IF

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Increment the value by the variance.
!
!  Returns the actual step taken. Used for partial derivatives
!-------------------------------------------------------------------------------
      FUNCTION param_increment_single(this)
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: param_increment_single
      TYPE (param_class), INTENT(inout) :: this

!  local variables
      REAL (rprec)                      :: value
      REAL (rprec)                      :: step
      INTEGER                           :: itry

!  Start of executable code
      value = param_get_value(this)
!  Take the absolute value here to control the direction of the step.
      step = ABS(value*this%vrnc)

!  There are 5 attempts to increment the parameter before giving up.
      DO itry = 1, 5

!  Test to make sure there is a step
         IF (step .eq. 0.0) THEN
            CALL err_fatal('param_increment_single: tried to take ' //         &
     &                     'a step size of 0')
         END IF

!  Try a step in the positive direction. If the parameter is still in range set
!  otherwise.
         param_increment_single = value + step
         IF (param_is_in_upper_range(this,                                     &
     &                               param_increment_single)) THEN
            CALL model_set_param(this%eq_model, this%param_id,                 &
     &                           this%indicies(1), this%indicies(2),           &
     &                           param_increment_single)
            RETURN
         END IF

!  Try a step in the negative direction.
         param_increment_single = value - step
         IF (param_is_in_lower_range(this,                                     &
     &                              param_increment_single)) THEN
            CALL model_set_param(this%eq_model, this%param_id,                 &
     &                           this%indicies(1), this%indicies(2),           &
     &                           param_increment_single)
            RETURN
         END IF

!  Decrease the step and try again.
         step = step/10.0
      END DO

      CALL err_fatal('param_increment_single: failed to change ' //            &
     &               'reconstruction. Try decreaing rp_vrnc or ' //            &
     &               'expanding the range.')

      END FUNCTION

!-------------------------------------------------------------------------------
!  Increment the value by the variance and multiple number of steps.
!
!  Returns the actual step taken. USed for taking reconstruction steps.
!-------------------------------------------------------------------------------
      FUNCTION param_increment_size(this, step_size)
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                      :: param_increment_size
      TYPE (param_class), INTENT(inout) :: this
      REAL (rprec), INTENT(in)          :: step_size

!  local variables
      REAL (rprec)                      :: value

!  Start of executable code
      value = param_get_value(this)

      param_increment_size =                                                   &
     &   param_set_value(this, value + value*this%vrnc*step_size)

      END FUNCTION

      END MODULE
