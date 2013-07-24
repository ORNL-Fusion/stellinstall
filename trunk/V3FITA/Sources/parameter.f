!*******************************************************************************
!>  @file parameter.f
!>  @brief Contains module @ref param.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref param_class
!*******************************************************************************

      MODULE param

      USE stel_kinds, only: rprec
      USE data_parameters
      USE model

      IMPLICIT NONE

!*******************************************************************************
!  param module parameters
!*******************************************************************************
!>  No bounding type specified.
      INTEGER, PARAMETER :: param_range_no_type = -1
!>  Parameter value is bounded by a value.
      INTEGER, PARAMETER :: param_range_value_type = 0
!>  Parameter value is unbounded.
      INTEGER, PARAMETER :: param_range_infinity_type = 1
!>  Parameter value is bounded by another parameter.
      INTEGER, PARAMETER :: param_range_parameter_type = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) parameter base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a reconstructed parameter. An upper and lower bound
!>  may be set for the parameter value.
!>  @todo Add signal effectiveness.
!-------------------------------------------------------------------------------
      TYPE param_class
!>  Id number of the parameter. Id's are provided by either the @ref model or
!>  the @ref equilibrium.
         INTEGER                     :: param_id = data_no_id
!>  The i and j indices of the parameter. These are only used if the parameter
!>  is the element of an array.
         INTEGER, DIMENSION(2)       :: indicies = 0
!>  The maximum increment size of the parameter for calculating the jacobian.
!>  @see param_increment
         REAL (rprec)                :: vrnc = 0.0

!  Parameter ranges
!>  Type descriptor of the boundry type for the lower(1) and upper(2) ranges.
!>  @par Possible values are:
!>  * @ref param_range_no_type
!>  * @ref param_range_value_type
!>  * @ref param_range_infinity_type
!>  * @ref param_range_parameter_type
         INTEGER, DIMENSION(2) :: range_type = param_range_no_type
!>  If the @ref param_class::range_type is @ref param_range_parameter_type, this
!>  holds the parameter id's of the lower(1) and upper(2) boundary. Id's are
!>  provided by either the @ref model or the @ref equilibrium.
         INTEGER, DIMENSION(2)               :: range_id = data_no_id
!>  If the @ref param_class::range_type is @ref param_range_parameter_type, this
!>  holds the parameter indices of the lower(1,:) and upper(2,:) boundary.
         INTEGER, DIMENSION(2,2)             :: range_indicies = 0
!>  If the @ref param_class::range_type is @ref param_range_value_type, this
!>  holds the value of the lower(1) and upper(2) boundary.
         REAL (rprec), DIMENSION(2)          :: range_value = 0.0

!>  The contraints can change the normalization size of the patarameter when
!>  incrementing it. Store the actual change in parameter here for
!>  normalization. This value is set by @ref param_increment.
         REAL (rprec)                        :: delta

!>  Stored value of the parameter uncertainty for the current reconstruction
!>  step. This value is calculated by
!>  @ref reconstruction::reconstruction_eval_pos_covar.
         REAL (rprec)                        :: sigma = 0.0
!>  Stored row of the signal effectiveness matrix.
         REAL (rprec), DIMENSION(:), POINTER :: sem => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Pointer to a parameter object. Used for creating arrays of signal pointers.
!>  This is needed because fortran does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE param_pointer
         TYPE (param_class), POINTER :: p => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref param_class object.
!>
!>  Allocates memory and initializes a @ref param_class object. Parameters are
!>  converted from strings to and internal id prepresentation by either the
!>  @ref model or the @ref equilibrium.
!>
!>  Parameter boundaries are specified by a range type. This range type may be
!>  one of three types. 'infinity' represents an unbounded case. 'value'
!>  represents a bounds by a fixed value. By specifying another parameter name,
!>  the boundary will be bounded by a fixed or varying parameter. The bounding
!>  parameter string is converted into an internal id representation.
!>  @param[in] a_model        An instance of a @ref model object.
!>  @param[in] param_name     The name the parameter.
!>  @param[in] indicies       The array indices of the parameter.
!>  @param[in] vrnc           The maximum incremental step size.
!>  @param[in] range_type     Range type of the bounds. Values can be
!>                            'infinity', 'value' or the name of a parameter.
!>  @param[in] range_indicies The array indices of the boundary parameter.
!>  @param[in] range_value    Value of a boundary.
!>  @param[in] num_signals    Number of signals. Used to determine the size of
!>                            the @ref param_class::sem arrays.
!>  @returns A pointer to a constructed @ref param_class object.
!-------------------------------------------------------------------------------
      FUNCTION param_construct(a_model, param_name, indicies, vrnc,            &
     &                         range_type, range_indicies, range_value,        &
     &                         num_signals)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), POINTER                 :: param_construct
      TYPE (model_class), INTENT(in)              :: a_model
      CHARACTER (len=*), INTENT(in)               :: param_name
      INTEGER, DIMENSION(2), INTENT(in)           :: indicies
      REAL (rprec), INTENT(in)                    :: vrnc
      CHARACTER (len=*), DIMENSION(2), INTENT(in) :: range_type
      INTEGER, DIMENSION(2,2), INTENT(in)         :: range_indicies
      REAL (rprec), DIMENSION(2), INTENT(in)      :: range_value
      INTEGER, INTENT(in)                         :: num_signals

!  local variables
      INTEGER                                     :: id

!  Start of executable code
      id = model_get_param_id(a_model, TRIM(param_name))
      IF (.not.model_is_recon_param(a_model, id)) THEN
         CALL err_fatal('param_construct: invalid ' //                         &
     &                  'reconstruction parameter')
      END IF

      ALLOCATE(param_construct)

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
     &         model_get_param_id(a_model, TRIM(range_type(1)))
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
     &         model_get_param_id(a_model, TRIM(range_type(2)))
            param_construct%range_indicies(2,:) = range_indicies(2,:)

      END SELECT

      ALLOCATE(param_construct%sem(num_signals))
      param_construct%sem = 0.0

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref param_class object.
!>
!>  Deallocates memory and uninitializes a @ref param_class object.
!>
!>  @param[inout] this A @ref param_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE param_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%param_id = data_no_id
      this%indicies = 0
      this%vrnc = 0.0
      this%range_type = param_range_no_type
      this%range_id = data_no_id
      this%range_indicies = 0
      this%range_value = 0.0

      DEALLOCATE(this%sem)
      this%sem => null()

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets the parameter value.
!>
!>  Sets the value of a reconstruction parameter. If the parameter is out side
!>  one of the bounding parameters, the value if clamped to the bounding value.
!>  The resulting bounded parameter value is set by ethier the @ref model or
!>  the @ref equilibrium.
!>
!>  @param[in]    this    A @ref param_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[in]    value   The value to set the parameter to.
!-------------------------------------------------------------------------------
      SUBROUTINE param_set_value(this, a_model, value)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), INTENT(in)    :: this
      TYPE (model_class), INTENT(inout) :: a_model
      REAL (rprec), INTENT(in)          :: value

!  Local variables
      REAL (rprec)                      :: set_value
      REAL (rprec)                      :: upper_value
      REAL (rprec)                      :: lower_value

!  Start of executable code
      set_value = value

!  Check if the value is within the lower range. If the value is outside,
!  truncate the value to the lower range.
      IF (this%range_type(1) .ne. param_range_infinity_type) THEN
         lower_value = param_get_lower_range_value(this, a_model)

         IF (set_value .lt. lower_value) THEN
            CALL model_set_param(a_model, this%param_id,                       &
     &              this%indicies(1), this%indicies(2), lower_value)
            set_value = lower_value
            RETURN
         END IF
      END IF

!  The value is in the lower range. Check upper range. If the value is outside,
!  truncate the value to the upper range.
      IF (this%range_type(2) .ne. param_range_infinity_type) THEN
         upper_value = param_get_lower_range_value(this, a_model)

         IF (set_value .gt. upper_value) THEN
            CALL model_set_param(a_model, this%param_id,                       &
     &              this%indicies(1), this%indicies(2), upper_value)
            set_value = upper_value
            RETURN
         END IF
      END IF

!  The value is in range. Set the value to this.
      CALL model_set_param(a_model, this%param_id,                             &
     &                     this%indicies(1), this%indicies(2),                 &
     &                     set_value)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the parameter value.
!>
!>  Gets the current parameter value from either the @ref model or the
!>  @ref equilibrium.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @returns The parameter value.
!-------------------------------------------------------------------------------
      FUNCTION param_get_value(this, a_model)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: param_get_value
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model

!  Start of executable code
      param_get_value = model_get_param_value(a_model,                         &
     &                                        this%param_id,                   &
     &                                        this%indicies(1),                &
     &                                        this%indicies(2))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the parameter name.
!>
!>  Converts the parameter id from either the @ref model or the
!>  @ref equilibrium.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @returns The parameter name.
!-------------------------------------------------------------------------------
      FUNCTION param_get_name(this, a_model)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: param_get_name
      TYPE (param_class), INTENT(in)   :: this
      TYPE (model_class), INTENT(in)   :: a_model

!  Start of executable code
      param_get_name = model_get_param_name(a_model, this%param_id)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the lower boundary value.
!>
!>  Gets lower boundary value. The results of this function are invalid if the
!>  bounding type is @ref param_range_infinity_type or @ref param_range_no_type.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @returns The lower bounding value.
!-------------------------------------------------------------------------------
      FUNCTION param_get_lower_range_value(this, a_model)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: param_get_lower_range_value
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model

!  Start of executable code
      SELECT CASE (this%range_type(1))

         CASE DEFAULT
            param_get_lower_range_value = this%range_value(1)

         CASE (param_range_parameter_type)
            param_get_lower_range_value =                                      &
     &         model_get_param_value(a_model, this%range_id(1),                &
     &                               this%range_indicies(1,1),                 &
     &                               this%range_indicies(1,2))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the upper boundary value.
!>
!>  Gets upper boundary value. The results of this function are invalid if the
!>  bounding type is @ref param_range_infinity_type or @ref param_range_no_type.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @returns The upper bounding value.
!-------------------------------------------------------------------------------
      FUNCTION param_get_upper_range_value(this, a_model)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                   :: param_get_upper_range_value
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model

!  Start of executable code
      SELECT CASE (this%range_type(2))

         CASE DEFAULT
            param_get_upper_range_value = this%range_value(2)

         CASE (param_range_parameter_type)
            param_get_upper_range_value =                                      &
     &         model_get_param_value(a_model, this%range_id(2),                &
     &                               this%range_indicies(2,1),                 &
     &                               this%range_indicies(2,2))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the lower boundary type description.
!>
!>  Gets lower boundary type description.
!>
!>  @param[in] this A @ref param_class instance.
!>  @returns A string description of the lower boundary type.
!-------------------------------------------------------------------------------
      FUNCTION param_get_lower_range_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: param_get_lower_range_type
      TYPE (param_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE (this%range_type(1))

         CASE DEFAULT
            param_get_lower_range_type = 'No Type'

         CASE (param_range_value_type)
            param_get_lower_range_type = 'value'

         CASE (param_range_infinity_type)
            param_get_lower_range_type = 'infinity'

         CASE (param_range_parameter_type)
            param_get_lower_range_type = 'parameter'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets the upper boundary type description.
!>
!>  Gets upper boundary type description.
!>
!>  @param[in] this A @ref param_class instance.
!>  @returns A string description of the upper boundary type.
!-------------------------------------------------------------------------------
      FUNCTION param_get_upper_range_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: param_get_upper_range_type
      TYPE (param_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE (this%range_type(2))

         CASE DEFAULT
            param_get_upper_range_type = 'No Type'

         CASE (param_range_value_type)
            param_get_upper_range_type = 'value'

         CASE (param_range_infinity_type)
            param_get_upper_range_type = 'infinity'

         CASE (param_range_parameter_type)
            param_get_upper_range_type = 'parameter'

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Checks if the value is in the lower range.
!>
!>  If the value is less than the bounding value, the parameter value is out of
!>  range.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @param[in] value   A parameter value to check.
!>  @returns True for in range and false for out of range.
!-------------------------------------------------------------------------------
      FUNCTION param_is_in_lower_range(this, a_model, value)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                        :: param_is_in_lower_range
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model
      REAL (rprec), INTENT(in)       :: value

!  Start of executable code
      param_is_in_lower_range = .true.
      IF (this%range_type(1) .ne. param_range_infinity_type) THEN
         param_is_in_lower_range =                                             &
     &      value .ge. param_get_lower_range_value(this, a_model)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Checks if the value is in the upper range.
!>
!>  If the value is greater than the bounding value, the parameter value is out
!>  of range.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!>  @param[in] value   A parameter value to check.
!>  @returns True for in range and false for out of range.
!-------------------------------------------------------------------------------
      FUNCTION param_is_in_upper_range(this, a_model, value)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                        :: param_is_in_upper_range
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model
      REAL (rprec), INTENT(in)       :: value

!  Start of executable code
      param_is_in_upper_range = .true.
      IF (this%range_type(2) .ne. param_range_infinity_type) THEN
         param_is_in_upper_range =                                             &
     &      value .le. param_get_upper_range_value(this, a_model)
      END IF

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Increments the parameter value.
!>
!>  The value is incremented by the vrnc. The incremented value is then checked
!>  if the new value is in range. If the value is out of range, the value is
!>  incremented in the opposite direction. If the incremented value is outside
!>  of both ranges, the step size value is halved and the procedure is repeated.
!>  Five attempts are made until the function quits with an error. If
!>  successful, the actual step size used is cached in @ref param_class::delta
!>  for normalization purposes.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] a_model A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE param_increment(this, a_model)
      USE v3_utilities, only: err_fatal

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout) :: a_model

!  local variables
      REAL (rprec)                      :: new_value
      REAL (rprec)                      :: value
      REAL (rprec)                      :: step
      INTEGER                           :: itry

!  Start of executable code
      value = param_get_value(this, a_model)
!  Take the absolute value here to control the direction of the step.
      step = ABS(this%vrnc)

!  There are 5 attempts to increment the parameter before giving up.
      DO itry = 1, 5

!  Test to make sure there is a step
         IF (step .eq. 0.0) THEN
            CALL err_fatal('param_increment_single: tried to take ' //         &
     &                     'a step size of 0')
         END IF

!  Try a step in the positive direction. If the parameter is still in range set
!  otherwise.
         new_value = value + step
         IF (param_is_in_upper_range(this, a_model, new_value)) THEN
            CALL model_set_param(a_model, this%param_id,                       &
     &                           this%indicies(1), this%indicies(2),           &
     &                           new_value)
!  Store the actual change in value.
            this%delta = step
            RETURN
         END IF

!  Try a step in the negative direction.
         new_value = value - step
         IF (param_is_in_lower_range(this, a_model, new_value)) THEN
            CALL model_set_param(a_model, this%param_id,                       &
     &                           this%indicies(1), this%indicies(2),           &
     &                           new_value)
!  Store the actual change in value.
            this%delta = -step
            RETURN
         END IF

!  Decrease the step and try again.
         step = step/10.0
      END DO

      CALL err_fatal('param_increment: failed to change ' //                   &
     &               'reconstruction. Try decreaing rp_vrnc or ' //            &
     &               'expanding the range.')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Writes out a parameter to an output file.
!>
!>  Parameter information is formated as a the index, name, index 1, index2,
!>  value and sigma.
!>
!>  @param[in] this    A @ref param_class instance.
!>  @param[in] iou     Input/output unit representing the file to write to.
!>  @param[in] index   The index of a parameter.
!>  @param[in] a_model A @ref model instance.
!>  @note Eventually parameter range informarion will be added as well.
!-------------------------------------------------------------------------------
      SUBROUTINE param_write(this, iou, index, a_model)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: iou
      INTEGER, INTENT(in)            :: index
      TYPE (model_class), INTENT(in) :: a_model

!  Start of executable code

      WRITE (iou,1000) index,                                                  &
     &                 param_get_name(this, a_model),                          &
     &                 this%indicies(1), this%indicies(2),                     &
     &                 param_get_value(this, a_model),                         &
     &                 this%sigma, this%vrnc,                                  &
     &                 param_get_upper_range_type(this),                       &
     &                 param_get_lower_range_type(this),                       &
     &                 param_get_upper_range_value(this, a_model),             &
     &                 param_get_lower_range_value(this, a_model),             &
     &                 this%range_indicies(2,:),                               &
     &                 this%range_indicies(1,:)
1000  FORMAT(i3,2x,a18,2(1x,i4),3(2x,es12.5),2(2x,a12),2(2x,es12.5),           &
     &       2(2x,i8,1x,i8))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Writes out parameter header information to an output file.
!>
!>  Parameter information is formated as a the index, name, index 1, index2,
!>  value and sigma. This is implmented as a static method and does not require
!>  a parameter instance. This should only be called once to produce a single
!>  header.
!>
!>  @param[in] iou A input/output representing the file to write to.
!>  @note Eventually parameter range informarion will be added as well.
!-------------------------------------------------------------------------------
      SUBROUTINE param_write_header(iou)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: iou

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,*) ' *** Reconstruction parameters'
      WRITE (iou,1000)
1000  FORMAT ('irp',2x,'p_type',13x,'inx1',1x,'inx2',2x,'value',9x,            &
     &        'sigma',9x,'vrnc',10x,'range_type_h',2x,'range_type_l',2x,       &
     &        'r_value_h',5x,'r_value_l',5x,                                   &
     &        'r_inx1_h',1x,'r_inx2_h',2x,'r_inx1_l',1x,'r_inx2_l')
      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out the parameter data for a step to the result netcdf file.
!>
!>  Writes out the parameter value and sigma to the result file.
!>
!>  @param[in] this           A @ref param_class instance.
!>  @param[in] a_model        The equilibrium model.
!>  @param[in] result_ncid    A NetCDF id of the result file.
!>  @param[in] current_step   The surrent reconstruction step.
!>  @param[in] index          A index of a parameter.
!>  @param[in] param_value_id NetCDF variable id of the parameter value.
!>                            @see v3fit_context
!>  @param[in] param_sigma_id NetCDF variable id of the parameter sigma.
!>                            @see v3fit_context
!>  @param[in] param_sem_id   NetCDF variable id of the parameter signal
!>                            effectiveness.
!>  @note This assumes that the parameter sigma has already been calculated.
!>  Name and index information is written when the result file is first created
!>  by @ref v3fit_context::v3fit_context_write_step_data.
!-------------------------------------------------------------------------------
      SUBROUTINE param_write_step_data(this, a_model, result_ncid,             &
     &                                 current_step, index,                    &
     &                                 param_value_id, param_sigma_id,         &
     &                                 param_sem_id)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (param_class), INTENT(in) :: this
      TYPE (model_class), INTENT(in) :: a_model
      INTEGER, INTENT(in)            :: result_ncid
      INTEGER, INTENT(in)            :: current_step
      INTEGER, INTENT(in)            :: index
      INTEGER, INTENT(in)            :: param_value_id
      INTEGER, INTENT(in)            :: param_sigma_id
      INTEGER, INTENT(in)            :: param_sem_id

!  local variables
      INTEGER                        :: status

!  Start of executable code
      status = nf_put_var1_double(result_ncid, param_value_id,                 &
     &                            (/ index, current_step /),                   &
     &                            param_get_value(this, a_model))
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_var1_double(result_ncid, param_sigma_id,                 &
     &                            (/ index, current_step /),                   &
     &                            this%sigma)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      status = nf_put_vara_double(result_ncid, param_sem_id,                   &
     &                            (/ 1, index, current_step /),                &
     &                            (/ SIZE(this%sem), 1, 1 /),                  &
     &                            this%sem)

      END SUBROUTINE

      END MODULE
