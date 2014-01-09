!*******************************************************************************
!>  @file prior.f
!>  @brief Contains module @ref prior
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref prior_class.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  Sub Classes: @ref prior_guassian
!*******************************************************************************

      MODULE prior

      USE prior_gaussian

      IMPLICIT NONE

!*******************************************************************************
!  prior module parameters
!*******************************************************************************
!>  Type descriptor for the @ref prior_class no sub class.
      INTEGER, PARAMETER :: prior_no_type = -1
!>  Type descriptor for the @ref diagnostic_class subclass @ref prior_guassian
      INTEGER, PARAMETER :: prior_gaussian_type = 0

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a prior signal.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref prior_guassian
!-------------------------------------------------------------------------------
      TYPE prior_class
!>  Type descirptor of the prior subclass type.
!>  @par Possible values are:
!>  * @ref prior_no_type
!>  * @ref prior_gaussian_type
         INTEGER                             :: type = prior_no_type
!>  An instance of a @ref sxrem object.
         TYPE(prior_gaussian_class), POINTER :: gaussian => null()
      END TYPE prior_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the construction of @ref prior_class types using
!>  @ref prior_construct_gaussian
!-------------------------------------------------------------------------------
      INTERFACE prior_construct
         MODULE PROCEDURE prior_construct_gaussian
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref prior_class containing a @ref prior_guassian object.
!>
!>  Allocates memory and initializes a @ref prior_class object.
!>
!>  @param[in] gaussian_object An instance of a @ref prior_guassian subclass.
!>  @returns A pointer to a constructed @ref prior_class object.
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
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref prior_class object.
!>
!>  Deallocates memory and uninitializes a @ref prior_class object.
!>
!>  @param[inout] this A @ref prior_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%gaussian)) THEN
         CALL prior_gaussian_destruct(this%gaussian)
         this%gaussian => null()
      END IF

      this%type = prior_no_type

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  Getters SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_modeled_signal method of the subclass instance
!>  @ref prior_class was constructed with.
!>
!>  @see prior_guassian::prior_guassian_get_modeled_signal
!>
!>  @param[in]    this    A @ref prior_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION prior_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec) :: prior_get_modeled_signal
      TYPE(prior_class), INTENT(in)           :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma

!  Start of executable code
      SELECT CASE(this%type)

         CASE (prior_gaussian_type)
            prior_get_modeled_signal =                                         &
     &         prior_gaussian_get_modeled_signal(this%gaussian,                &
     &                                           a_model, sigma)

         CASE DEFAULT
            prior_get_modeled_signal = 0.0
            sigma = 0.0

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the prior type.
!>
!>  This method is meant to be overwritten by a subclass method. As a result,
!>  returns the get_signal_type method of the subclass instance
!>  @ref prior_class was constructed with.
!>
!>  @see prior_gaussian::prior_gaussia_get_signal_type
!>
!>  @param[in] this A @ref prior_class instance.
!>  @returns A string describing the prior type.
!-------------------------------------------------------------------------------
      FUNCTION prior_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: prior_get_signal_type
      TYPE (prior_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (prior_gaussian_type)
            prior_get_signal_type = TRIM('prior ' //                           &
     &         prior_gaussian_get_signal_type(this%gaussian))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the model and model sigma array indices.
!>
!>  This method is meant to be overwritten by a subclass method if need be. As a
!>  result, returns the get_header method of the subclass instance
!>  @ref prior_class was constructed with. 
!>
!>  @param[in] this A @ref prior_class instance.
!>  @returns A string describing the model and model sigma array indices.
!-------------------------------------------------------------------------------
      FUNCTION prior_get_header(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(7) :: prior_get_header
      TYPE (prior_class), INTENT(in)                 :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE DEFAULT
            prior_get_header(1:3) = 'N/A'
            prior_get_header(4) = 'model_sig(1)'
            prior_get_header(5) = 'model_sig(2)'
            prior_get_header(6) = 'model_sig(3)'
            prior_get_header(7) = 'model_sig(4)'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Gets a discription of the prior type.
!>
!>  Returns a string describing the prior subclass.
!>
!>  @param[in] this A @ref prior_class instance.
!>  @returns A string describing the prior subclass.
!-------------------------------------------------------------------------------
      FUNCTION prior_get_prior_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: prior_get_prior_type
      TYPE (prior_class), INTENT(in)   :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (prior_gaussian_type)
            prior_get_prior_type = 'gaussian'
      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out auxiliary signal information to an output file.
!>
!>  Writes out the s_name and coefficient of the combined signals.
!>
!>  @param[in] this    A @ref prior_class instance.
!>  @param[in] iou     A input/output representing the file to write to.
!>  @param[in] index   A index of this signal.
!>  @param[in] a_model The equilibrium model.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_write_auxiliary(this, iou, index, a_model)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_class), INTENT(in) :: this
      INTEGER, INTENT(in)            :: iou
      INTEGER, INTENT(in)            :: index
      TYPE (model_class), INTENT(in) :: a_model

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,1000) index, prior_get_prior_type(this)

      SELECT CASE(this%type)

         CASE (prior_gaussian_type)
            CALL prior_gaussian_write_auxiliary(this%gaussian, iou,            &
     &                                          a_model)

      END SELECT

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a prior signal, type: ',a)

      END SUBROUTINE

      END MODULE
