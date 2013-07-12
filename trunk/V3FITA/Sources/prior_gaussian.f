!*******************************************************************************
!>  @file prior_gaussian.f
!>  @brief Contains module prior_gaussian
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref prior_gaussian_class. This class
!>  implements priors of the type
!>
!>  P(B) = Exp(-((B - mu)/sigma)^2)
!>
!>  @par Super Class:
!>  @ref prior
!*******************************************************************************

      MODULE prior_gaussian

      USE data_parameters
      USE model
      USE stel_kinds, only : rprec, dp

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior_gaussian base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a prior_guassian signal.
!>  @par Super Class:
!>  @ref signal
!>  @par Sub Classes:
!>  @ref prior_guassian
!-------------------------------------------------------------------------------
      TYPE prior_gaussian_class
!>  Parameter id to check the guassian against.
         INTEGER               :: param_id = data_no_id
!>  Indices of the parameter.
         INTEGER, DIMENSION(2) :: indices = 0
      END TYPE prior_gaussian_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref prior_gaussian_class object.
!>
!>  Allocates memory and initializes a @ref prior_class object.
!>  @param[in] a_model    An instance of a @ref model object.
!>  @param[in] param_name The name the parameter.
!>  @param[in] indices    Indicies of the parameter.
!>  @returns A pointer to a constructed @ref prior_gaussian_class object.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_construct(a_model, param_name, indices)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER :: prior_gaussian_construct
      TYPE (model_class), INTENT(in)       :: a_model
      CHARACTER (len=*), INTENT(in)        :: param_name
      INTEGER, DIMENSION(2), INTENT(in)    :: indices

!  Start of executable code
      ALLOCATE(prior_gaussian_construct)

      prior_gaussian_construct%param_id =                                      &
     &   model_get_param_id(a_model, TRIM(param_name))
      prior_gaussian_construct%indices = indices

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref prior_gaussian_class object.
!>
!>  Deallocates memory and uninitializes a @ref prior_gaussian_class object.
!>
!>  @param[inout] this A @ref prior_gaussian_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE prior_gaussian_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%param_id = data_no_id
      this%indices = 0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  Getters SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the modeled signal.
!>
!>  The modeled signal is just the value of the parameter.
!>
!>  @param[in]    this    A @ref prior_gaussian_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @param[out]   sigma   The modeled sigma.
!>  @returns The model value.
!-------------------------------------------------------------------------------
      FUNCTION prior_gaussian_get_modeled_signal(this, a_model, sigma)

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec) :: prior_gaussian_get_modeled_signal
      TYPE(prior_gaussian_class), INTENT(in)  :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma

!  Start of executable code
      prior_gaussian_get_modeled_signal =                                      &
     &   model_get_param_value(a_model, this%param_id,                         &
     &                         this%indices(1), this%indices(2))
      sigma = 0.0

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
      FUNCTION prior_gaussian_get_signal_type(this)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: prior_gaussian_get_signal_type
      TYPE (prior_gaussian_class), INTENT(in) :: this

!  Start of executable code
      prior_gaussian_get_signal_type = 'guassian'

      END FUNCTION

      END MODULE
