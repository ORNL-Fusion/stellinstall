!*******************************************************************************
!  File prior_gaussian.f
!  Contains module prior_gaussian
!  Defines the base class of the type prior_gaussian:
!
!  Super Class: prior
!  Sub Classes: None
!
!  This class implements priors of the type that
!
!  P(B) = Exp(-((B - mu)/sigma)^2)
!
!*******************************************************************************

      MODULE prior_gaussian

      USE stel_kinds, only : rprec, dp

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) prior_gaussian base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  prior_gaussian_class
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
      TYPE prior_gaussian_class
         REAL(rprec) :: mu
         REAL(rprec) :: sigma
      END TYPE prior_gaussian_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION prior_gaussian_construct(mu, sigma)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER :: prior_gaussian_construct
      REAL(rprec), INTENT(in)              :: mu
      REAL(rprec), INTENT(in)              :: sigma

!  Local Variables
      INTEGER                              :: status

!  Start of executable code
      ALLOCATE(prior_gaussian_construct)

      prior_gaussian_construct%mu = mu
      prior_gaussian_construct%sigma = sigma

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE prior_gaussian_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (prior_gaussian_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%mu = 0.0
      this%sigma = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  Getters SUBROUTINES
!*******************************************************************************
      FUNCTION prior_gaussian_get_g2(this, value)
      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec)                            :: prior_gaussian_get_g2
      TYPE(prior_gaussian_class), INTENT(in) :: this
      REAL(rprec), INTENT(in)                :: value

!  Start of executable code
      prior_gaussian_get_g2 = ((value - this%mu)/this%sigma)**2.0

      END FUNCTION

      END MODULE
