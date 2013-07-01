!*******************************************************************************
!  File magnetic.f
!  Contains module extcurz
!  Defines the base class of the type extcurz:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE extcurz

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) extcurz base class
!
!*******************************************************************************

      TYPE extcurz_class
         REAL (rprec) :: r
         REAL (rprec) :: theta
      END TYPE extcurz_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
      FUNCTION extcurz_construct(r, theta)
      USE eq_T

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), POINTER  :: extcurz_construct
      REAL(rprec), INTENT(in)        :: r
      REAL(rprec), INTENT(in)        :: theta

!  Start of executable code
      ALLOCATE(extcurz_construct)
      extcurz_construct%r = r
      extcurz_construct%theta = theta

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE extcurz_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (extcurz_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%r = 0.0
      this%theta = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION extcurz_get_modeled_signal(this, a_model, sigma)
      USE model
      USE stel_constants, only: mu0

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: extcurz_get_modeled_signal
      TYPE (extcurz_class), INTENT(in)        :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
      extcurz_get_modeled_signal = 0.0

      extcurz_get_modeled_signal(1) =                                          &
     &   equilibrium_get_Int_B_dphi(a_model%equilibrium, this%r,               &
     &                              this%theta)/mu0

      sigma = 0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the extcurz type
!-------------------------------------------------------------------------------
      FUNCTION extcurz_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: extcurz_get_signal_type
      TYPE (extcurz_class), INTENT(in) :: this

!  Start of executable code
      extcurz_get_signal_type = 'extcur'

      END FUNCTION

      END MODULE
