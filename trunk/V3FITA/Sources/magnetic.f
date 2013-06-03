!*******************************************************************************
!  File magnetic.f
!  Contains module magnetic
!  Defines the base class of the type magnetic:
!
!  Super Class: diagnostic
!  Sub Classes: None
!*******************************************************************************

      MODULE magnetic

      USE stel_kinds, only: rprec

      IMPLICIT NONE

!*******************************************************************************
!  magnetic module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: magnetic_no_type = -1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) magnetic base class
!
!*******************************************************************************

      TYPE magnetic_class
         INTEGER      :: type
         REAL (rprec) :: flux_factor
         LOGICAL      :: l_mdcoil_def
      END TYPE magnetic_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a magnetic.
!-------------------------------------------------------------------------------
!      FUNCTION magnetic_construct(l_mdcoil_def, flux_factor)
      FUNCTION magnetic_construct(flux_factor)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER    :: magnetic_construct
      REAL(rprec), INTENT(in), OPTIONAL :: flux_factor

!  Start of executable code
      ALLOCATE(magnetic_construct)

!      magnetic_construct%l_mdcoil_def = l_mdcoil_def
      IF (PRESENT(flux_factor)) THEN
         magnetic_construct%flux_factor = flux_factor
      ELSE
         magnetic_construct%flux_factor = 1.0  
      END IF

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE magnetic_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (magnetic_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%l_mdcoil_def = .false.
      this%flux_factor = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION magnetic_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: magnetic_get_modeled_signal
      TYPE (magnetic_class), INTENT(in)       :: this
      TYPE (model_class), INTENT(in)          :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Start of executable code
      magnetic_get_modeled_signal = 0.0
      sigma = 0.0

      END FUNCTION

      END MODULE
