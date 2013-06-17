!*******************************************************************************
!  File boundary.f
!  Contains module boundary
!  Defines the base class of the type boundary:
!
!  Super Class: geometric
!  Sub Classes: None
!*******************************************************************************

      MODULE boundary

      USE stel_kinds, only: rprec

      IMPLICIT NONE
!*******************************************************************************
!  boundary module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: boundary_no_type = -1
      INTEGER, PARAMETER :: boundary_rbc_type = 0
      INTEGER, PARAMETER :: boundary_zbs_type = 1
      INTEGER, PARAMETER :: boundary_rbs_type = 2
      INTEGER, PARAMETER :: boundary_zbc_type = 3

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) boundary base class
!
!*******************************************************************************

      TYPE boundary_class
         INTEGER     :: type = boundary_no_type
         INTEGER     :: m_index
         INTEGER     :: n_index
      END TYPE boundary_class

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION boundary_construct(coefficient_type, m_index, n_index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (boundary_class), POINTER :: boundary_construct
      CHARACTER (len=3), INTENT(in)  :: coefficient_type
      INTEGER, INTENT(in)            :: m_index
      INTEGER, INTENT(in)            :: n_index

!  Start of executable code
      ALLOCATE(boundary_construct)

      boundary_construct%m_index = m_index
      boundary_construct%n_index = n_index

      SELECT CASE (coefficient_type)
         CASE ('rbc')
            boundary_construct%type = boundary_rbc_type

         CASE ('zbs')
            boundary_construct%type = boundary_zbs_type

         CASE ('rbs')
            boundary_construct%type = boundary_rbs_type

         CASE ('zbc')
            boundary_construct%type = boundary_zbc_type

      END SELECT

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE boundary_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (boundary_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%m_index = 0
      this%n_index = 0
      this%type = boundary_no_type

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION boundary_get_modeled_signal(this, a_model, sigma)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: boundary_get_modeled_signal
      TYPE (boundary_class), INTENT(in)       :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  local variables
      INTEGER                           :: m, n, s

!  Start of executable code
      boundary_get_modeled_signal = 0.0
      sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the boundary type
!-------------------------------------------------------------------------------
      FUNCTION boundary_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length)   :: boundary_get_signal_type
      TYPE (boundary_class), INTENT(in) :: this

!  Start of executable code
      boundary_get_signal_type = 'boundary'

      END FUNCTION

      END MODULE
