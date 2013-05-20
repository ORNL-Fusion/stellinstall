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
      FUNCTION boundary_get_modeled_signal(this, a_model)
      USE model_T

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)        :: boundary_get_modeled_signal
      TYPE (boundary_class), INTENT(in) :: this
      TYPE (model), INTENT(inout)       :: a_model

!  local variables
      INTEGER                           :: m, n, s

!  Start of executable code
      IF (.not. a_model%eqstate%l_def_aux1) THEN
         CALL eq_aux_1_define(a_model%eqstate)
      END IF

      s = a_model%eqstate%fixp%ns_array(1)
      m = this%m_index
      n = this%n_index

      SELECT CASE (this%type)
         CASE (boundary_rbc_type)
            boundary_get_modeled_signal(1) =                                   &
     &         a_model%eqstate%aux1%rmnc(m*n + m + 1, s)

		 CASE (boundary_zbs_type)
            boundary_get_modeled_signal(1) =                                   &
     &         a_model%eqstate%aux1%zmns(m*n + m + 1, s)

         CASE (boundary_rbs_type)
            boundary_get_modeled_signal(1) =                                   &
     &         a_model%eqstate%aux1%rmns(m*n + m + 1, s)

         CASE (boundary_zbc_type)
            boundary_get_modeled_signal(1) =                                   &
     &         a_model%eqstate%aux1%zmnc(m*n + m + 1, s)

      END SELECT

      END FUNCTION

      END MODULE
