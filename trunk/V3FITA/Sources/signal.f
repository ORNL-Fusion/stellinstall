!*******************************************************************************
!  File signal.f
!  Contains module signal
!  Defines the base class of the type signal:
!
!  Super Class: None
!  Sub Classes: diagnostic, geometric, combination
!
!  ---------------------------------------------------------------------------
!
!  Defines the base class of the type combination:
!
!  Super Class: signal
!  Sub Classes: None
!
!  NOTE: combination is defined in the signal module to avoid a circular
!  module dependency.
!
!*******************************************************************************

      MODULE signal

      USE stel_kinds, only: rprec
      USE diagnostic
      USE geometric

      IMPLICIT NONE

!*******************************************************************************
!  signal module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: signal_name_short_length = 30
      INTEGER, PARAMETER :: signal_name_long_length = 80
      INTEGER, PARAMETER :: signal_units_length = 30

      INTEGER, PARAMETER :: signal_no_type = -1
      INTEGER, PARAMETER :: signal_diagnostic_type = 0
      INTEGER, PARAMETER :: signal_geometric_type = 1
      INTEGER, PARAMETER :: signal_combination_type = 2

!*******************************************************************************
!  combination module parameters
!*******************************************************************************
      INTEGER, PARAMETER :: combination_none = -1
      INTEGER, PARAMETER :: combination_sum = 0
      INTEGER, PARAMETER :: combination_max = 1
      INTEGER, PARAMETER :: combination_min = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) signal base class
!  2) signal pointer type
!  3) combination base class
!
!*******************************************************************************

      TYPE signal_class
         INTEGER                           :: type = signal_no_type
         CHARACTER(len=signal_name_short_length) :: s_name
         CHARACTER(len=signal_name_long_length)  :: l_name
         REAL(rprec)                       :: observed
         REAL(rprec)                       :: sigma
         REAL(rprec)                       :: weight
         TYPE(diagnostic_class), POINTER   :: diagnostic => null()
         TYPE(geometric_class), POINTER    :: geometric => null()
         TYPE(combination_class), POINTER  :: combination => null()
      END TYPE signal_class

      TYPE signal_pointer
         TYPE (signal_class), POINTER :: p => null()
      END TYPE

      TYPE combination_class
         INTEGER :: type = combination_none
         TYPE (signal_pointer), DIMENSION(:), ALLOCATABLE :: signals
         REAL(rprec), DIMENSION(:), ALLOCATABLE :: a
      END TYPE combination_class

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE signal_construct
         MODULE PROCEDURE signal_construct_diagnostic,                         &
     &                    signal_construct_geometric,                          &
     &                    signal_construct_combination                         
      END INTERFACE

      CONTAINS 
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a signal containing a diagnostic.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_diagnostic(diagnostic_object, s_name,          &
     &                                     l_name, observed, sigma,            &
     &                                     weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class)                 :: signal_construct_diagnostic
      TYPE (diagnostic_class), INTENT(in) :: diagnostic_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      REAL(rprec), INTENT(in)             :: observed
      REAL(rprec), INTENT(in)             :: sigma
      REAL(rprec), INTENT(in)             :: weight

!  Start of executable code
      signal_construct_diagnostic%type = signal_geometric_type
      signal_construct_diagnostic%s_name = s_name
      signal_construct_diagnostic%l_name = l_name
      signal_construct_diagnostic%observed = observed
      signal_construct_diagnostic%sigma = sigma
      signal_construct_diagnostic%weight = weight
      signal_construct_diagnostic%diagnostic = diagnostic_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a signal containing a geometric.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_geometric(geometric_object, s_name,            &
     &                                    l_name, observed, sigma,             &
     &                                    weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class)                :: signal_construct_geometric
      TYPE (geometric_class), INTENT(in) :: geometric_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  Start of executable code
      signal_construct_geometric%type = signal_geometric_type
      signal_construct_geometric%s_name = s_name
      signal_construct_geometric%l_name = l_name
      signal_construct_geometric%observed = observed
      signal_construct_geometric%sigma = sigma
      signal_construct_geometric%weight = weight
      signal_construct_geometric%geometric = geometric_object

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a signal containing a combination.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_combination(combination_object, s_name,        &
     &                                      l_name, observed, sigma,           &
     &                                      weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class)               :: signal_construct_combination
      TYPE (combination_class), INTENT(in) :: combination_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      REAL(rprec), INTENT(in)              :: observed
      REAL(rprec), INTENT(in)              :: sigma
      REAL(rprec), INTENT(in)              :: weight

!  Start of executable code
      signal_construct_combination%type = signal_combination_type
      signal_construct_combination%s_name = s_name
      signal_construct_combination%l_name = l_name
      signal_construct_combination%observed = observed
      signal_construct_combination%sigma = sigma
      signal_construct_combination%weight = weight
      signal_construct_combination%combination = combination_object

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE signal_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(inout) :: this

!  Start of executable code
      IF (ASSOCIATED(this%diagnostic)) THEN
         CALL diagnostic_destruct(this%diagnostic)
         this%diagnostic => null()
      END IF

      IF (ASSOCIATED(this%geometric)) THEN
         CALL geometric_destruct(this%geometric)
         this%geometric => null()
      END IF

      IF (ASSOCIATED(this%combination)) THEN
         CALL combination_destruct(this%combination)
         this%combination => null()
      END IF

      this%type = signal_no_type
      this%s_name = ''
      this%l_name = ''
      this%observed = 0.0
      this%sigma = 0.0
      this%weight = 0.0

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION signal_get_modeled_signal(this, a_model)
      USE model_T
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2)         :: signal_get_modeled_signal
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model), INTENT(inout)        :: a_model

!  Start of executable code
      SELECT CASE(this%type)
         CASE (signal_diagnostic_type)
            signal_get_modeled_signal =                                        &
     &         diagnostic_get_modeled_signal(this%diagnostic, a_model)

         CASE (signal_geometric_type)
            signal_get_modeled_signal =                                        &
     &         geometric_get_modeled_signal(this%geometric, a_model)

         CASE (signal_combination_type)
            signal_get_modeled_signal =                                        &
     &         combination_get_modeled_signal(this%combination, a_model)

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the g2 contribution.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_g2(this, a_model)
      USE model_T
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                       :: signal_get_g2
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model), INTENT(inout)        :: a_model

!  Local Variables
      REAL (rprec), DIMENSION(2)         :: modeled_signal

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model)
      signal_get_g2 =                                                          &
     &   this%weight*((this%observed - modeled_signal(1))                      &
     &                /this%sigma)**2.0

      END FUNCTION

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!  Combination Signals
!
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a combination object.
!-------------------------------------------------------------------------------
      FUNCTION combination_construct(n_signals, combination_type)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class)      :: combination_construct
      INTEGER, INTENT(in)           :: n_signals
      CHARACTER (len=3), INTENT(in) :: combination_type

!  Local Variables
      INTEGER                       :: status

!  Start of executable code
      ALLOCATE(combination_construct%signals(n_signals),                       &
     &         combination_construct%a(n_signals), STAT=status)

      SELECT CASE (combination_type)
         CASE ('sum')
            combination_construct%type = combination_sum

         CASE ('max')
            combination_construct%type = combination_max

         CASE ('min')
            combination_construct%type = combination_min

      END SELECT

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE combination_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(inout) :: this

!  Local Variables
      INTEGER                                 :: i

!  Start of executable code

!  Null all pointers in the signals array.
      DO i = 1, UBOUND(this%signals, 1)
         this%signals(i)%p => null()
      END DO

      DEALLOCATE(this%signals)
      this%type = combination_none

      END SUBROUTINE
!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Set the signal
!-------------------------------------------------------------------------------
      SUBROUTINE combination_set_signal(this, signal, a, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(inout)  :: this
      TYPE (signal_class), INTENT(in), POINTER :: signal
      INTEGER, INTENT(in)                      :: index
      REAL (rprec), INTENT(in)                 :: a

!  Local Variables

!  Start of executable code
      this%a(index) = a
      this%signals(index)%p = signal

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION combination_get_modeled_signal(this, a_model)
      USE model_T
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(2)      :: combination_get_modeled_signal
      TYPE (combination_class), INTENT(inout) :: this
      TYPE (model), INTENT(inout)             :: a_model

!  Local Variables
      INTEGER                                 :: i

!  Start of executable code
      combination_get_modeled_signal =                                         &
     &   this%a(1)*signal_get_modeled_signal(this%signals(1)%p, a_model)

      SELECT CASE(this%type)
         CASE (combination_sum)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal =                                &
     &            combination_get_modeled_signal +                             &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model)
            END DO

         CASE (combination_max)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal = MAX(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model))
            END DO

         CASE (combination_min)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal = MIN(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model))
            END DO

      END SELECT

      END FUNCTION

      END MODULE
