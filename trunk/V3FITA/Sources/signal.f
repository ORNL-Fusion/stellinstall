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
      INTEGER, PARAMETER :: combination_type_length= 3

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
!-------------------------------------------------------------------------------
!  Signal base class
!
!  Superclass:
!     None
!
!  Subclasses:
!     diagnostic_class
!     geometric_class
!     combination_class
!-------------------------------------------------------------------------------
      TYPE signal_class
         INTEGER                            :: type = signal_no_type
         CHARACTER (len=signal_name_short_length) :: s_name
         CHARACTER (len=signal_name_long_length)  :: l_name
         CHARACTER (len=signal_units_length)      :: units
         REAL (rprec)                       :: observed
         REAL (rprec)                       :: observed_sigma
         REAL (rprec)                       :: weight
         REAL (rprec), DIMENSION(4)         :: modeled
         REAL (rprec), DIMENSION(4)         :: modeled_sigma
         TYPE (diagnostic_class), POINTER   :: diagnostic => null()
         TYPE (geometric_class), POINTER    :: geometric => null()
         TYPE (combination_class), POINTER  :: combination => null()
      END TYPE signal_class

!-------------------------------------------------------------------------------
!  Derived type to hold a pointer to a signal. This is needed because fortran
!  does not allow arrays of pointers directly.
!-------------------------------------------------------------------------------
      TYPE signal_pointer
         TYPE (signal_class), POINTER :: p => null()
      END TYPE

!-------------------------------------------------------------------------------
!  Combination base class
!
!  Superclass:
!     signal_class
!
!  Subclasses:
!     None
!-------------------------------------------------------------------------------
      TYPE combination_class
         INTEGER :: type = combination_none
         TYPE (signal_pointer), DIMENSION(:), POINTER :: signals
         REAL(rprec), DIMENSION(:), POINTER           :: a
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
     &                                     l_name, units, observed,            &
     &                                     sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER        :: signal_construct_diagnostic
      TYPE (diagnostic_class), POINTER, INTENT(in) :: diagnostic_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      CHARACTER(len=signal_units_length), INTENT(in)      :: units
      REAL(rprec), INTENT(in)             :: observed
      REAL(rprec), INTENT(in)             :: sigma
      REAL(rprec), INTENT(in)             :: weight

!  Start of executable code
      ALLOCATE(signal_construct_diagnostic)

      signal_construct_diagnostic%type = signal_geometric_type
      signal_construct_diagnostic%s_name = s_name
      signal_construct_diagnostic%l_name = l_name
      signal_construct_diagnostic%units = units
      signal_construct_diagnostic%observed = observed
      signal_construct_diagnostic%observed_sigma = sigma
      signal_construct_diagnostic%weight = weight
      signal_construct_diagnostic%diagnostic => diagnostic_object
      signal_construct_diagnostic%modeled = 0.0
      signal_construct_diagnostic%modeled_sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a signal containing a geometric.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_geometric(geometric_object, s_name,            &
     &                                    l_name, units, observed,             &
     &                                    sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER       :: signal_construct_geometric
      TYPE (geometric_class), POINTER, INTENT(in) :: geometric_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      CHARACTER(len=signal_units_length), INTENT(in)      :: units
      REAL(rprec), INTENT(in)            :: observed
      REAL(rprec), INTENT(in)            :: sigma
      REAL(rprec), INTENT(in)            :: weight

!  Start of executable code
      ALLOCATE(signal_construct_geometric)

      signal_construct_geometric%type = signal_geometric_type
      signal_construct_geometric%s_name = s_name
      signal_construct_geometric%l_name = l_name
      signal_construct_geometric%units = units
      signal_construct_geometric%observed = observed
      signal_construct_geometric%observed_sigma = sigma
      signal_construct_geometric%weight = weight
      signal_construct_geometric%geometric => geometric_object
      signal_construct_geometric%modeled = 0.0
      signal_construct_geometric%modeled_sigma = 0.0

      END FUNCTION

!-------------------------------------------------------------------------------
!  Construct a signal containing a combination.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_combination(combination_object, s_name,        &
     &                                      l_name, units, observed,           &
     &                                      sigma, weight)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_combination
      TYPE (combination_class), POINTER, INTENT(in)                            &
     &                             :: combination_object
      CHARACTER(len=signal_name_short_length), INTENT(in) :: s_name
      CHARACTER(len=signal_name_long_length), INTENT(in)  :: l_name
      CHARACTER(len=signal_units_length), INTENT(in)      :: units
      REAL(rprec), INTENT(in)              :: observed
      REAL(rprec), INTENT(in)              :: sigma
      REAL(rprec), INTENT(in)              :: weight

!  Start of executable code
      ALLOCATE(signal_construct_combination)

      signal_construct_combination%type = signal_combination_type
      signal_construct_combination%s_name = s_name
      signal_construct_combination%l_name = l_name
      signal_construct_combination%units = units
      signal_construct_combination%observed = observed
      signal_construct_combination%observed_sigma = sigma
      signal_construct_combination%weight = weight
      signal_construct_combination%combination => combination_object
      signal_construct_combination%modeled = 0.0
      signal_construct_combination%modeled_sigma = 0.0

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE signal_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER, INTENT(inout) :: this

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
      this%units = ''
      this%observed = 0.0
      this%observed_sigma = 0.0
      this%weight = 0.0
      this%modeled = 0.0
      this%modeled_sigma = 0.0

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION signal_get_modeled_signal(this, a_model, sigma,                 &
     &                                   use_cache)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: signal_get_modeled_signal
      TYPE (signal_class), INTENT(inout)      :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(:), INTENT(out) :: sigma
      LOGICAL, INTENT(in)                     :: use_cache

!  Start of executable code
      IF (use_cache) THEN
         sigma = this%modeled_sigma
         signal_get_modeled_signal = this%modeled
      ELSE

         SELECT CASE(this%type)
            CASE (signal_diagnostic_type)
               signal_get_modeled_signal =                                     &
     &            diagnostic_get_modeled_signal(this%diagnostic,               &
     &                                          a_model, sigma)

            CASE (signal_geometric_type)
               signal_get_modeled_signal =                                     &
     &            geometric_get_modeled_signal(this%geometric,                 &
     &                                         a_model, sigma)

            CASE (signal_combination_type)
               signal_get_modeled_signal =                                     &
     &            combination_get_modeled_signal(this%combination,             &
     &                                           a_model, sigma)

         END SELECT

         this%modeled = signal_get_modeled_signal
         this%modeled_sigma = sigma

      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the g2 contribution.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_g2(this, a_model)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                       :: signal_get_g2
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model

!  Local Variables
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .false.)
      signal_get_g2 =                                                          &
     &   this%weight*((this%observed - modeled_signal(1))                      &
     &               /(this%observed_sigma - modeled_sigma(1)))**2.0

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
      TYPE (combination_class), POINTER :: combination_construct
      INTEGER, INTENT(in)               :: n_signals
      CHARACTER (len=combination_type_length), INTENT(in)                      &
     &                                  :: combination_type

!  Local Variables
      INTEGER                           :: status

!  Start of executable code
      ALLOCATE(combination_construct)

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
      TYPE (combination_class), POINTER, INTENT(inout) :: this

!  Local Variables
      INTEGER                                          :: i

!  Start of executable code

!  Null all pointers in the signals array. Do not deallocate the signals in the
!  array since this object is not the owner.
      DO i = 1, UBOUND(this%signals, 1)
         this%signals(i)%p => null()
      END DO

      DEALLOCATE(this%signals)
      this%type = combination_none

      DEALLOCATE(this)

      END SUBROUTINE
!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Set the signal
!
!  Add a pointer to the signal to be used for the combination.
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
      this%signals(index)%p => signal

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get the modeled signal
!-------------------------------------------------------------------------------
      FUNCTION combination_get_modeled_signal(this, a_model, sigma)
      USE model
      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)      :: combination_get_modeled_signal
      TYPE (combination_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)       :: a_model
      REAL (rprec), DIMENSION(4), INTENT(out) :: sigma

!  Local Variables
      INTEGER                                 :: i
      REAL (rprec), DIMENSION(4)              :: sigma_local

!  Start of executable code
      combination_get_modeled_signal =                                         &
     &   this%a(1)*signal_get_modeled_signal(this%signals(1)%p, a_model,       &
     &                                       sigma, .true.)

      SELECT CASE(this%type)
         CASE (combination_sum)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal =                                &
     &            combination_get_modeled_signal +                             &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.)
               sigma = sigma + sigma_local
            END DO

         CASE (combination_max)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal = MAX(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.))
               sigma = MAX(sigma, sigma_local)
            END DO

         CASE (combination_min)
            DO i = 2, UBOUND(this%signals, 1)
               combination_get_modeled_signal = MIN(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.))
               sigma = MIN(sigma, sigma_local)
            END DO

      END SELECT

      END FUNCTION

      END MODULE
