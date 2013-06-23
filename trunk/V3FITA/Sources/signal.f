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
     &                    signal_construct_combination,                        &
     &                    signal_construct_diagnostic_netcdf
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
      CHARACTER(len=*), INTENT(in)        :: s_name
      CHARACTER(len=*), INTENT(in)        :: l_name
      CHARACTER(len=*), INTENT(in)        :: units
      REAL(rprec), INTENT(in)             :: observed
      REAL(rprec), INTENT(in)             :: sigma
      REAL(rprec), INTENT(in)             :: weight

!  Start of executable code
      ALLOCATE(signal_construct_diagnostic)

      signal_construct_diagnostic%type = signal_diagnostic_type
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
      CHARACTER(len=*), INTENT(in)       :: s_name
      CHARACTER(len=*), INTENT(in)       :: l_name
      CHARACTER(len=*), INTENT(in)       :: units
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
      CHARACTER(len=*), INTENT(in) :: s_name
      CHARACTER(len=*), INTENT(in) :: l_name
      CHARACTER(len=*), INTENT(in) :: units
      REAL(rprec), INTENT(in)      :: observed
      REAL(rprec), INTENT(in)      :: sigma
      REAL(rprec), INTENT(in)      :: weight

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

!-------------------------------------------------------------------------------
!  Construct a signal containing a diagnostic from a netcdf file.
!-------------------------------------------------------------------------------
      FUNCTION signal_construct_diagnostic_netcdf(diagnostic_object,           &
     &                                            mdsig_iou, observed,         &
     &                                            sigma, weight)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), POINTER :: signal_construct_diagnostic_netcdf
      TYPE (diagnostic_class), POINTER, INTENT(in) :: diagnostic_object
      INTEGER, INTENT(in)                      :: mdsig_iou
      REAL(rprec), INTENT(in)                  :: observed
      REAL(rprec), INTENT(in)                  :: sigma
      REAL(rprec), INTENT(in)                  :: weight

!  local variables
      CHARACTER (len=signal_name_short_length) :: s_name
      CHARACTER (len=signal_name_long_length)  :: l_name
      CHARACTER (len=signal_units_length)      :: units

!  Start of executable code
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_s_name', s_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_l_name', l_name)
      CALL cdf_read(mdsig_iou, 'diagnostic_desc_units', units)

      signal_construct_diagnostic_netcdf =>                                     &
     &   signal_construct(diagnostic_object, s_name, l_name, units,             &
     &                   observed, sigma, weight)

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
      REAL (rprec)                       :: e
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .false.)

! g^2 = W*(O - M)^2/(sigma_o^2 + sigma_m^2) = e*e
      e = signal_get_e(this, a_model)
      signal_get_g2 = e*e

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the e contribution.
!-------------------------------------------------------------------------------
      FUNCTION signal_get_e(this, a_model)
      USE model

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                       :: signal_get_e
      TYPE (signal_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model

!  Local Variables
      REAL (rprec), DIMENSION(4)         :: modeled_signal
      REAL (rprec), DIMENSION(4)         :: modeled_sigma

!  Start of executable code
      modeled_signal = signal_get_modeled_signal(this, a_model,                &
     &                                           modeled_sigma, .false.)

! e = SQRT(W)*(O - M)/SQRT(sigma_o^2 + sigma_m^2)
      signal_get_e = (this%observed - modeled_signal(1))                       &
     &             * SQRT(this%weight/(this%observed_sigma**2.0 +              &
     &                                 modeled_sigma(1)**2.0))

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the signal type
!-------------------------------------------------------------------------------
      FUNCTION signal_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: signal_get_signal_type
      TYPE (signal_class), INTENT(in)  :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_signal_type =                                           &
     &         TRIM(diagnostic_get_signal_type(this%diagnostic))

         CASE (signal_geometric_type)
            signal_get_signal_type =                                           &
     &         TRIM(geometric_get_signal_type(this%geometric))

         CASE (signal_combination_type)
            signal_get_signal_type =                                           &
     &         TRIM(combination_get_signal_type(this%combination))

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the signal header
!-------------------------------------------------------------------------------
      FUNCTION signal_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(3) ::                        &
     &   signal_get_header
      TYPE (signal_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE(this%type)

         CASE (signal_diagnostic_type)
            signal_get_header = diagnostic_get_header(this%diagnostic)

         CASE (signal_geometric_type)
            signal_get_header = geometric_get_header(this%geometric)

         CASE (signal_combination_type)
            signal_get_header = combination_get_header(this%combination)

      END SELECT

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Generate a short name with an index
!-------------------------------------------------------------------------------
      FUNCTION signal_make_short_name(name, index)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=signal_name_short_length) :: signal_make_short_name
      CHARACTER (len=*), INTENT(in)            :: name
      INTEGER, INTENT(in)                      :: index

!  Start of executable code
!  Find the number of digits in the index
      IF (index .lt. 10) THEN
         WRITE(signal_make_short_name,'(a,i1)') TRIM(name), index
      ELSE IF (index .lt. 100) THEN
         WRITE(signal_make_short_name,'(a,i2)') TRIM(name), index
      ELSE IF (index .lt. 1000) THEN
         WRITE(signal_make_short_name,'(a,i3)') TRIM(name), index
      ELSE
         WRITE(signal_make_short_name,'(a,i4)') TRIM(name), index
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!  Write out the signal
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou
      INTEGER, INTENT(in)             :: index

!  Start of executable code
!  FIXME: It is too convoluted to keep the recout identical at this point.
!  Write out the model 1 and 2 for every signal. Eventually all for should be
!  written out.
      WRITE (iou,1000) index, signal_get_signal_type(this),                    &
     &                 this%s_name, this%observed, this%observed_sigma,        &
     &                 this%modeled

1000  FORMAT(i4,2x,a23,2x,a20,2(2x,es12.5),4(2x,es12.5))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Write out the signal header
!-------------------------------------------------------------------------------
      SUBROUTINE signal_write_header(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,1000) signal_get_header(this)

1000  FORMAT(3x,'#',2x,'type',21x,'s_name',16x,'observed',6x,'sigma',          &
     &       9x,'model',7x,2(2x,a12),3(2x,a12))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Write out the signal auxiliary information.
!-------------------------------------------------------------------------------
      SUBROUTINE signals_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_class), INTENT(in) :: this
      INTEGER, INTENT(in)             :: iou
      INTEGER, INTENT(in)             :: index

!  Start of executable code
      SELECT CASE (this%type)

         CASE (signal_combination_type)
            CALL combination_write_auxiliary(this%combination,                 &
     &                                       iou, index)

      END SELECT

      END SUBROUTINE

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
      DO i = 1, SIZE(this%signals)
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
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal =                                &
     &            combination_get_modeled_signal +                             &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.)
               sigma = sigma + sigma_local
            END DO

         CASE (combination_max)
            DO i = 2, SIZE(this%signals)
               combination_get_modeled_signal = MAX(                           &
     &            combination_get_modeled_signal,                              &
     &            this%a(i)*signal_get_modeled_signal(this%signals(i)%p,       &
     &                                                a_model,                 &
     &                                                sigma_local,             &
     &                                                .true.))
               sigma = MAX(sigma, sigma_local)
            END DO

         CASE (combination_min)
            DO i = 2, SIZE(this%signals)
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

!-------------------------------------------------------------------------------
!  Get the combination type
!-------------------------------------------------------------------------------
      FUNCTION combination_get_signal_type(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length) :: combination_get_signal_type
      TYPE (combination_class), INTENT(in) :: this

!  Start of executable code
      SELECT CASE (this%type)
         CASE (combination_sum)
            combination_get_signal_type = 'coosig sum'

         CASE (combination_max)
            combination_get_signal_type = 'coosig max'

         CASE (combination_min)
            combination_get_signal_type = 'coosig min'

      END SELECT

      END FUNCTION

!-------------------------------------------------------------------------------
!  Get the combination header
!-------------------------------------------------------------------------------
      FUNCTION combination_get_header(this)
      USE data_parameters

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=data_name_length), DIMENSION(3) ::                        &
     &   combination_get_header
      TYPE (combination_class), INTENT(in) :: this

!  Start of executable code
      combination_get_header(1) = 'model(2)'
      combination_get_header(2) = 'model(3)'
      combination_get_header(3) = 'model(4)'

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the combination auxiliary information.
!-------------------------------------------------------------------------------
      SUBROUTINE combination_write_auxiliary(this, iou, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (combination_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: iou
      INTEGER, INTENT(in)                  :: index

!  local variables
      INTEGER                              :: i

!  Start of executable code
      WRITE (iou,*)
      WRITE (iou,1000) index, combination_get_signal_type(this)
      WRITE (iou,1001)

      DO i = 1, SIZE(this%signals)
         WRITE (iou,1002) i, this%signals(i)%p%s_name, this%a(i)
      END DO

1000  FORMAT('Signal',1x,i4,1x,                                                &
     &       'is a combination of other signals, type: ',a)
1001  FORMAT('term #',2x,'s_name',16x,'Coefficient')
1002  FORMAT(2x,i4,2x,a20,2x,es12.5)

      END SUBROUTINE

      END MODULE
