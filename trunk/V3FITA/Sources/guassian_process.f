!*******************************************************************************
!>  @file guassian_process.f
!>  @brief Contains module @ref guassian_process.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref guassian_process_class. The
!>  guassian_process contains code to compuet guassian process profiles.
!*******************************************************************************
      MODULE guassian_process
      USE signal

!*******************************************************************************
!  gaussian process module parameters
!*******************************************************************************
!>  Maximum length of the gaussp2 type descritpion string.
      INTEGER, PARAMETER :: gaussp2_type_length = 10

!>  Type descriptor for  the gaussian process type no type.
      INTEGER, PARAMETER :: gaussp2_no_type     = -1
!>  Type descriptor for the gaussian process sxrem type.
      INTEGER, PARAMETER :: gaussp2_sxrem_type  = 0
!>  Type descriptor for the gaussian process te type.
      INTEGER, PARAMETER :: gaussp2_te_type     = 1
!>  Type descriptor for the gaussian process ne type.
      INTEGER, PARAMETER :: gaussp2_ne_type     = 2

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) gaussian process base class
!
!*******************************************************************************

!-------------------------------------------------------------------------------
!>  Base class representing a gaussian process.
!-------------------------------------------------------------------------------
      TYPE gaussp2_class
!>  Type descriptor of the gaussian procees type.
!>  @par Possible values are:
!>  * @ref gaussp2_no_type
!>  * @ref gaussp2_sxrem_type
!>  * @ref gaussp2_te_type
!>  * @ref gaussp2_ne_type
!  Keep for now, but this might not be needed later.
         INTEGER :: type = gaussp2_no_type

!>  Array of @ref signal_pointer
         TYPE (signal_pointer), DIMENSION(:), POINTER ::                       &
     &      signals => null()

!>  Matrix of signal profile.
      REAL (rprec), DIMENSION(:,:), POINTER :: kls => null()
!>  Matrix of signal signal.
      REAL (rprec), DIMENSION(:,:), POINTER :: kll => null()
!>  Work matrix for the cholesky factorization.
      REAL (rprec), DIMENSION(:,:), POINTER :: work => null()
!>  Positive factor to add to the diagonal entries of the matrix before
!>  factoring.
         REAL (rprec)                   :: cholesky_fac

!>  Profile index for sxrem profiles.
         INTEGER  :: profile_index = -1

!>  Offset to add to the log evidence.
         REAL (rprec) :: chi2_offset

!>  Pointer to the profile af array.
         REAL (rprec), DIMENSION(:), POINTER :: fpoints => null()
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref gaussp2_class.
!>
!>  Allocates memory and initializes a @ref gaussp2_class object. This just
!>  allocates the @ref gaussp2_class::signals, @ref gaussp2_class::a
!>  array, and matrixices for the GP computataion.  The indices of these arrays
!>  are set by @ref gaussp2_set_signal, and the matrices are calculated in
!>  @ref gaussp2_get_modeled_signal
!>
!>  @param[in] n_signals     Number of signals in this gaussian process.
!>  @param[in] eq            Pointer to the
!>                           @ref equilibrium::equilibrium_class model.
!>  @param[in] gaussp2_type   The type of profile modeled by the gaussian
!>                           process.
!>  @param[in] profile_index Index of the profile of the sxrem profile.
!>  @param[in] choleskey_fac Value added to diagnoal element to ensure K is PSD.
!>  @returns A pointer to a constructed @ref gaussp2_class object.
!-------------------------------------------------------------------------------
      FUNCTION gaussp2_construct(n_signals, eq, gaussp2_type,                  &
     &                           profile_index, cholesky_fac,                  &
     &                           chi2_offset)
      USE equilibrium

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp2_class), POINTER :: gaussp2_construct
      INTEGER, INTENT(in)                             :: n_signals
      TYPE (equilibrium_class), POINTER, INTENT(in)   :: eq
      INTEGER, INTENT(in)                             :: profile_index
      CHARACTER (len=gaussp2_type_length), INTENT(in) :: gaussp2_type
      REAL (rprec)                                    :: cholesky_fac
      REAL (rprec), INTENT(in)                        :: chi2_offset

!  Local Variables
      INTEGER                                         :: n_points
      REAL (rprec)                                    :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(gaussp2_construct)

      SELECT CASE (gaussp2_type)

         CASE ('sxrem')
            gaussp2_construct%type = gaussp2_sxrem_type
            gaussp2_construct%profile_index = profile_index
            gaussp2_construct%fpoints =>                                       &
     &         equilibrium_get_sxrem_af(eq, profile_index)

         CASE ('te')
            gaussp2_construct%type = gaussp2_te_type
            gaussp2_construct%fpoints => equilibrium_get_te_af(eq)

         CASE ('ne')
            gaussp2_construct%type = gaussp2_ne_type
            gaussp2_construct%fpoints => equilibrium_get_ne_af(eq)

      END SELECT

      n_points = SIZE(gaussp2_construct%fpoints)

      ALLOCATE(gaussp2_construct%signals(n_signals))

      ALLOCATE(gaussp2_construct%kls(n_signals,n_points))
      ALLOCATE(gaussp2_construct%kll(n_signals,n_signals))

      ALLOCATE(gaussp2_construct%work(n_signals,n_points + 1))

      gaussp2_construct%cholesky_fac = cholesky_fac
      gaussp2_construct%chi2_offset = chi2_offset

      CALL profiler_set_stop_time('gaussp2_construct', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref gaussp2_class object.
!>
!>  Deallocates memory and uninitializes a @ref gaussp2_class object.
!>
!>  @param[inout] this A @ref gaussp2_class instance.
!>  @note The destructors of the signals in the signal array are not called her
!>  because the gaussian process signals do not own the memory for those signals.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp2_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp2_class), POINTER :: this

!  Local Variables
      INTEGER                       :: i

!  Start of executable code

!  Null all pointers in the signals array. Do not deallocate the signals in the
!  array since this object is not the owner.
      DO i = 1, SIZE(this%signals)
         this%signals(i)%p => null()
      END DO

      IF (ASSOCIATED(this%signals)) THEN
         DEALLOCATE(this%signals)
         this%signals => null()
      END IF

      IF (ASSOCIATED(this%kls)) THEN
         DEALLOCATE(this%kls)
         this%kls => null()
      END IF

      IF (ASSOCIATED(this%kll)) THEN
         DEALLOCATE(this%kll)
         this%kll => null()
      END IF

      IF (ASSOCIATED(this%work)) THEN
         DEALLOCATE(this%work)
         this%work => null()
      END IF

      this%fpoints => null()

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Set the object and coefficient for an index.
!>
!>  @param[inout] this   A @ref gaussp2_class instance.
!>  @param[in]    signal A @ref signal_class instance.
!>  @param[in]    a      A weight for the signal
!>  @param[in]    index  The index to place the gaussp2 @ref signal.
!>  @note Only signals that were created before this gaussian process
!>  can be added. Otherwise it will creat an infinitely recursive loop in
!>  @ref gaussp2_get_modeled_signal
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp2_set_signal(this, signal, index)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (gaussp2_class), INTENT(inout) :: this
      TYPE (signal_class), POINTER        :: signal
      INTEGER, INTENT(in)                 :: index

!  local Variables
      REAL (rprec)                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      this%signals(index)%p => signal

      CALL profiler_set_stop_time('gaussp2_set_signal', start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set the object and coefficient for an index.
!>
!>  This method uses a gaussian process to set a model using the signals pointed
!>  in the  @ref gaussp_class::signals array. Hyper parameters are determined
!>  by maximizing the evidence.
!>
!>  @param[in]    this         A @ref gaussp_class instance.
!>  @param[inout] a_model      A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE gaussp2_set_profile(this, a_model)
      USE model
      USE model_state

!  Declare Arguments
      TYPE (gaussp_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)  :: a_model

!  Local Variables
      INTEGER                            :: i
      INTEGER                            :: j
      REAL (rprec)                       :: start_time
      INTEGER                            :: ierr
!      REAL (rprec)                       :: log_evidence
      REAL (rprec)                       :: evidence
      REAL (rprec)                       :: devidence

!  local parameters
      REAL(rprec), PARAMETER             :: tolerance = 1.0E-8

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (BTEST(a_model%state_flags, model_state_vmec_flag)   .or.             &
     &    BTEST(a_model%state_flags, model_state_siesta_flag) .or.             &
     &    BTEST(a_model%state_flags, model_state_shift_flag)  .or.             &
     &    BTEST(a_model%state_flags, model_state_signal_flag) .or.             &
     &    (this%type .eq. gaussp_ne_type)    .and.                             &
     &     BTEST(a_model%state_flags, model_state_ne_flag)    .or.             &
     &    (this%type .eq. gaussp_te_type)    .and.                             &
     &     BTEST(a_model%state_flags, model_state_te_flag)    .or.             &
     &    (this%type .eq. gaussp_sxrem_type) .and.                             &
     &     BTEST(a_model%state_flags, model_state_sxrem_flag +                 &
     &                                this%profile_index - 1)) THEN

!  Use gradient ascent to fine the optimal hyper parameters by maximizing the
!  the evidence.
         DO WHILE (devidence > tolerance)

         END DO

      END IF

      CALL profiler_set_stop_time('gaussp_get_modeled_signal',                 &
     &                            start_time)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculates the evidence.
!>
!>  This method finds a gaussian process evidence. The matrices computed can be
!>  later used to set the profile.
!>
!>  @param[in]    this    A @ref gaussp_class instance.
!>  @param[inout] a_model A @ref model instance.
!>  @returns The log evidence for the current hyper parameters.
!-------------------------------------------------------------------------------
      FUNCTION gaussp2_get_evidence(this, a_model)
      USE model
      USE stel_constants, ONLY: twopi

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4)          :: gaussp2_get_evidence
      TYPE (gaussp2_class), INTENT(inout) :: this
      TYPE (model_class), INTENT(inout)   :: a_model

!  Local Variables
      INTEGER                             :: i
      INTEGER                             :: j
      REAL (rprec)                        :: start_time
      INTEGER                             :: ierr

!  Local Parameters
      REAL (rprec), PARAMETER             :: log2pi = LOG(twopi)

!  Start of executable code
      start_time = profiler_get_start_time()

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j)

!$OMP DO
!$OMP& SCHEDULE(DYNAMIC)
!  FIXME: There are duplicate matrix elements here. We could avoid unnecessary
!         computation.
      DO j = 1, SIZE(this%fpoints)
         DO i = 1, SIZE(this%signals)
            this%kls(i,j) = signal_get_gp(this%signals(i)%p, a_model, j)
            this%work(i,j + 1) = this%kls(i,j)
         END DO
      END DO
!$OMP END DO

!$OMP DO
!$OMP& SCHEDULE(DYNAMIC)
      DO j = 1, SIZE(this%signals)
         DO i = j, SIZE(this%signals)
            this%kll(i,j) = signal_get_gp(this%signals(i)%p, a_model,          &
     &                                    this%signals(j)%p)
            IF (j .ne. i) THEN
               this%kll(j,i) = this%kll(i,j)
            END IF
         END DO
         this%kll(j,j) = this%kll(j,j)                                         &
     &                 + signal_get_sigma2(this%signals(j)%p)                  &
     &                 + this%cholesky_fac
         this%work(j,1) = signal_get_observed_signal(this%signals(j)%p,        &
     &                                               a_model)
      END DO
!$OMP END DO
!$OMP END PARALLEL

!  CALL LAPACK to do Cholesky factorization of the A=(K_LL + Sigma_y)
      CALL DPOTRF('L', SIZE(this%signals), this%kll,                           &
     &            SIZE(this%signals), ierr)

!  The elements in the cholesky matrix for j > i are meaning less.
      IF (ierr .lt. 0) THEN
         CALL err_fatal('gaussp_get_modeled_signal: DROTRF ' //                &
     &                  'cannot factor matrix. ierr is negative')
      ELSE IF (ierr .gt. 0) THEN
         CALL err_fatal('gaussp_get_modeled_signal: DROTRF ' //                &
     &                  'cannot factor matrix. ierr is '     //                &
     &                  'positive. Try increasing gp_cholesky_fact')
      END IF

!  IF statement correctly throws an error if I cannot factor the matrix!

!  CALL LAPACK to solve A^-1y and A^-1 K_LS simultaneously. The work matrix
!  stores the RHS on input, and returns the solution matrix.
      CALL DPOTRS('L', SIZE(this%signals), SIZE(this%fpoints) + 1,             &
     &            this%kll, SIZE(this%signals), this%work,                     &
     &            SIZE(this%signals), ierr)
      IF (ierr .lt. 0) THEN
         CALL err_fatal('gaussp_get_modeled_signal: DROTRS ' //                &
     &                  'failed to solve the equation')
      END IF

!  Calculate log evidence.
!
!    P(y|sigma,x) = 1/Sqrt((2Pi)^N*|A|)Exp(-y^T.A^-1.y/2)                    (1)
!
!  Taking the ln of the evidence results in.
!
!    ln(P) = -N/2*ln(2Pi) - ln(|A|)/2 - y^T.A^-1.y/2                         (2)
!
!  The determinant of A is found by the Cholesky factorization and equals the
!  product of the squares of the diagonal.
!
!    |A| = Prod (A_ii)^2                                                     (3)
!
!  Taking the ln of this results in.
!
!    ln(|A|) = Sum 2*ln(A_ii)                                                (4)
!

      gaussp2_get_evidence = -0.5*SIZE(this%signals)*log2pi

      DO j = 1, SIZE(this%signals)
         gaussp2_get_evidence = gaussp2_get_evidence                           &
     &                        - LOG(this%kll(j,j))                             &
     &                        - 0.5*signal_get_observed_signal(                &
     &                                 this%signals(j)%p, a_model)*            &
     &                          this%work(j,1)
      END DO

      CALL profiler_set_stop_time('gaussp2_get_evidence', start_time)

      END FUNCTION

      END MODULE
