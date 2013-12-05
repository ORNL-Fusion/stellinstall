!*******************************************************************************
!>  @file reconstruction.f
!>  @brief Contains module @ref reconstruction.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref reconstruction_class. This class
!>  contains the minimization algorithm.
!*******************************************************************************

      MODULE reconstruction
      USE stel_kinds
      USE param
      USE signal

      IMPLICIT NONE

!*******************************************************************************
!  reconstruction module parameters
!*******************************************************************************
!>  No reconstruction step type.
      INTEGER, PARAMETER :: reconstruction_no_step_type = -1
!>  Straight line reconstruction step type.
      INTEGER, PARAMETER :: reconstruction_sl_step_type = 0
!>  Levenberg-Marquardt reconstruction step type.
      INTEGER, PARAMETER :: reconstruction_lm_step_type = 1
!>  Segmented path reconstruction step type.
      INTEGER, PARAMETER :: reconstruction_seg_step_type = 2

!>  Maximum number of reconstruction step attemps.
      INTEGER, PARAMETER :: reconstruction_max_step_try = 5

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) reconstruction base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing containing all the data needed to reconstruct a
!>  model.
!-------------------------------------------------------------------------------
      TYPE reconstruction_class
!>  Type descriptor of the boundry type for the lower(1) and upper(2) ranges.
!>  @par Possible values are:
!>  * @ref reconstruction_no_step_type
!>  * @ref reconstruction_sl_step_type
!>  * @ref reconstruction_lm_step_type
!>  * @ref reconstruction_seg_step_type
         INTEGER :: step_type = reconstruction_no_step_type
!>  Maximum number of reconstruction step.
!>  @see v3fit_input::nrstep
         REAL (rprec)                          :: step_max

!>  Counter to track the step number.
         INTEGER                               :: current_step = 0
!>  Array of all e vectors for each reconstruction step.
!>  @see signal::signal_get_e
         REAL (rprec), DIMENSION(:,:), POINTER :: e => null()
!>  Array of all f vectors for each reconstruction step. The f vector contains
!>  the initial value of the derived parameter.
         REAL (rprec), DIMENSION(:,:), POINTER :: f => null()
!>  The normalized jacobian for the current step.
         REAL (rprec), DIMENSION(:,:), POINTER :: jacobian => null()
!>  The normalized derived jacobian for the current step.
         REAL (rprec), DIMENSION(:,:), POINTER ::                              &
     &      derived_jacobian => null()
!>  The normalized hessian for the current step.
         REAL (rprec), DIMENSION(:,:), POINTER :: hessian => null()
!>  The normalized gradient for the current step.
         REAL (rprec), DIMENSION(:), POINTER   :: gradient => null()
!>  The normalized parameter step for each sigular value.
         REAL (rprec), DIMENSION(:,:), POINTER :: delta_a => null()
!>  The matrix of singular value decomposition singular values.
         REAL (rprec), DIMENSION(:), POINTER   :: j_svd_w => null()
!>  The U matrix of singular value decomposition.
         REAL (rprec), DIMENSION(:,:), POINTER :: j_svd_u => null()
!>  The tansposed V matrix of singular value decomposition.
         REAL (rprec), DIMENSION(:,:), POINTER :: j_svd_vt => null()
!>  The normalized step length for each sigular value.
         REAL (rprec), DIMENSION(:), POINTER   :: delta_a_len => null()

!>  The expected g^2 for each reconstruction step.
         REAL (rprec), DIMENSION(:), POINTER   :: exp_g2 => null()
!>  The normalized step size for each reconstruction step.
         REAL (rprec), DIMENSION(:), POINTER   :: step_size => null()
!>  The number of sigular values used for each reconstruction step.
         INTEGER, DIMENSION(:), POINTER        :: num_sv => null()

!  Cut off value controls.
!>  Cutoff value for relative singular values.
!>  @see v3fit_input::cut_svd
         REAL (rprec) :: cut_svd
!>  Cutoff value for expected step efficiency.
!>  @see v3fit_input::cut_eff
         REAL (rprec) :: cut_eff
!>  Cutoff value for expected marginal step efficiency.
!>  @see v3fit_input::cut_marg_eff
         REAL (rprec) :: cut_marg_eff
!>  Cutoff value for expected step size.
!>  @see v3fit_input::cut_delta_a
         REAL (rprec) :: cut_delta_a
!>  Cutoff value for expected change in g^2.
!>  @see v3fit_input::cut_dg2
         REAL (rprec) :: cut_dg2
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref reconstruction_class object.
!>
!>  Allocates memory and initializes a @ref param_class object.
!>
!>  @param[in] num_steps              Number of reconstruction steps.
!>                                    @ref signal::signal_get_e
!>  @param[in] num_signals            Number of signals.
!>  @param[in] num_derived_parameters Number of drived parameters.
!>  @param[in] num_parameters         Number of reconstruction parameters.
!>  @param[in] step_type              Reconstruction step type.
!>                                    @ref signal::step_type
!>  @param[in] step_max               Reconstruction step size.
!>                                    @ref signal::astep_max
!>  @param[in] cut_svd                Cutoff value for relative singular values.
!>                                    @ref v3fit_input::cut_svd
!>  @param[in] cut_eff                Cutoff value for expected step efficiency.
!>                                    @ref v3fit_input::cut_eff
!>  @param[in] cut_marg_eff           Cutoff value for expected marginal step
!>                                    efficiency. @ref v3fit_input::cut_marg_eff
!>  @param[in] cut_delta_a            Cutoff value for expected step size.
!>                                    @ref v3fit_input::cut_delta_a
!>  @param[in] cut_dg2                Cutoff value for expected change in g^2.
!>                                    @ref v3fit_input::cut_dg2
!>  @returns A pointer to a constructed @ref reconstruction_class object.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_construct(num_steps, num_signals,                &
     &                                  num_derived_parameters,                &
     &                                  num_parameters, step_type,             &
     &                                  step_max, cut_svd, cut_eff,            &
     &                                  cut_marg_eff, cut_delta_a,             &
     &                                  cut_dg2)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), POINTER :: reconstruction_construct
      INTEGER, INTENT(in)                  :: num_steps
      INTEGER, INTENT(in)                  :: num_signals
      INTEGER, INTENT(in)                  :: num_derived_parameters
      INTEGER, INTENT(in)                  :: num_parameters
      CHARACTER (len=*), INTENT(in)        :: step_type
      REAL (rprec), INTENT(in)             :: step_max
      REAL (rprec), INTENT(in)             :: cut_svd
      REAL (rprec), INTENT(in)             :: cut_eff
      REAL (rprec), INTENT(in)             :: cut_marg_eff
      REAL (rprec), INTENT(in)             :: cut_delta_a
      REAL (rprec), INTENT(in)             :: cut_dg2

!  local variables
      INTEGER                              :: num_svd

!  Start of executable code
      ALLOCATE(reconstruction_construct)

      SELECT CASE (TRIM(step_type))

         CASE ('lm')
            reconstruction_construct%step_type =                               &
     &         reconstruction_lm_step_type

         CASE ('seg')
            reconstruction_construct%step_type =                               &
     &         reconstruction_seg_step_type

         CASE DEFAULT
            reconstruction_construct%step_type =                               &
     &         reconstruction_sl_step_type

      END SELECT

      reconstruction_construct%step_max = step_max

      reconstruction_construct%cut_svd = cut_svd
      reconstruction_construct%cut_eff = cut_eff
      reconstruction_construct%cut_marg_eff = cut_marg_eff
      reconstruction_construct%cut_delta_a = cut_delta_a
      reconstruction_construct%cut_dg2 = cut_dg2

      num_svd = MIN(num_signals, num_parameters)

      ALLOCATE(reconstruction_construct%e(0:num_steps,num_signals))
      ALLOCATE(reconstruction_construct%f(0:num_steps,                         &
     &                                    num_derived_parameters))
      ALLOCATE(reconstruction_construct%jacobian(num_signals,                  &
     &                                           num_parameters))
      ALLOCATE(reconstruction_construct%derived_jacobian(                      &
     &            num_derived_parameters,num_parameters))
      ALLOCATE(reconstruction_construct%hessian(num_parameters,                &
     &                                          num_parameters))
      ALLOCATE(reconstruction_construct%gradient(num_parameters))
      ALLOCATE(reconstruction_construct%delta_a(0:num_svd,                     &
     &                                          num_parameters))
      ALLOCATE(reconstruction_construct%j_svd_w(num_svd))
      ALLOCATE(reconstruction_construct%j_svd_u(num_signals,                   &
     &                                          num_signals))
      ALLOCATE(reconstruction_construct%j_svd_vt(num_parameters,               &
     &                                           num_parameters))
      ALLOCATE(reconstruction_construct%delta_a_len(0:num_svd))
      reconstruction_construct%e = 0.0

      ALLOCATE(reconstruction_construct%exp_g2(num_steps))
      ALLOCATE(reconstruction_construct%step_size(num_steps))
      ALLOCATE(reconstruction_construct%num_sv(num_steps))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref reconstruction_class object.
!>
!>  Deallocates memory and uninitializes a @ref reconstruction_class object.
!>
!>  @param[inout] this A @ref reconstruction_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), POINTER, INTENT(inout) :: this

!  Start of executable code
      this%step_type = reconstruction_no_step_type

      this%cut_svd = 0.0
      this%cut_eff = 0.0
      this%cut_marg_eff = 0.0
      this%cut_delta_a = 0.0
      this%cut_dg2 = 0.0

      this%current_step = 0

      IF (ASSOCIATED(this%e)) THEN
         DEALLOCATE(this%e)
         this%e => null()
      END IF

      IF (ASSOCIATED(this%f)) THEN
         DEALLOCATE(this%f)
         this%f => null()
      END IF

      IF (ASSOCIATED(this%jacobian)) THEN
         DEALLOCATE(this%jacobian)
         this%jacobian => null()
      END IF

      IF (ASSOCIATED(this%derived_jacobian)) THEN
         DEALLOCATE(this%derived_jacobian)
         this%derived_jacobian => null()
      END IF

      IF (ASSOCIATED(this%hessian)) THEN
         DEALLOCATE(this%hessian)
         this%hessian => null()
      END IF

      IF (ASSOCIATED(this%gradient)) THEN
         DEALLOCATE(this%gradient)
         this%gradient => null()
      END IF

      IF (ASSOCIATED(this%delta_a)) THEN
         DEALLOCATE(this%delta_a)
         this%delta_a => null()
      END IF

      IF (ASSOCIATED(this%j_svd_w)) THEN
         DEALLOCATE(this%j_svd_w)
         this%j_svd_w => null()
      END IF

      IF (ASSOCIATED(this%j_svd_u)) THEN
         DEALLOCATE(this%j_svd_u)
         this%j_svd_u => null()
      END IF

      IF (ASSOCIATED(this%j_svd_vt)) THEN
         DEALLOCATE(this%j_svd_vt)
         this%j_svd_vt => null()
      END IF

      IF (ASSOCIATED(this%delta_a_len)) THEN
         DEALLOCATE(this%delta_a_len)
         this%delta_a_len => null()
      END IF

      IF (ASSOCIATED(this%exp_g2)) THEN
         DEALLOCATE(this%exp_g2)
         this%exp_g2 => null()
      END IF

      IF (ASSOCIATED(this%step_size)) THEN
         DEALLOCATE(this%step_size)
         this%step_size => null()
      END IF

      IF (ASSOCIATED(this%num_sv)) THEN
         DEALLOCATE(this%num_sv)
         this%num_sv => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get the number of singular values to use.
!>
!>  Number of sigular values to uses is determined by the various cutoffs.
!>
!>  @param[inout] this   A @ref reconstruction_class instance.
!>  @param[in]    num_sv Total number of sigular values.
!>  @returns The number of sigular values to use.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_k_use(this, num_sv)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: reconstruction_get_k_use
      TYPE (reconstruction_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: num_sv

!  local variables
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: dg2exp
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: exp_eff
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: marg_exp_eff
      LOGICAL, DIMENSION(:), ALLOCATABLE         :: l_all
      INTEGER                                    :: k
      REAL (rprec)                               :: denom
      REAL (rprec)                               :: largest_w

!  Start of executable code
      ALLOCATE(dg2exp(0:num_sv))
      ALLOCATE(exp_eff(0:num_sv))
      ALLOCATE(marg_exp_eff(0:num_sv))

!  Esstimate the expected changes in g^2. Equation 22 in Hanson et. al.
!  doi: 10.1088/0029-5515/49/7/075031
      DO k = 0, num_sv
         dg2exp(k) = reconstruction_get_exp_dg2(this,                          &
     &                                          this%delta_a(k,:))

         this%delta_a_len(k) = SQRT(DOT_PRODUCT(this%delta_a(k,:),             &
     &                                          this%delta_a(k,:)))

!  Equation 23 in Hanson et. al. doi: 10.1088/0029-5515/49/7/075031
         exp_eff(k) = ABS(dg2exp(k))/this%delta_a_len(k)
      END DO

!  Although marginal efficiencies are not strictly defined for the Steepest
!  Descent (index 0) and just one singular value cases, for convenience define
!  them here.
      marg_exp_eff(0:1) = exp_eff(0:1)
      DO k = 2, num_sv
         denom = this%delta_a_len(k) - this%delta_a_len(k - 1)
         IF (denom .eq. 0.0) THEN
            IF (this%delta_a_len(k) .eq. 0.0) THEN
               denom = 1.0E-10_rprec
            ELSE
               denom = 1.0E-10_rprec*this%delta_a_len(k)
            END IF
         END IF
         marg_exp_eff(k) = ABS(dg2exp(k)) - ABS(dg2exp(k - 1))/denom
      END DO

!  Check for cutoffs.
      IF (this%j_svd_w(1) .le. 0.0) THEN
         largest_w = 1.0
      ELSE
         largest_w = this%j_svd_w(1)
      END IF

      ALLOCATE(l_all(num_sv))

      l_all = (this%j_svd_w/largest_w .ge. this%cut_svd) .and.                 &
     &        (exp_eff(1:num_sv) .ge. this%cut_eff) .and.                      &
     &        (marg_exp_eff(1:num_sv) .ge. this%cut_marg_eff) .and.            &
     &        (this%delta_a_len(1:num_sv) .ge. this%cut_delta_a) .and.         &
     &        (ABS(dg2exp(1:num_sv)) .ge. this%cut_dg2)

!  Find the largest number of sinular values to use.
      DO k = num_sv, 1, -1
         IF (l_all(k)) THEN
            reconstruction_get_k_use = k
            EXIT
         END IF
      END DO

      DEALLOCATE(dg2exp)
      DEALLOCATE(exp_eff)
      DEALLOCATE(marg_exp_eff)
      DEALLOCATE(l_all)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the expected dg^2 size
!>
!>  Estimates the change in g^2.
!>  Equation 22 in Hanson et. al. doi:10.1088/0029-5515/49/7/075031
!>  Note: there is a type in equation 22. The - sign in the last term should be
!>  a plus sign.
!>
!>  @param[in] this    A @ref reconstruction_class instance.
!>  @param[in] delta_a The normalized parameter step.
!>  @returns Change expected change in g^2.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_exp_dg2(this, delta_a)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_get_exp_dg2
      TYPE (reconstruction_class), INTENT(in) :: this
      REAL (rprec), DIMENSION(:), INTENT(in)  :: delta_a

!  local arguments
      REAL (rprec)                               :: dg2exp_lin
      REAL (rprec)                               :: dg2exp_quad

!  Start of executable code

!  2e * A * da  
      dg2exp_lin =                                                             &
     &   2.0*DOT_PRODUCT(this%e(this%current_step,:),                          &
     &                   MATMUL(this%jacobian, delta_a))

!  da * A^T * A * da = da * alpha * da where alpha is the hessian matrix.
!  Equation 14 in Hanson et. al. doi: 10.1088/0029-5515/49/7/075031
      dg2exp_quad = DOT_PRODUCT(delta_a, MATMUL(this%hessian, delta_a))
      reconstruction_get_exp_dg2 = dg2exp_lin - dg2exp_quad

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the expected g^2 size
!>
!>  Estimates the g^2.
!>  Equation 22 in Hanson et. al. doi:10.1088/0029-5515/49/7/075031
!>
!>  @param[in] this    A @ref reconstruction_class instance.
!>  @param[in] delta_a The normalized parameter step.
!>  @returns Expected g^2.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_exp_g2(this, delta_a)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_get_exp_g2
      TYPE (reconstruction_class), INTENT(in) :: this
      REAL (rprec), DIMENSION(:), INTENT(in)  :: delta_a

!  Start of executable code
      reconstruction_get_exp_g2 = reconstruction_get_g2(this)                  &
     &                          - reconstruction_get_exp_dg2(this,             &
     &                                                       delta_a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the current g^2
!>
!>  Gets the g^2 value for the current reconstruction step.
!>
!>  @param[in] this A @ref reconstruction_class instance.
!>  @returns Current g^2.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_g2(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_get_g2
      TYPE (reconstruction_class), INTENT(in) :: this

!  Start of executable code
      reconstruction_get_g2 =                                                  &
     &   DOT_PRODUCT(this%e(this%current_step,:),                              &
     &               this%e(this%current_step,:))
      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the last g^2
!>
!>  Gets the g^2 value for the last reconstruction step.
!>
!>  @param[in] this A @ref reconstruction_class instance.
!>  @returns The last g^2.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_lastg2(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_get_lastg2
      TYPE (reconstruction_class), INTENT(in) :: this

!  Start of executable code
      reconstruction_get_lastg2 =                                              &
     &   DOT_PRODUCT(this%e(this%current_step - 1,:),                          &
     &               this%e(this%current_step - 1,:))
      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the change in g^2
!>
!>  Gets the change in g^2 from the last reconstruction step to the current
!>  step.
!>
!>  @param[in] this A @ref reconstruction_class instance.
!>  @returns The last g^2.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_get_dg2(this)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_get_dg2
      TYPE (reconstruction_class), INTENT(in) :: this

!  Start of executable code
      reconstruction_get_dg2 = reconstruction_get_lastg2(this)                 &
     &                       - reconstruction_get_g2(this)
      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Evaluate the e vector for a reconstruction step.
!>
!>  Evaluates the e vector after soving the equilibrium.
!>  @see signal::signal_get_e
!>
!>  @param[inout] this     A @ref reconstruction_class instance.
!>  @param[inout] signals  Array of @ref signal objects.
!>  @param[inout] a_model  A @ref model instance.
!>  @param[inout] eq_steps Number of steps for the equilibrium to take.
!>  @param[in]    iou      Input/output file to write log messages to.
!>  @returns True if the equilibrium converged and false if otherwise.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_eval_e(this, signals, a_model, eq_steps,         &
     &                               iou)
      USE signal
      USE model

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: reconstruction_eval_e
      TYPE (reconstruction_class), INTENT(inout) :: this
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      TYPE (model_class), INTENT(inout)          :: a_model
      INTEGER, INTENT(inout)                     :: eq_steps
      INTEGER, INTENT(in)                        :: iou

!  local arguments
      INTEGER                                    :: i

!  Start of executable code
!  Converge the equilibrium. this%e is initialized to zero in the constructor.
      IF (equilibrium_converge(a_model%equilibrium, eq_steps,                  &
     &                         iou)) THEN
         DO i = 1, SIZE(signals)
            this%e(this%current_step,i) =                                      &
     &         signal_get_e(signals(i)%p, a_model, .false.)
         END DO
         !  Write out the step information.
         reconstruction_eval_e = .true.
      ELSE
         reconstruction_eval_e = .false.
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Evaluate the f vector for a reconstruction step.
!>
!>  Evaluates the f vector after soving the equilibrium.
!>  @see param::param_get_value
!>
!>  @param[inout] this           A @ref reconstruction_class instance.
!>  @param[in]    derived_params Array of @ref param objects.
!>  @param[inout] a_model        A @ref model instance.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_eval_f(this, derived_params, a_model)
      USE param
      USE model

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), INTENT(inout)     :: this
      TYPE (param_pointer), DIMENSION(:), INTENT(in) :: derived_params
      TYPE (model_class), INTENT(inout)              :: a_model

!  local arguments
      INTEGER                                        :: i

!  Start of executable code
      DO i = 1, SIZE(derived_params)
         this%f(this%current_step,i) =                                         &
     &      param_get_value(derived_params(i)%p, a_model)
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Evaluate the proper quasi-Newton step.
!>
!>  Evaluates the jacobian, performs the SVD and computes the quasi-Newton step.
!>
!>  @param[inout] this           A @ref reconstruction_class instance.
!>  @param[inout] signals        Array of @ref signal objects.
!>  @param[in]    derived_params Array of derived @ref param objects.
!>  @param[inout] a_model        A @ref model instance.
!>  @param[inout] params         Array of @ref parameter objects.
!>  @param[inout] eq_steps       Number of steps for the equilibrium to take.
!>  @param[in]    iou            Input/output file to write log messages to.
!>  @returns The number of singular values to use.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_eval_step(this, signals, derived_params,         &
     &                                  a_model, params, eq_steps, iou)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: reconstruction_eval_step
      TYPE (reconstruction_class), INTENT(inout)         :: this
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      TYPE (param_pointer), DIMENSION(:), INTENT(in)     ::                    &
     &   derived_params
      TYPE (model_class), INTENT(inout)                  :: a_model
      TYPE (param_pointer), DIMENSION(:), INTENT(inout)  :: params
      INTEGER, INTENT(inout)                             :: eq_steps
      INTEGER, INTENT(in)                                :: iou

!  local arguments
      INTEGER                                    :: i, j
      REAL (rprec)                               :: param_value
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE  :: temp_jacobian
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: temp_work
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: temp_work2
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: temp_work3
      INTEGER                                    :: svd_status

!  Start of executable code
      CALL equilibrium_save_state(a_model%equilibrium)

      DO j = 1, SIZE(params)
!  Save the initial parameter value and the equilirbium state.
         param_value = param_get_value(params(j)%p, a_model)

!  Increment the parameter. param_increment checks the bounds and adjust the
!  step size to keep the parameter in range.
         CALL param_increment(params(j)%p, a_model)

         CALL assert(equilibrium_converge(a_model%equilibrium,                 &
     &                                    eq_steps, iou),                      &
     &               'reconstruction_eval_jac: failed to converge ' //         &
     &               'equilibrium. Try decreasing rp_vrnc')

!  Compute the normalized jacobian. d e_i/d a_j . Equation 12 in Hanson et. al.
!  doi:10.1088/0029-5515/49/7/075031
!
!    a_j = p_j/|pi_j|                                                        (1)
!
!  Where p is the jth parameter, and pi is the jth parameter vrnc. For the
!  finite difference, the partial derivative becomes
!
!    d e_i/d a_j = (e'_i - e_i)/(a'_j - a_j)                                 (2)
!
!  Where the primed components are values for the incremented value.
!
!    a'_j = (p_j + pi_j)/|pi_j|                                              (3)
!
!  Subsituting equation 2 into 3, reduces to
!
!    d e_i/d a_j = (e'_i - e_i)*|pi_j|/pi_j                                  (4)
!
!  In this coding e_i = Sqrt(signal_get_g2)
!  and |pi_j|/pi_j = SIGN(1,delta_param)

         DO i = 1, SIZE(signals)
            this%jacobian(i,j) =                                               &
     &         -(signal_get_e(signals(i)%p, a_model, .false.) -                &
     &           this%e(this%current_step,i))
         END DO
         this%jacobian(:,j) = this%jacobian(:,j)                               &
     &                      * SIGN(1.0_rprec, params(j)%p%delta)

!  Compute the derived parameter jacobian. d f_i/d a_j . This will parallel the
!  signal jacobian with the exception that the parameters aren't normalized.
!  Using equation (3) above the derived jacobian becomes
!
!    d f_i/d a_j = (f'_i - f_i)*|pi_j|/pi_j                                  (5)
!
!  Again in this coding |pi_j|/pi_j = SIGN(1,delta_param)
         DO i = 1, SIZE(derived_params)
            this%derived_jacobian(i,j) =                                       &
     &         -(param_get_value(derived_params(i)%p, a_model) -               &
     &           this%f(this%current_step,i))
         END DO
         this%derived_jacobian(:,j) = this%derived_jacobian(:,j)               &
     &                              * SIGN(1.0_rprec, params(j)%p%delta)

!  Restore the parameter and the equilibirum to it's initial value.
         CALL param_set_value(params(j)%p, a_model, param_value)
         CALL equilibrium_reset_state(a_model%equilibrium)

      END DO

!  Compute the positive definite Hessian, matrix. Equation 14 in Hanson et. al.
!  doi:10.1088/0029-5515/49/7/075031
!
!    alpha = A^T * A                                                         (6)
!
!  Numerical Recipes Third edition Equation 15.5.7 . In here, the second
!  derivatives are assumed to be small and ignored.
      this%hessian = MATMUL(TRANSPOSE(this%jacobian), this%jacobian)

!  Compute the gradient vector in parameter space. Equation 13 in Hanson et. al.
!  doi:10.1088/0029-5515/49/7/075031
!
!    beta = A^T * e                                                          (7)
      this%gradient = MATMUL(TRANSPOSE(this%jacobian),                         &
     &                       this%e(this%current_step,:))

!  Compute the steepest descent step size in normalized parameter space.
!
!    da = beta * beta * beta / (beta * alpha * beta)                         (8)
      this%delta_a(0,:) = this%gradient*DOT_PRODUCT(this%gradient,             &
     &                                              this%gradient)             &
     &                  / DOT_PRODUCT(this%gradient,                           &
     &                                MATMUL(this%hessian,                     &
     &                                       this%gradient))

!  Singular value decomposition of the jacobian. Need to allocate a new jacobian
!  matrix because dgesvd will destroy the one passed to it. dgesvd is from
!  LAPACK
      ALLOCATE(temp_jacobian(SIZE(signals),SIZE(params)))
      temp_jacobian = this%jacobian

      ALLOCATE(temp_work(5*MAX(SIZE(signals), SIZE(params))))
      temp_work = 0.0

      CALL dgesvd('All', 'All', SIZE(signals), SIZE(params),                   &
     &            temp_jacobian, SIZE(signals), this%j_svd_w,                  &
     &            this%j_svd_u, SIZE(signals), this%j_svd_vt,                  &
     &            SIZE(params), temp_work, SIZE(temp_work), svd_status)
      CALL assert_eq(0, svd_status, 'reconstruction_eval_step: ' //            &
     &               'dgesvd problem')

      DEALLOCATE(temp_work)
      DEALLOCATE(temp_jacobian)

!  Reuse the temp work array.
      ALLOCATE(temp_work(SIZE(this%j_svd_w)))
      ALLOCATE(temp_work2(SIZE(signals)))
      ALLOCATE(temp_work3(SIZE(params)))

!  Define the inverse singular values.
      WHERE (this%j_svd_w .gt. 0.0)
         temp_work = 1.0/this%j_svd_w
      ELSE WHERE
         temp_work = 0.0
      END WHERE

!  U^T * e
      temp_work2 = MATMUL(TRANSPOSE(this%j_svd_u),                             &
     &                    this%e(this%current_step,:))

!  Peform pseudo-inversion for successive numbers of singular values retained
      DO i = 1, SIZE(this%j_svd_w)
         temp_work3(1:i) = temp_work(1:i)*temp_work2(1:i)
         temp_work3(i + 1:SIZE(params)) = 0.0
         this%delta_a(i,:) = MATMUL(TRANSPOSE(this%j_svd_vt),                  &
     &                              temp_work3)
      END DO

      DEALLOCATE(temp_work)
      DEALLOCATE(temp_work2)
      DEALLOCATE(temp_work3)

      reconstruction_eval_step =                                               &
     &   reconstruction_get_k_use(this, SIZE(this%j_svd_w))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Take a straight line step.
!>
!>  Determines the normalized change in parameters for a straight line step.
!>
!>  @param[inout] this      A @ref reconstruction_class instance.
!>  @param[in]    k_use     Number of singular values to use.
!>  @param[in]    step_size Current step size limit.
!>  @param[out]   delta_a   Normalized change in parameter.
!>  @returns The step size used.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_sl_step(this, k_use, step_size, delta_a)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_sl_step
      TYPE (reconstruction_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: k_use
      REAL (rprec), INTENT(in)                   :: step_size
      REAL (rprec), DIMENSION(:), INTENT(out)    :: delta_a

!  Start of executable code

      IF ((step_size .ge. 0.0) .and.                                           &
     &    (this%delta_a_len(k_use) .gt. step_size)) THEN
         delta_a = this%step_max/this%delta_a_len(k_use)                       &
     &           * this%delta_a(k_use,:)
         reconstruction_sl_step = SQRT(DOT_PRODUCT(delta_a, delta_a))
      ELSE
         delta_a = this%delta_a(k_use,:)
         reconstruction_sl_step = this%delta_a_len(k_use)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Take a segmented path step.
!>
!>  Determines the normalized change in parameters for a segmented path step.
!>
!>  @param[inout] this      A @ref reconstruction_class instance.
!>  @param[in]    k_use     Number of singular values to use.
!>  @param[in]    step_size Current step size limit.
!>  @param[out]   delta_a   Normalized change in parameter.
!>  @returns The step size used.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_seg_step(this, k_use, step_size, delta_a)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_seg_step
      TYPE (reconstruction_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: k_use
      REAL (rprec), INTENT(in)                   :: step_size
      REAL (rprec), DIMENSION(:), INTENT(out)    :: delta_a

!  local variables
      INTEGER                                    :: k, k_min
      REAL (rprec)                               :: ysq, zsq, ydz
      REAL (rprec)                               :: aa, bb, cc, xx
      REAL (rprec)                               :: discriminant

!  Start of executable code

      IF ((step_size .ge. 0.0) .and.                                           &
     &    (this%delta_a_len(k_use) .gt. step_size)) THEN
!  First, determine k_min - the largest k value where the step is smaller than
!  step_size.
         k_min = 0
         DO k = k_use - 1, 1, -1
            IF (this%delta_a_len(k) .le. step_size) THEN
               k_min = k
               EXIT
            END IF
         END DO

!  Set up the quadratic
         ysq = DOT_PRODUCT(this%delta_a(k_min,:), this%delta_a(k_min,:))
         zsq = DOT_PRODUCT(this%delta_a(k_min + 1,:),                          &
     &                     this%delta_a(k_min + 1,:))
         ydz = DOT_PRODUCT(this%delta_a(k_min,:),                              &
     &                     this%delta_a(k_min + 1,:))

!  Special case when (k_min .eq. 0) - reconstruction arrays contain SD step, but
!  use the zero step for (kmin .eq. 0)
         IF (k_min .eq. 0) THEN
            ysq = 0.0
            ydz = 0.0
         END IF

         aa = ysq + zsq - 2.0*ydz
         bb = -2.0*ysq + 2.0*ydz
         cc = ysq - step_size**2.0

         discriminant = bb**2.0 - 4.0*aa*cc

         IF (discriminant .le. 0.0) THEN
            xx = 0.0
         ELSE
            xx = (-bb + SQRT(discriminant))/(2.0*aa)
         END IF

         IF ((xx .lt. 0.0) .or. (xx .gt. 1.0)) THEN
            xx = 0.0
         END IF

         delta_a = xx*this%delta_a(k_min + 1,:)
         IF (k_min .ne. 0) THEN
            delta_a = this%delta_a(k_min,:)*(1.0 - xx) + delta_a
         END IF

         reconstruction_seg_step = SQRT(DOT_PRODUCT(delta_a, delta_a))
      ELSE
         delta_a = this%delta_a(k_use,:)
         reconstruction_seg_step = this%delta_a_len(k_use)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Take a Levenberg-Marquardt step.
!>
!>  Determines the normalized change in parameters for a Levenberg-Marquardt
!>  step.
!>
!>  @param[inout] this      A @ref reconstruction_class instance.
!>  @param[in]    k_use     Number of singular values to use.
!>  @param[in]    step_size Current step size limit.
!>  @param[out]   delta_a   Normalized change in parameter.
!>  @returns The step size used.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_lm_step(this, k_use, step_size, delta_a)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_lm_step
      TYPE (reconstruction_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: k_use
      REAL (rprec), INTENT(in)                   :: step_size
      REAL (rprec), DIMENSION(:), INTENT(out)    :: delta_a

!  local variables
      REAL (rprec)                               :: lambda
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: utdote

!  Start of executable code

      IF ((step_size .ge. 0.0) .and.                                           &
     &    (this%delta_a_len(k_use) .gt. step_size)) THEN

!  Define a default value of lambda in case the root find does not work well.
         lambda =                                                              &
     &      (this%j_svd_w(1)*this%delta_a_len(k_use)/step_size)**2.0

         ALLOCATE(utdote(SIZE(this%e, 2)))
         utdote = MATMUL(this%e(this%current_step,:), this%j_svd_u)

!  Find the L-M parameter lambda that corresponds to a step length of step_size.
         lambda = reconstruction_lm_rootfind(this, k_use, lambda,              &
     &                                       step_size, utdote)

!  Find the step
         utdote(1:k_use) =                                                     &
     &      utdote(1:k_use)*this%j_svd_w/(this%j_svd_w**2.0 + lambda)

         delta_a = MATMUL(utdote(1:k_use), this%j_svd_vt(1:k_use,:))
         reconstruction_lm_step = SQRT(DOT_PRODUCT(delta_a, delta_a))

         DEALLOCATE(utdote)
      ELSE
         delta_a = this%delta_a(k_use,:)
         reconstruction_lm_step = this%delta_a_len(k_use)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Root find for Levenberg-Marquardt step.
!>
!>  Find the root of the @ref reconstruction_lm_function for a
!>  Levenberg-Marquardt step.
!>
!>  @param[inout] this           A @ref reconstruction_class instance.
!>  @param[in]    k_use          Number of singular values to use.
!>  @param[in]    lambda_default The default lambda.
!>  @param[in]    step_size      Current step size limit.
!>  @param[in]    utdote         Vector of U^T * e
!>  @returns The root lambda.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_lm_rootfind(this, k_use, lambda_default,         &
     &                                    step_size, utdote)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_lm_rootfind
      TYPE (reconstruction_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                        :: k_use
      REAL (rprec), INTENT(in)                   :: lambda_default
      REAL (rprec), INTENT(in)                   :: step_size
      REAL (rprec), DIMENSION(:), INTENT(in)     :: utdote

!  local variables
      INTEGER                                    :: i
      REAL (rprec), DIMENSION(:), ALLOCATABLE    :: f_sqrd
      REAL (rprec)                               :: step_size_sqrd
      REAL (rprec)                               :: f_0, f_2
      REAL (rprec)                               :: lambda_0, lambda_2
      REAL (rprec)                               :: lambda_try
      REAL (rprec)                               :: lambda_new
      REAL (rprec)                               :: f_try, f_new
      REAL (rprec)                               :: s

!  local parameters
      REAL (rprec), PARAMETER :: d_lambda_min = 1.0E-6_rprec

!  Start of executable code

!  Set up stuff for the function reconstruction_lm_function for root find.
      ALLOCATE(f_sqrd(k_use))

!  Put (U^t*e)^2 in the function.
      f_sqrd = utdote(1:k_use)**2.0
      step_size_sqrd = step_size**2.0

      reconstruction_lm_rootfind = lambda_default

!  Find the bracketing values of lambda. Look for a small value of lambda, to
!  give a positive function value. f_0 should be greater than zero, as the step
!  size at k_use is larger than step_size.
      f_0 = reconstruction_lm_function(this, f_sqrd, step_size_sqrd,           &
     &                                 0.0_rprec)
      IF (f_0 .lt. 0.0) THEN
         CALL err_warn('reconstruction_lm_rootfind: f_0 < 0.0')
         DEALLOCATE(f_sqrd)
         RETURN
      END IF
      lambda_0 = 0.0

!  Look for a large value of lambda, to give a negative function value.
      lambda_2 = this%j_svd_w(1)**2.0
      f_2 = reconstruction_lm_function(this, f_sqrd, step_size_sqrd,           &
     &                                 lambda_2)
      DO i = 1, 20
         IF (f_2 .le. 0.0) THEN
            EXIT
         ELSE
            lambda_2 = 4.0*lambda_2
            f_2 = reconstruction_lm_function(this, f_sqrd,                     &
     &                                       step_size_sqrd,                   &
     &                                       lambda_2)
         END IF
      END DO
      IF (f_2 .gt. 0.0) THEN
         CALL err_warn('reconstruction_lm_rootfind: f_2 > 0.0')
         DEALLOCATE(f_sqrd)
         RETURN
      END IF

!  Check that f_0 and f_2 actually bracket.
      IF (f_0*f_2 .gt. 0.0) THEN
         CALL err_warn('reconstruction_lm_rootfind: fail to bracket')
         DEALLOCATE(f_sqrd)
         RETURN
      END IF

!  Find the root using the Ridder method. Numerical Recipes.
      reconstruction_lm_rootfind = -9.99E99_rprec
      DO i = 1, 99
         lambda_try = 0.5*(lambda_0 + lambda_2)
         f_try = reconstruction_lm_function(this, f_sqrd,                      &
     &                                      step_size_sqrd,                    &
     &                                      lambda_try)
         s = SQRT(f_try**2.0 - f_0*f_2)
         IF (s .eq. 0.0) THEN
            EXIT
         END IF

         lambda_new = lambda_try + (lambda_try - lambda_0)*                    &
     &                             (SIGN(1.0_rprec, f_0 - f_2)*f_try/s)

         IF (ABS(lambda_new - reconstruction_lm_rootfind) .le.                 &
     &       d_lambda_min) THEN
            EXIT
         END IF

         reconstruction_lm_rootfind = lambda_new
         f_new = reconstruction_lm_function(this, f_sqrd,                      &
     &                                      step_size_sqrd,                    &
     &                                      lambda_new)
         IF (f_new .eq. 0.0) THEN
            EXIT
         ELSE IF (SIGN(f_try, f_new) .ne. f_try) THEN
            lambda_0 = lambda_try
            f_0 = f_try
            lambda_2 = lambda_new
            f_2 = f_new
         ELSE IF (SIGN(f_0, f_new) .ne. f_0) THEN
            lambda_2 = lambda_new
            f_2 = f_new
         ELSE IF (SIGN(f_2, f_new) .ne. f_2) THEN
            lambda_0 = lambda_new
            f_0 = f_new
         ELSE
            STOP 'reconstruction_lm_rootfind: should never get here'
         END IF

         IF (ABS(lambda_0 - lambda_2) .le. d_lambda_min) THEN
            EXIT
         END IF
      END DO

      DEALLOCATE(f_sqrd)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Levenberg-Marquardt function.
!>
!>  Function used by @ref reconstruction_lm_rootfind for a
!>  Levenberg-Marquardt step.
!>
!>  @param[inout] this           A @ref reconstruction_class instance.
!>  @param[in]    f_sqrd         f^2
!>  @param[in]    step_size_sqrd Square of the current step size.
!>  @param[in]    lambda         Lambda value.
!>  @returns F(lambda)
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_lm_function(this, f_sqrd, step_size_sqrd,        &
     &                                    lambda)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: reconstruction_lm_function
      TYPE (reconstruction_class), INTENT(inout) :: this
      REAL (rprec), DIMENSION(:), INTENT(in)     :: f_sqrd
      REAL (rprec), INTENT(in)                   :: step_size_sqrd
      REAL (rprec), INTENT(in)                   :: lambda

!  local variables
      INTEGER                                    :: i

!  Start of executable code
      reconstruction_lm_function = 0.0
      DO i = 1, SIZE(f_sqrd)
         reconstruction_lm_function = reconstruction_lm_function               &
     &      + f_sqrd(i)*(this%j_svd_w(i)/(this%j_svd_w(i)**2.0 +               &
     &                                    lambda))**2.0
      END DO

      reconstruction_lm_function = reconstruction_lm_function                  &
     &                           - step_size_sqrd

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Step the reconstruction.
!>
!>  Advance the reconstruction by a single step. This is the main V3FIT
!>  algorithm.
!>
!>  @param[inout] this     A @ref reconstruction_class instance.
!>  @param[inout] signals  Array of @ref signal objects.
!>  @param[inout] a_model  A @ref model instance.
!>  @param[inout] params   Array of @ref parameter objects.
!>  @param[inout] eq_steps Number of steps for the equilibrium to take.
!>  @param[in]    iou      Input/output file to write log messages to.
!>  @returns True if the step was successfull and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION reconstruction_step(this, signals, derived_params,              &
     &                             a_model, params, eq_steps, iou)
      USE signal
      USE model
      USE param

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL :: reconstruction_step
      TYPE (reconstruction_class), INTENT(inout)         :: this
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      TYPE (param_pointer), DIMENSION(:), INTENT(in)     ::                    &
     &   derived_params
      TYPE (model_class), INTENT(inout)                  :: a_model
      TYPE (param_pointer), DIMENSION(:), INTENT(inout)  :: params
      INTEGER, INTENT(inout)                             :: eq_steps
      INTEGER, INTENT(in)                                :: iou

!  local variables
      INTEGER                                            :: k_use
      INTEGER                                            :: i, j
      REAL (rprec), DIMENSION(:), ALLOCATABLE            :: param_value
      REAL (rprec), DIMENSION(:), ALLOCATABLE            :: delta_a
      REAL (rprec)                                       :: step_use

!  Start of executable code
      reconstruction_step = .false.

!  Evaluate the quasi-Newton step.
      k_use = reconstruction_eval_step(this, signals, derived_params,          &
     &                                 a_model, params, eq_steps, iou)

      ALLOCATE(param_value(SIZE(params)))
      ALLOCATE(delta_a(SIZE(params)))

!  Save the parameter values.
      DO i = 1, SIZE(params)
         param_value(i) = param_get_value(params(i)%p, a_model)
      END DO

!  Set the initial step size.
      step_use = this%step_max

      DO i = 1, reconstruction_max_step_try
!  Step the parameters using the specified step method.
         SELECT CASE (this%step_type)

            CASE (reconstruction_sl_step_type)
               step_use = reconstruction_sl_step(this, k_use, step_use,        &
     &                                           delta_a)
               this%exp_g2(this%current_step + 1) =                            &
     &            reconstruction_get_exp_g2(this, delta_a)
               this%step_size(this%current_step + 1) = step_use
               this%num_sv(this%current_step + 1) = k_use
               CALL reconstruction_write_step(this, 'sl', iou)

            CASE (reconstruction_lm_step_type)
               step_use = reconstruction_lm_step(this, k_use, step_use,        &
     &                                           delta_a)
               this%exp_g2(this%current_step + 1) =
     &            reconstruction_get_exp_g2(this, delta_a)
               this%step_size(this%current_step + 1) = step_use
               this%num_sv(this%current_step + 1) = k_use
               CALL reconstruction_write_step(this, 'lm', iou)

            CASE (reconstruction_seg_step_type)
               step_use = reconstruction_seg_step(this, k_use, step_use,       &
     &                                            delta_a)
               this%exp_g2(this%current_step + 1) =
     &            reconstruction_get_exp_g2(this, delta_a)
               this%step_size(this%current_step + 1) = step_use
               this%num_sv(this%current_step + 1) = k_use
               CALL reconstruction_write_step(this, 'seg', iou)

         END SELECT

!  Set the new parameter values. ABS(params(j)%p%delta) is the normalization
!  factor. This is not vrnc because the increment function can take a non vrnc
!  step size.
         DO j = 1, SIZE(params)
            CALL param_set_value(params(j)%p, a_model, param_value(j) +        &
     &                           delta_a(j)*ABS(params(j)%p%delta))
         END DO

!  Reconverge the equilibrium. Increment the current step so
!  reconstruction_eval_e writes the e array to the correct index. If the
!  reconstruction step should fail for any reason. decrement the step back down.
         this%current_step = this%current_step + 1
         IF (reconstruction_eval_e(this, signals, a_model,                     &
     &                             eq_steps, iou)) THEN
!  Equilibrium converged
            CALL reconstruction_eval_f(this, derived_params, a_model)

            IF (reconstruction_get_dg2(this) .ge. 0.0) THEN
               reconstruction_step = .true.
               WRITE (*,*)
               WRITE (*,*) ' *** Step succeeded'
               WRITE (*,*) ' Result of step'
               WRITE (*,1000)
               WRITE (*,1001) reconstruction_get_lastg2(this),                 &
     &                        reconstruction_get_g2(this),                     &
     &                        reconstruction_get_dg2(this)

               WRITE (iou,*)
               WRITE (iou,*) ' *** Step succeeded'
               WRITE (iou,*) ' Result of step'
               WRITE (iou,1000)
               WRITE (iou,1001) reconstruction_get_lastg2(this),               &
     &                          reconstruction_get_g2(this),                   &
     &                          reconstruction_get_dg2(this)
               EXIT
            ELSE
!  Equilibrium converged but g^2 increased
               WRITE (*,*) 'Equilibrium converged but g2 increased'
               WRITE (iou,*) 'Equilibrium converged but g2 increased'
            END IF
         ELSE
!  Equilibrium did not converge
            WRITE (*,*) 'Equilibrium did not converge'
            WRITE (iou,*) 'Equilibrium did not converge'
         END IF

!  If this section was reached, the reconstruction step failed. Reset the
!  current step to the last good step.
         this%current_step = this%current_step - 1
         IF (i .eq. reconstruction_max_step_try) THEN
            WRITE (*,*) 'Reconstruction step failed.'
            WRITE (iou,*) 'Reconstruction step failed.'
            EXIT
         END IF

!  Retry step with a smaller stepsize.
         step_use = step_use/2.0
         WRITE (*,*) 'Changing max step size to ', step_use
         WRITE (iou,*) 'Changing max step size to ', step_use

!  Reset the equilibrium state back to the last step. The state of the last step
!  was saved before the jacobian was computed.
         CALL equilibrium_reset_state(a_model%equilibrium)
      END DO

      DEALLOCATE(param_value)
      DEALLOCATE(delta_a)

1000  FORMAT('g^2 old',7x,'g^2 new',7x,'g^2 diff')
1001  FORMAT(es12.5,2x,es12.5,2x,es12.5)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Evaluate the parameter covariance and signal effectiveness.
!>
!>  Calculates the signal effectiveness matrix. This value is stored in
!>  @ref parameter::param_class::sem. To fine the signal effectiveness, an
!>  intermediate step compute the parameter covariance.
!>
!>  @param[inout] this           A @ref reconstruction_class instance.
!>  @param[inout] params         Array of @ref parameter objects.
!>  @param[inout] signals        Array of @ref signal objects.
!>  @param[inout] derived_params Array of @ref parameter objects.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_eval_sem(this, params, signals,                &
     &                                   derived_params)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), INTENT(inout)         :: this
      TYPE (param_pointer), DIMENSION(:), INTENT(inout)  :: params
      TYPE (signal_pointer), DIMENSION(:), INTENT(inout) :: signals
      TYPE (param_pointer), DIMENSION(:), INTENT(inout)  ::                    &
     &   derived_params

!  local arguments
      INTEGER                                            :: j
      REAL (rprec), DIMENSION(:,:), POINTER              :: cp, cd
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE          :: j_dot_cp

!  local parameters
      CHARACTER (len=*), PARAMETER ::                                          &
     &   sub_name = 'reconstruction_eval_sem'

!  Start of executable code
!  Since the jacobian that will be used is normalized, the signal covaiance
!  must be normalized as well. This results in an idenity matrix of
!  SIZE(signals) X SIZE(signals).
!
!    C_p^-1 = J^T * C_s^-1 * J -> A^T * I * A = A^T * A                      (1)
!
!  Equation 3 in Hanson et. al. doi:10.1088/0029-5515/49/7/075031
!  Note: This is the same as the hessian that was computed for the
!  reconstruction. Use that value to avoid a matrix multiplication.
      ALLOCATE(cp(SIZE(params),SIZE(params)))
      cp = this%hessian

!  Invert C_p^-1 to get the parameter covariance matrix.
      CALL reconstruction_invert_matrix(cp, sub_name, 'C_p^-1')

!  From the parameter covariance matrix, map the reconstructed uncertainty to
!  the derived parameters by the derived jacobian.
!
!    C_d = J' * C_p * J'^T                                                   (2)
      IF (SIZE(derived_params) .gt. 0) THEN
         ALLOCATE(cd(SIZE(derived_params), SIZE(derived_params)))
         cd = MATMUL(this%derived_jacobian,                                    &
     &               MATMUL(cp, TRANSPOSE(this%derived_jacobian)))

         DO j = 1, SIZE(derived_params)
            derived_params(j)%p%sigma = SQRT(cd(j,j))
         END DO

         DEALLOCATE(cd)
      END IF

!  Store the parameter covariance. Since the normalized signal covariance matrix
!  was used, the parameter covariance needs to be denormalized.
      DO j = 1, SIZE(params)
         params(j)%p%sigma = SQRT(cp(j,j)*params(j)%p%delta**2.0)
      END DO

!  Compute signal effectiveness.
      ALLOCATE(j_dot_cp(SIZE(signals),SIZE(params)))
      j_dot_cp = MATMUL(this%jacobian, cp)

!  The signal effectiveness for each parameter is stored in the parameter
!  itself. Equation 6 in Hanson et. al. doi:10.1002/ctpp.200900059
!
!    (J * cp)^2 / (sigma_p^2 * sigma_s^2)                                    (3)
!
!  Due to use of the normalized jacobian and normalized parameter covariance
!  matrix, the 1/sigma_s^2 term is found in the J^2 term. An extra Pi term must
!  be multiplied out. In terms of the normalized quantities, becomes
!
!    (A * cp_n)^2 * pi^2 / sigma_p^2                                         (4)
      DO j = 1, SIZE(params)
         params(j)%p%sem = (j_dot_cp(:,j)**2.0)                                &
     &                   * params(j)%p%delta**2.0                              &
     &                   / params(j)%p%sigma**2.0
      END DO

      DEALLOCATE(j_dot_cp, cp)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Invert an M x N matrix.
!>
!>  Perform a pseudo-inversion using a singular value decomposition.
!>
!>  @param[inout] matrix The matrix to invert.
!>  @param[in]    sub    The subroutine of the matrix this was called from.
!>  @param[in]    name   The name of the matrix to be inverted.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_invert_matrix(matrix, sub, name)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(:,:) ,INTENT(inout), POINTER :: matrix
      CHARACTER (len=*), INTENT(in)               :: sub
      CHARACTER (len=*), INTENT(in)               :: name

!  local arguments
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE   :: svd_vt
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE   :: svd_u
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: svd_w
      REAL (rprec), DIMENSION(:), ALLOCATABLE     :: svd_work
      INTEGER                                     :: num_svd
      INTEGER                                     :: work_size
      INTEGER                                     :: i, status
      INTEGER                                     :: m, n

!  Start of executable code
      m = SIZE(matrix, 1)
      n = SIZE(matrix, 2)

      num_svd = MIN(m, n)

      IF (num_svd .gt. 0) THEN

         ALLOCATE(svd_vt(n,n))
         ALLOCATE(svd_u(m,m))
         ALLOCATE(svd_w(num_svd))
         ALLOCATE(svd_work(1))

         svd_vt = 0.0
         svd_u = 0.0
         svd_w = 0.0
         svd_work = 0.0

!  Find the optimal work size.
         CALL dgesvd('All', 'All', m, n, matrix, m, svd_w, svd_u,              &
     &               m, svd_vt, n, svd_work, -1, status)
         work_size = INT(svd_work(1))
         DEALLOCATE(svd_work)
         ALLOCATE(svd_work(work_size))

!  Factor the matrix to M = U * W * V^T
         CALL dgesvd('All', 'All', m, n, matrix, m, svd_w, svd_u,              &
     &               m, svd_vt, n, svd_work, work_size, status)
         CALL assert_eq(0, status, sub // ': dgesvd problem when ' //          &
     &                  'inverting ' // name)

!  If the matrix is rectangular the inverse matrix will be transposed.
         DEALLOCATE(matrix)
         ALLOCATE(matrix(n,m))

!  Invert the nonzero elements of the W matrix.
         matrix = 0
         DO i = 1, num_svd
            IF (svd_w(i) .gt. 0.0) THEN
               matrix(i,i) = 1.0/svd_w(i)
            END IF
         END DO

!  Invert the matric M^-1 = V * W^-1 * U^T
         matrix = MATMUL(TRANSPOSE(svd_vt),                                    &
     &                   MATMUL(matrix, TRANSPOSE(svd_u)))
         DEALLOCATE(svd_vt, svd_u, svd_w, svd_work)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out a reconstruction.
!>
!>  Write out the final information about a reconstruction. This writes out a
!>  table of the g^2, expected g^2, number of sigular values used and the
!>  normalized step size used.
!>
!>  @param[in] this A @ref reconstruction_class instance.
!>  @param[in] iou  Input/output file to write final result to.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_write(this, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), POINTER, INTENT(in) :: this
      INTEGER, INTENT(in)                              :: iou

!  local arguments
      INTEGER                                          :: i
      REAL (rprec)                                     :: temp_g2

!  Start of executable code
      WRITE (*,*)
      WRITE (iou,*)
      WRITE (*,*) ' *** Reconstruction Steps'
      WRITE (iou,*) ' *** Reconstruction Steps'
      WRITE (*,1000)
      WRITE (iou,1000)

      temp_g2 = DOT_PRODUCT(this%e(0,:), this%e(0,:))
      WRITE (*,1001) 0, temp_g2
      WRITE (iou,1001) 0, temp_g2

      DO i = 1, this%current_step
         temp_g2 = DOT_PRODUCT(this%e(i,:), this%e(i,:))
         WRITE (*,1001) i, temp_g2, this%exp_g2(i), this%num_sv(i),            &
     &                  this%step_size(i)
         WRITE (iou,1001) i, temp_g2, this%exp_g2(i), this%num_sv(i),          &
     &                    this%step_size(i)
      END DO

1000  FORMAT('irstep',2x,'g^2',11x,'Expected g^2',4x,'sv ',2x,                 &
     &       'normalized step size')
1001  FORMAT(i6,2(2x,es12.5),2x,i3,2x,es12.5)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out a reconstruction step.
!>
!>  Write out information about a reconstruction step. This writes the number of
!>  singular values use, the step type, the normalized step size and the
!>  expected g^2.
!>
!>  @param[in] this      A @ref reconstruction_class instance.
!>  @param[in] step_type The step type used.
!>  @param[in] iou       Input/output file to write final result to.
!-------------------------------------------------------------------------------
      SUBROUTINE reconstruction_write_step(this, step_type, iou)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (reconstruction_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)           :: step_type
      INTEGER, INTENT(in)                     :: iou

!  local variables.
      CHARACTER (len=21)                      :: jac_format
      INTEGER                                 :: i

!  Start of executable code
      WRITE (*,*)
      WRITE (*,1000) this%num_sv(this%current_step + 1)
      WRITE (*,1001) step_type, this%step_size(this%current_step + 1)
      WRITE (*,1002) this%exp_g2(this%current_step + 1)

      WRITE (iou,*)
      WRITE (iou,1000) this%num_sv(this%current_step + 1)
      WRITE (iou,1001) step_type, this%step_size(this%current_step + 1)
      WRITE (iou,1002) this%exp_g2(this%current_step + 1)

      WRITE (iou,*)
      WRITE (iou,*) ' *** Jacobian'
      WRITE (jac_format,1003) SIZE(this%jacobian, 2)
      DO i = 1, SIZE(this%jacobian, 1)
         WRITE(iou,jac_format) i, this%jacobian(i,:)
      END DO
      WRITE (iou,*)

1000  FORMAT('Number of singular values used: 'i3)
1001  FORMAT(a4,': Step size limited to ',es12.5)
1002  FORMAT('  Expected g^2 = ',es12.5)
1003  FORMAT('(i4,',i3,'(2x,es12.5))')

      END SUBROUTINE

      END MODULE
