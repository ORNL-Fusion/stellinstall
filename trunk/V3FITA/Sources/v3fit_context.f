!*******************************************************************************
!  File v3fit_context.f
!  Contains module v3fit_context
!  Defines the base class of the type v3fit_context_class:
!
!  Super Class: None
!  Sub Classes: None
!
!*******************************************************************************

      MODULE v3fit_context
      USE v3fit_input
      USE param
      USE commandline_parser
      USE reconstruction
      USE ezcdf

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3fit context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  This object holds all the memory for running v3fit.
!-------------------------------------------------------------------------------
      TYPE v3fit_context_class
         REAL (rprec)                                 :: recon_stop

         INTEGER                                      :: runlog_iou
         INTEGER                                      :: recout_iou

         TYPE (commandline_parser_class), POINTER :: cl_parser => null()
         TYPE (model_class), POINTER                  :: model => null()
         TYPE (signal_pointer), DIMENSION(:), POINTER ::                       &
     &      signals => null()
         TYPE (param_pointer), DIMENSION(:), POINTER :: params => null()
         TYPE (reconstruction_class), POINTER         :: recon => null()

!  The default for the following index is -1 indicating that this signal hasn't
!  been created.
         INTEGER :: magnetic_index = -1
         INTEGER :: sxrem_index = -1
         INTEGER :: intpol_index = -1
         INTEGER :: thomson_index = -1
         INTEGER :: extcurz_index = -1
         INTEGER :: limiter_index = -1
         INTEGER :: combination_index = -1

!  Result file NetCDF variables
         INTEGER :: result_ncid
         INTEGER :: nsteps_id
         INTEGER :: g2_id

         INTEGER :: param_value_id
         INTEGER :: param_sigma_id

         INTEGER :: signal_observed_value_id
         INTEGER :: signal_model_value_id
         INTEGER :: signal_sigma_value_id
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
      FUNCTION v3fit_context_construct(cl_parser)
      USE safe_open_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER :: v3fit_context_construct
      TYPE (commandline_parser_class), POINTER :: cl_parser

!  local variables
      CHARACTER (len=data_filename_length)     :: filename
      INTEGER                                  :: status

!  Start of executable code
      filename = commandline_parser_get_string(cl_parser, '-file')

      ALLOCATE(v3fit_context_construct)

!  Setup the runlog file.
      CALL safe_open(v3fit_context_construct%runlog_iou, status,               &
     &               TRIM('runlog.' // filename_base(filename)),               &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0, status, 'v3fit_context_construct: ' //                 &
     &               'Safe_open of runlog. ' //                                &
     &               TRIM(filename_base(filename)) // 'failed')

!  Setup the recout file.
      CALL safe_open(v3fit_context_construct%recout_iou, status,               &
     &               TRIM('recout.' // filename_base(filename)),               &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0, status, 'v3fit_context_construct: ' //                 &
     &               'Safe_open of recout. ' //                                &
     &               TRIM(filename_base(filename)) // 'failed')

!  Create a result file.
      status = nf_create('runlog.' // TRIM(filename_base(filename)) //         &
     &                   '.nc', nf_clobber,                                    &
     &                   v3fit_context_construct%result_ncid)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      WRITE (v3fit_context_construct%runlog_iou, *) 'V3FITA RUN'
      WRITE (v3fit_context_construct%runlog_iou, *)                            &
     &   '  Namelist Input from file ', TRIM(filename)

      WRITE (v3fit_context_construct%recout_iou, *) 'V3FITA RUN'
      WRITE (v3fit_context_construct%recout_iou, *)                            &
     &   '  Namelist Input from file ', TRIM(filename)

      v3fit_context_construct%cl_parser => cl_parser

!  Allocate the number of of signals needed. This is allocated to the size
!  na_s_desc. na_s_desc will default to the maximum number of signals needed.
      ALLOCATE(v3fit_context_construct%signals(na_s_desc))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
      SUBROUTINE v3fit_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER, INTENT(inout) :: this

!  local variables
      INTEGER                                            :: i
      INTEGER                                            :: status

!  Start of executable code
      CLOSE(this%runlog_iou)
      CLOSE(this%recout_iou)

      status = nf_close(this%result_ncid)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%cl_parser)) THEN
         CALL commandline_parser_destruct(this%cl_parser)
         this%cl_parser => null()
      END IF

      IF (ASSOCIATED(this%model)) THEN
         CALL model_destruct(this%model)
         this%model => null()
      END IF

!  Deconstruct and deallocate all the signals.
      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            IF (ASSOCIATED(this%signals(i)%p)) THEN
               CALL signal_destruct(this%signals(i)%p)
               this%signals(i)%p => null()
            END IF
         END DO
         DEALLOCATE(this%signals)
         this%signals => null()
      END IF

!  Deconstruct and deallocate all the parameters.
      IF (ASSOCIATED(this%params)) THEN
         DO i = 1, SIZE(this%params)
            IF (ASSOCIATED(this%params(i)%p)) THEN
               CALL param_destruct(this%params(i)%p)
               this%params(i)%p => null()
            END IF
         END DO
        DEALLOCATE(this%params)
        this%params => null()
      END IF

      IF (ASSOCIATED(this%recon)) THEN
         CALL reconstruction_destruct(this%recon)
         this%recon => null()
      END IF

!  Reset teh signal index
      this%magnetic_index = -1
      this%sxrem_index = -1
      this%intpol_index = -1
      this%thomson_index = -1
      this%extcurz_index = -1
      this%limiter_index = -1
      this%combination_index = -1

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Reduce memory useage by resizing the arrays to the smallest necessary size.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_resize(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout)    :: this

!  local variables
      INTEGER                                      :: i
      INTEGER                                      :: minsize = 0
      TYPE (signal_pointer), DIMENSION(:), POINTER :: temp_signal

!  Start of executable code
      WRITE (*,*) ' *** Resizing context arrays'
      WRITE (this%runlog_iou,*) ' *** Resizing context arrays'

!  Search an array for the first instance of a null pointer.
      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            IF (.not.ASSOCIATED(this%signals(i)%p)) THEN
               minsize = i - 1
               EXIT
            END IF
         END DO

!  Reallocate the array if it is smaller
         IF (minsize .lt. SIZE(this%signals)) THEN
            ALLOCATE(temp_signal(minsize))
            DO i = 1, minsize
                temp_signal(i)%p => this%signals(i)%p
            END DO

            DEALLOCATE(this%signals)
            this%signals => temp_signal
         END IF

      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Write out the context
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_write(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout)    :: this

!  local variables
      INTEGER                                      :: i

!  Start of executable code
      WRITE (*,*) ' *** Writting context to disk'
      WRITE (this%runlog_iou,*) ' *** Writting context to disk'

!  Write out the reconstruction
      IF (ASSOCIATED(this%recon)) THEN
         CALL reconstruction_write(this%recon, this%recout_iou)
      END IF

!  Write out the parameters
      IF (ASSOCIATED(this%params)) THEN
         CALL param_write_header(this%recout_iou)
         DO i = 1, SIZE(this%params)
            CALL param_write(this%params(i)%p, this%recout_iou, i,             &
     &                       this%model)
         END DO
      END IF

!  Write out the model.
      IF (ASSOCIATED(this%model)) THEN
         CALL model_write(this%model, this%recout_iou)
      END IF

!  Write out the signals
      DO i = 1, SIZE(this%signals)
!  Check if the index has reached a new signal type
         IF ((this%magnetic_index .eq. i) .or.                                 &
     &       (this%sxrem_index .eq. i) .or.                                    &
     &       (this%intpol_index .eq. i) .or.                                   &
     &       (this%thomson_index .eq. i) .or.                                  &
     &       (this%extcurz_index .eq. i) .or.                                  &
     &       (this%limiter_index .eq. i) .or.                                  &
     &       (this%combination_index .eq. i)) THEN

            CALL signal_write_header(this%signals(i)%p,                        &
     &                               this%recout_iou)
         END IF

         CALL signal_write(this%signals(i)%p, this%recout_iou, i,              &
     &                     this%model)
      END DO

!  All the signals should be written out before auxiliary information.
      DO i = 1, SIZE(this%signals)
         CALL signals_write_auxiliary(this%signals(i)%p,                      &
     &                                this%recout_iou, i)
      END DO

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!
!  Dimensions:
!    maxnsteps    The max number of reconstruction steps.
!    nparam       Number of reconstruction parameters.
!    nparamindex  Number of parameter indices.
!    nsignal      Number of signals.
!    string_len   Max length of strings.
!
!    *** Equilibrium Dimesnions ***
!    See equilibrium objects.
!
!  Multi dimensional arrays need to be transposed so arrays appear in the
!  correct order in non fortran languages.
!  Variables:
!    nsteps                          Actual number of steps taken.
!    g2 (maxnsteps)                  g^2 value for each reconstruction step.
!
!    param_name (string_len,         Name of the reconstruction parameter.
!                nparam)
!    param_index (nparamindex,       Indices of the reconstruction parameter.
!                 nparam)      
!    param_value (nparam,            Parameter values for each reconstruction step.
!                 maxnsteps)
!    param_sigma (nparam,            Parameter sigmas for each reconstruction step.
!                 maxnsteps)
!
!    signal_name (string_len,        s_name of the signal
!                 nsignal)
!    signal_type (string_len,        Description of the signal.
!                 nsignal)
!    signal_weight (nsignal)         Weight value of the signal.
!    signal_observed_value (nsignal, Observed value of the signal at each
!                        maxnsteps)  reconstruction step.
!    signal_model_value (nsignal,    Modeled value of the signal at each
!                        maxnsteps)  reconstruction step.
!    signal_sigma_value (nsignal,    Sigma value of the signal at each
!                        maxnsteps)  reconstruction step.
!
!    *** Equilibrium Variables ***
!    See equilibrium objects.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Define and write inital data to the defined variables.
!-------------------------------------------------------------------------------
      SUBROUTINE constext_init_data(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: status
      INTEGER                                   :: i
      INTEGER                                   :: maxnsteps_dim_id
      INTEGER                                   :: nparam_dim_id
      INTEGER                                   :: nparamindex_dim_id
      INTEGER                                   :: nsignal_dim_id
      INTEGER                                   :: string_len_dim_id
      INTEGER                                   :: param_name_id
      INTEGER                                   :: param_index_id
      INTEGER                                   :: signal_name_id
      INTEGER                                   :: signal_type_id
      INTEGER                                   :: signal_weight_id

!  Start of executable code
!  Define dimensions
      IF (ASSOCIATED(this%recon)) THEN
         status = nf_def_dim(this%result_ncid, 'maxnsteps',                    &
     &                       SIZE(this%recon%e), maxnsteps_dim_id)
      ELSE
         status = nf_def_dim(this%result_ncid, 'maxnsteps',                    &
     &                       1, maxnsteps_dim_id)
      END IF
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%params) .or. ASSOCIATED(this%signals)) THEN
         status = nf_def_dim(this%result_ncid, 'string_len',                   &
     &                       data_name_length, string_len_dim_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%params)) THEN
         status = nf_def_dim(this%result_ncid, 'nparam',                       &
     &                       SIZE(this%params), nparam_dim_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_dim(this%result_ncid, 'nparamindex',                  &
     &                       2, nparamindex_dim_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%signals)) THEN
         status = nf_def_dim(this%result_ncid, 'nsignal',                      &
     &                       SIZE(this%signals), nsignal_dim_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

!  Define variables
      status = nf_def_var(this%result_ncid, 'nsteps', nf_int, 0, 0,            &
     &                    this%nsteps_id)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%recon)) THEN
         status = nf_def_var(this%result_ncid, 'g2', nf_double, 1,             &
     &                       (/ maxnsteps_dim_id /), this%g2_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%params)) THEN
         status = nf_def_var(this%result_ncid, 'param_name', nf_char,          &
     &                       2, (/ string_len_dim_id,                          &
     &                             nparam_dim_id /), param_name_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'param_index', nf_int,          &
     &                       2, (/ nparamindex_dim_id,                         &
     &                             nparam_dim_id /), param_index_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'param_value', nf_double,       &
     &                       2, (/ nparam_dim_id,                              &
     &                             maxnsteps_dim_id /),                        &
     &                       this%param_value_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'param_sigma', nf_double,       &
     &                       2, (/ nparam_dim_id,                              &
     &                             maxnsteps_dim_id /),                        &
     &                       this%param_sigma_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%signals)) THEN
         status = nf_def_var(this%result_ncid, 'signal_name', nf_char,         &
     &                       2, (/ string_len_dim_id,                          &
     &                             nsignal_dim_id /), signal_name_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'signal_type', nf_char,         &
     &                       2, (/ string_len_dim_id,                          &
     &                             nsignal_dim_id /), signal_type_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'signal_weight',                &
     &                       nf_double, 1, (/ nsignal_dim_id /),               &
     &                       signal_weight_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'signal_observed_value',        &
     &                       nf_double, 2, (/ nsignal_dim_id,                  &
     &                                        maxnsteps_dim_id /),             &
     &                       this%signal_observed_value_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'signal_model_value',           &
     &                       nf_double, 2, (/ nsignal_dim_id,                  &
     &                                        maxnsteps_dim_id /),             &
     &                       this%signal_model_value_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

         status = nf_def_var(this%result_ncid, 'signal_sigma',                 &
     &                       nf_double, 2, (/ nsignal_dim_id,                  &
     &                                        maxnsteps_dim_id /),             &
     &                       this%signal_sigma_value_id)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%model)) THEN
         CALL equilibrium_def_result(this%model%equilibrium,                   &
     &                               this%result_ncid,                         &
     &                               maxnsteps_dim_id)
      END IF

! Finished defining netcdf file. Exit out of define mode.
      status = nf_enddef(this%result_ncid)
      CALL assert_eq(status, nf_noerr, nf_strerror(status))

      IF (ASSOCIATED(this%params)) THEN
         DO i = 1, SIZE(this%params)
            status = nf_put_vara_text(this%result_ncid, param_name_id,         &
     &                                (/ 1, i /),                              &
     &                                (/ data_name_length, 1 /),               &
     &                                param_get_name(this%params(i)%p,         &
     &                                               this%model))
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

            status = nf_put_vara_int(this%result_ncid,                         &
     &                               param_index_id,                           &
     &                               (/ 1, i /), (/ 2, 1 /),                   &
     &                               this%params(i)%p%indicies)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))
         END DO
      END IF

      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            status = nf_put_vara_text(this%result_ncid, signal_name_id,        &
     &                                (/ 1, i /),                              &
     &                                (/ data_name_length, 1 /),               &
     &                                this%signals(i)%p%s_name)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

            status = nf_put_vara_text(this%result_ncid, signal_type_id,        &
     &                                (/ 1, i /),                              &
     &                                (/ data_name_length, 1 /),               &
     &                                signal_get_signal_type(                  &
     &                                   this%signals(i)%p))
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

            status = nf_put_vara_double(this%result_ncid,                      &
     &                                  signal_weight_id,                      &
     &                                  (/ 1, i /), (/ 2, 1 /),                &
     &                                  this%signals(i)%p%weight)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))
         END DO
      END IF

      CALL constext_write_step_data(this, .true.)

      IF (ASSOCIATED(this%model)) THEN
         CALL equilibrium_write_init_data(this%model%equilibrium,              &
     &                                    this%result_ncid)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!  Write step data to the defined variables.
!-------------------------------------------------------------------------------
      SUBROUTINE constext_write_step_data(this, first_step)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this
      LOGICAL, INTENT(in)                       :: first_step

!  local variables
      INTEGER                                   :: status
      INTEGER                                   :: i
      INTEGER                                   :: current_step

!  Start of executable code
!  Update number of steps taken.
      IF (first_step) THEN
         status = nf_get_var_int(this%result_ncid, this%nsteps_id, 0)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         current_step = 1
      ELSE
         status = nf_get_var_int(this%result_ncid, this%nsteps_id,             &
     &                           current_step)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         current_step = current_step + 1
         status = nf_put_var_int(this%result_ncid, this%nsteps_id,             &
     &                        current_step)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%recon)) THEN
         status = nf_put_var1_double(this%result_ncid, this%g2_id,             &
     &                               current_step,                             &
     &                               reconstruction_get_g2(this%recon))
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
      END IF

      IF (ASSOCIATED(this%params)) THEN
         DO i = 1, SIZE(this%params)
            CALL param_write_step_data(this%params(i)%p, this%model,           &
     &                                 this%result_ncid, current_step,         &
     &                                 i, this%param_value_id,                 &
     &                                 this%param_sigma_id)
         END DO
      END IF

      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            CALL signal_write_step_data(this%signals(i)%p, this%model,         &
     &                                  this%result_ncid,                      &
     &                                  current_step, i,                       &
     &                                  this%signal_observed_value_id,         &
     &                                  this%signal_model_value_id,            &
     &                                  this%signal_sigma_value_id)
         END DO
      END IF

      END SUBROUTINE

      END MODULE
