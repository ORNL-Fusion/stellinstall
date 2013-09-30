!*******************************************************************************
!>  @file v3fit_context.f
!>  @brief Contains module @ref v3fit_context.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines a @ref v3fit_context_class object to contain all the memory for
!>  running v3fit. Contains methods to write memory to disk.
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
!>  Base class representing a v3fit context. This contains all memory needed to
!>  operate v3fit.
!-------------------------------------------------------------------------------
      TYPE v3fit_context_class
         REAL (rprec)                                 :: recon_stop

!>  File input output unit for the runlog file. This file is a human readable
!>  output history of a v3fit run.
         INTEGER                                      :: runlog_iou
!>  File input output unit for the recout file. This file is a human readable
!>  output of the final result.
         INTEGER                                      :: recout_iou

!>  The parsed command line options.
         TYPE (commandline_parser_class), POINTER :: cl_parser => null()
!>  The equilibrium model.
         TYPE (model_class), POINTER                  :: model => null()
!>  Arrays of signals.
!>  @see signal
         TYPE (signal_pointer), DIMENSION(:), POINTER ::                       &
     &      signals => null()
!>  Arrays of reconstruction parameters.
!>  @see param
         TYPE (param_pointer), DIMENSION(:), POINTER :: params => null()
         TYPE (reconstruction_class), POINTER         :: recon => null()

!  The default for the following index is -1 indicating that this signal hasn't
!  been created.
!>  Index of the first instance of a @ref magnetic signal. The default -1
!>  represents no @ref magnetic signal created.
         INTEGER :: magnetic_index = -1
!>  Index of the first instance of a @ref sxrem signal. The default -1
!>  represents no @ref sxrem signal created.
         INTEGER :: sxrem_index = -1
!>  Index of the first instance of a @ref intpol signal. The default -1
!>  represents no @ref intpol signal created.
         INTEGER :: intpol_index = -1
!>  Index of the first instance of a @ref thomson signal. The default -1
!>  represents no @ref thomson signal created.
         INTEGER :: thomson_index = -1
!>  Index of the first instance of a @ref extcurz signal. The default -1
!>  represents no @ref extcurz signal created.
         INTEGER :: extcurz_index = -1
!>  Index of the first instance of a @ref limiter signal. The default -1
!>  represents no @ref limiter signal created.
         INTEGER :: limiter_index = -1
!>  Index of the first instance of a @ref prior_guassian signal. The default -1
!>  represents no @ref prior_guassian signal created.
         INTEGER :: prior_guassian_index = -1
!>  Index of the first instance of a @ref combination_class signal. The default
!>  -1 represents no @ref combination_class signal created.
         INTEGER :: combination_index = -1

!  Result file NetCDF variables.
!>  NetCDF result file to write out the results in machine readable form.
         INTEGER :: result_ncid
!>  NetCDF variable id for the number of current steps taken in the
!>  reconstruction.
         INTEGER :: nsteps_id
!>  NetCDF variable id contain the g^2 value for each reconstruction step.
         INTEGER :: g2_id

!>  NetCDF variable id contain the parameter values for each reconstruction
!>  step.
!>  @see param
         INTEGER :: param_value_id
!>  NetCDF variable id contain the parameter sigmas for each reconstruction
!>  step.
!>  @see param
         INTEGER :: param_sigma_id
!>  NetCDF variable id contain the signal effectiveness matrix for each for each
!>  reconstruction step.
!>  @see param
         INTEGER :: param_sem_id

!>  NetCDF variable id contain the parameter values for each reconstruction
!>  step.
!>  @see signal
         INTEGER :: signal_model_value_id
!>  NetCDF variable id contain the parameter values for each reconstruction
!>  step.
!>  @see signal
         INTEGER :: signal_sigma_value_id
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref v3fit_context_class object.
!>
!>  Allocates memory and initializes a @ref v3fit_context_class object with an
!>  instance of a @ref commandline_parser object. This opens files for the
!>  recout, runlog, and result files. Over allocates space for the number of
!>  signals.
!>
!>  @param[in] cl_parser An instance of a commandline parser object.
!>  @returns A pointer to a constructed @ref v3fit_context_class object.
!-------------------------------------------------------------------------------
      FUNCTION v3fit_context_construct(cl_parser)
      USE safe_open_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), POINTER  :: v3fit_context_construct
      TYPE (commandline_parser_class), POINTER, INTENT(in) :: cl_parser

!  local variables
      CHARACTER (len=data_filename_length) :: filename
      INTEGER                              :: status

!  Start of executable code
      filename = commandline_parser_get_string(cl_parser, '-file')

      ALLOCATE(v3fit_context_construct)

!  Setup the runlog file.
!  Initalize a default value of the I\O unit. V3FIT increments from there.
      v3fit_context_construct%runlog_iou = 0
      CALL safe_open(v3fit_context_construct%runlog_iou, status,               &
     &               TRIM('runlog.' // filename_base(filename)),               &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0, status, 'v3fit_context_construct: ' //                 &
     &               'Safe_open of runlog. ' //                                &
     &               TRIM(filename_base(filename)) // 'failed')

!  Setup the recout file.
!  Initalize a default value of the I\O unit. V3FIT increments from there.
      v3fit_context_construct%recout_iou = 0
      CALL safe_open(v3fit_context_construct%recout_iou, status,               &
     &               TRIM('recout.' // filename_base(filename)),               &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0, status, 'v3fit_context_construct: ' //                 &
     &               'Safe_open of recout. ' //                                &
     &               TRIM(filename_base(filename)) // 'failed')

!  Create a result file.
      status = nf_create('result.' // TRIM(filename_base(filename)) //         &
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

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref v3fit_context_class object.
!>
!>  Deallocates memory and uninitializes a @ref v3fit_context_class object. This
!>  all the v3fit constructed objects are destroyed here.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!-------------------------------------------------------------------------------
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

!  Reset the signal index
      this%magnetic_index = -1
      this%sxrem_index = -1
      this%intpol_index = -1
      this%thomson_index = -1
      this%extcurz_index = -1
      this%limiter_index = -1
      this%prior_guassian_index = -1
      this%combination_index = -1

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Resize the arrays.
!>
!>  Resizes the signal arrays to contain only the minimun necessary to hold all
!>  signals. If no signals were constructed, the signal arrays is not allocated.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
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

!  Deallocate the array is no signals are used. Reallocate the array if it is
!  smaller.
         IF (minsize .eq. 0) THEN
            DEALLOCATE(this%signals)
            this%signals => null()
            RETURN
         ELSE IF (minsize .lt. SIZE(this%signals)) THEN
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
!>  @brief Write the @ref v3fit_context_class out to disk.
!>
!>  Writes out the results of v3fit to the recout file.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!>  @note The 'this' parameter is intent(inout) because @ref signal_write calls
!>  @ref signal_get_g2 which needs to be intent(inout).
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_write(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this

!  local variables
      INTEGER                                   :: i, j
      REAL (rprec), DIMENSION(:), ALLOCATABLE   :: sem_row
      INTEGER                                   :: sem_offset
!  Fortran 95 doesn't allow allocatable strings. Need to allocate a string large
!  enough to hold every possible parameter. Fortran 95 doesn't allow allocatable
!  scalar types so make this a pointer. This large string needs to be allocated
!  to avoid a stack overflow.
      CHARACTER (len=26 + 14*v3fit_max_parameters), POINTER ::                 &
     &   sem_header

!  Start of executable code
      WRITE (*,*) ' *** Writting context to disk'
      WRITE (this%runlog_iou,*) ' *** Writting context to disk'

!  Make the sem_header null so it may be allocated and deallocated.
      sem_header => null()

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

!  Write out the signal effectiveness matrix. This matrix will be written out
!  so that each parameter is one a column and each row is a a signal. Start by
!  generating the header. The header is formatted as
!
!     #  signal_name         (param_name[12]) ...
!
!  Start my allocating a string large enough to hold the entire header. The
!  length of this string needs to be 26 for the base an additional 14 for each
!  parameter.
         WRITE (this%recout_iou, *)
         WRITE (this%recout_iou, *) ' *** Signal Effectiveness Matrix'
         ALLOCATE(sem_header)

         WRITE (sem_header, 1000)
         sem_offset = 26
         DO i = 1, SIZE(this%params)
            WRITE (sem_header, 1001) sem_header(1:sem_offset),                 &
     &                               param_get_name(this%params(i)%p,          &
     &                                              this%model)
            sem_offset = sem_offset + 14
         END DO

!  Write out the header.
         WRITE (this%recout_iou, 1002)                                         &
     &      sem_header(1:26 + 14*SIZE(this%params))

         DEALLOCATE(sem_header)
1000  FORMAT (3x,'#',2x,'s_name',14x)
1001  FORMAT (a,2x,a12)
1002  FORMAT (a)

!  Write out each row of the sem matrix.
!  Reuse the sem_header string to generate the format string. The format string
!  will be
!
!  i4,2x,a20,<num_params>(2x,es12.5)
         IF (ASSOCIATED(this%signals)) THEN
            ALLOCATE(sem_row(SIZE(this%params)))
            ALLOCATE(sem_header)
            WRITE (sem_header, 1003) SIZE(this%params)

            DO i = 1, SIZE(this%signals)
               DO j = 1, SIZE(this%params)
                  sem_row(j) = this%params(j)%p%sem(i)
               END DO

               WRITE (this%recout_iou, sem_header(1:26))                       &
     &            i, this%signals(i)%p%s_name, sem_row
            END DO

!  Write out the total for each column. Sum over all the signals for a single
!  parameter. This should be 1.0.
            WRITE (sem_header, 1004) SIZE(this%params)
            DO j = 1, SIZE(this%params)
               sem_row(j) = SUM(this%params(j)%p%sem)
            END DO

            WRITE (this%recout_iou, *)
            WRITE (this%recout_iou, sem_header) sem_row

            DEALLOCATE(sem_header)
            DEALLOCATE(sem_row)
         END IF
1003  FORMAT ('(i4,2x,a20,',i3,'(2x,es12.5))')
1004  FORMAT ('(20x,''Total:''',i3,'(2x,es12.5))')
      END IF

!  Write out the model.
      IF (ASSOCIATED(this%model)) THEN
         CALL model_write(this%model, this%recout_iou)
      END IF

!  Write out the signals
      IF (ASSOCIATED(this%signals)) THEN
         WRITE (this%recout_iou, *)
         WRITE (this%recout_iou, *) ' *** Signals'

         DO i = 1, SIZE(this%signals)
!  Check if the index has reached a new signal type
            IF ((this%magnetic_index .eq. i) .or.                              &
     &          (this%sxrem_index .eq. i) .or.                                 &
     &          (this%intpol_index .eq. i) .or.                                &
     &          (this%thomson_index .eq. i) .or.                               &
     &          (this%extcurz_index .eq. i) .or.                               &
     &          (this%limiter_index .eq. i) .or.                               &
     &          (this%prior_guassian_index .eq. i) .or.                        &
     &          (this%combination_index .eq. i)) THEN

               CALL signal_write_header(this%signals(i)%p,                     &
     &                                  this%recout_iou)
            END IF

            CALL signal_write(this%signals(i)%p, this%recout_iou, i,           &
     &                        this%model)
         END DO

!  All the signals should be written out before auxiliary information.
         DO i = 1, SIZE(this%signals)
            CALL signals_write_auxiliary(this%signals(i)%p,                    &
     &                                   this%recout_iou, i,                   &
     &                                   this%model)
         END DO
      END IF

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!>  @page result_file_main Result File
!>
!>  @tableofcontents
!>  @section result_file_main_intro_sec Introduction
!>  This page doctuments the contents of a result NetCDF file.
!>
!>  @section result_file_main_dim_sec Dimensions
!>  @header{Dimension, Description, Code Reference}
!>  @begin_table
!>     @item{maxnsteps,   The max number of reconstruction steps., v3fit_input::nrstep}
!>     @item{nparam,      Number of reconstruction parameters.,    v3fit_input::n_rp}
!>     @item{nparamindex, Number of parameter indices.,            data_parameters::data_max_indicies}
!>     @item{nsignal,     Number of constructed signals.,          v3fit_context::v3fit_context_class::signals}
!>     @item{string_len,  Max lenth of strings.,                   data_parameters::data_name_length}
!>  @end_table
!>  @subsection result_file_main_eq_dim_sec Equilibrium Dimensions
!>  Equilibrium dimensions change depending on the equilibrium used.
!>  * VMEC @ref result_file_vmec_dim_sec
!>
!>  @section result_file_main_var_sec Variables
!>  @header{Variable(Dimensions), Description, Code Reference}
!>  @begin_table
!>     @item{nsteps,         Actual number of steps taken.,           reconstruction::reconstruction_class::current_step}
!>     @item{g2 (maxnsteps), g^2 value for each reconstruction step., reconstruction::reconstruction_get_g2}
!>  @end_table
!>  @table_section{result_file_main_param_var_sec, Parameter Variables}
!>     @item{param_name (string_len\, nparam),                 Name of the reconstruction parameter.,              param::param_get_name}
!>     @item{param_index (nparamindex\, nparam),               Indices of the reconstruction parameter.,           param::param_class::indicies}
!>     @item{param_value (nparam\, maxnsteps),                 Parameter values for each reconstruction step.,     param::param_get_value}
!>     @item{param_sigma (nparam\, maxnsteps),                 Parameter sigmas for each reconstruction step.,     param::param_class::sigma}
!>     @item{signal_eff_matrix (nsignal\, nparam\, maxnsteps), Signal effectiveness for each reconstruction step., param::param_class::sem}
!>  @end_table
!>  @table_section{result_file_main_sig_var_sec, Signal Variables}
!>     @item{signal_name (string_len\, nsignal),       Short name of the signal.,                                signal::signal_class::s_name}
!>     @item{signal_type (string_len\, nsignal),       Description of the signal.,                               signal::signal_get_signal_type}
!>     @item{signal_weight (nsignal),                  Weight value of the signal.,                              signal::signal_class::weight}
!>     @item{signal_observed_value (nsignal),          Observed value of the signal.,                            signal::signal_class::observed}
!>     @item{signal_model_value (nsignal\, maxnsteps), Modeled value of the signal at each reconstruction step., signal::signal_get_modeled_signal}
!>     @item{signal_sigma_value (nsignal\, maxnsteps), Sigma value of the signal at each reconstruction step.,   signal::signal_get_g2}
!>  @end_table
!>  @subsection result_file_main_eq_var_sec Equilibrium Variables
!>  Equilibrium variables change depending on the equilibrium used.
!>  * VMEC @ref result_file_vmec_var_sec
!-------------------------------------------------------------------------------
!>  @brief Initialize the dimensions and variables of the result file.
!>
!>  Defines dimensions, variables and writes out the initial step of v3fit to
!>  result netcdf file. Multi dimensional arrays need to be transposed so arrays
!>  appear in the correct order in non fortran languages.
!>
!>  @param[inout] this A @ref v3fit_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_init_data(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3fit_context_class), INTENT(inout) :: this

!  local variables
      INTEGER :: status
      INTEGER :: i
      INTEGER :: maxnsteps_dim_id
      INTEGER :: nparam_dim_id
      INTEGER :: nparamindex_dim_id
      INTEGER :: nsignal_dim_id
      INTEGER :: string_len_dim_id
      INTEGER :: param_name_id
      INTEGER :: param_index_id
      INTEGER :: signal_name_id
      INTEGER :: signal_type_id
      INTEGER :: signal_weight_id
      INTEGER :: signal_observed_value_id
!  Start of executable code
!  Define dimensions
      IF (ASSOCIATED(this%recon)) THEN
         status = nf_def_dim(this%result_ncid, 'maxnsteps',                    &
     &                       SIZE(this%recon%e, 1), maxnsteps_dim_id)
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

!  Define parameter variables.
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

         status = nf_def_var(this%result_ncid, 'signal_eff_matrix',            &
     &                       nf_double, 3, (/ nsignal_dim_id,                  &
     &                                        nparam_dim_id,                   &
     &                                        maxnsteps_dim_id /),             &
     &                       this%param_sem_id)
      END IF

!  Define signal variables.
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
     &                       nf_double, 1, (/ nsignal_dim_id /),               &
     &                       signal_observed_value_id)
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

            status = nf_put_var1_double(this%result_ncid,                      &
     &                                  signal_weight_id, i,                   &
     &                                  this%signals(i)%p%weight)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))

            status = nf_put_var1_double(this%result_ncid,                      &
     &                                  signal_observed_value_id, i,           &
     &                                  this%signals(i)%p%observed)
            CALL assert_eq(status, nf_noerr, nf_strerror(status))
         END DO
      END IF

      CALL v3fit_context_write_step_data(this, .true.)

      IF (ASSOCIATED(this%model)) THEN
         CALL equilibrium_write_init_data(this%model%equilibrium,              &
     &                                    this%result_ncid)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write step data to the defined variables.
!>
!>  Writes out the data for a reconstruction step to the result file. If this
!>  is being called for the first time. Do not read the number of steps.
!>
!>  @param[inout] this       A @ref v3fit_context_class instance.
!>  @param[in]    first_step Flags to determine if this is the first instance
!>                           this has been called.
!>  @note The only caller of this subrotine that should set first_step to true
!>  should be @ref v3fit_context_init_data. All others should callers should set
!>  first_step to false.
!-------------------------------------------------------------------------------
      SUBROUTINE v3fit_context_write_step_data(this, first_step)

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
         status = nf_put_var_int(this%result_ncid, this%nsteps_id, 0)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         current_step = 1
      ELSE
         status = nf_get_var_int(this%result_ncid, this%nsteps_id,             &
     &                           current_step)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))
         current_step = current_step + 1
         status = nf_put_var_int(this%result_ncid, this%nsteps_id,             &
     &                           current_step)
         CALL assert_eq(status, nf_noerr, nf_strerror(status))

!  Use current_step as an array index. The Netcdf arrays start at 1 so this
!  needs to be incremented to point to the correct index.
         current_step = current_step + 1
      END IF

      CALL equilibrium_write_step_data(this%model%equilibrium,                 &
     &                                 this%result_ncid, current_step)

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
     &                                 this%param_sigma_id,                    &
     &                                 this%param_sem_id)
         END DO
      END IF

      IF (ASSOCIATED(this%signals)) THEN
         DO i = 1, SIZE(this%signals)
            CALL signal_write_step_data(this%signals(i)%p, this%model,         &
     &                                  this%result_ncid,                      &
     &                                  current_step, i,                       &
     &                                  this%signal_model_value_id,            &
     &                                  this%signal_sigma_value_id)
         END DO
      END IF

      END SUBROUTINE

      END MODULE
