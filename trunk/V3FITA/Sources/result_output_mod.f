!*******************************************************************************
!  File result_output_mod.f
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  This modules containes code to output the results of a v3fit reconstruction 
!  into a computer readable netcdf file.
!
!*******************************************************************************

      MODULE result_output
      USE netcdf
      USE gsq_mod, ONLY: gsq_g2
      USE signal_T
      USE v3f_global, ONLY: rparam, s_desc

!*******************************************************************************
!  Object representing the results file. This is mean to be an opaque. All 
!  fields private and can only be accessed from inside this module.
!
!  NetCDF variables:
!    ncid    Reference to the netcdf file.
!    status  Status value of the netcdf file. If this evaluates to anything
!            that nf90_noerror, then an assert is triggered.
!
!  Dimensions: Reference IDs for the created dimensions.
!    maxnsteps_dim_id  The max number of reconstruction steps.
!    nparam_dim_id     Number of reconstruction parameters.
!    nsignal_dim_id    Number of signals.
!
!  Variables: Reference IDs for the created variables.
!    nsteps_id (scalar) Actual number of steps taken.
!    
!    g2_id (maxnsteps_dim_id) g^2 value for each reconstruction step.
!
!    param_name_id (nparam_dim_id) Name of the reconstruction parameters.
!    param_index_id (nparam_dim_id) First indices of the reconstruction 
!                   parameters.
!    param_value_id (maxnsteps_dim_id, nparam_dim_id) Parameter values for each
!                   reconstruction step.
!    param_sigma_id (maxnsteps_dim_id, nparam_dim_id) Parameter sigmas for each
!                   reconstruction step.
!
!    signal_name_id (nsignal_dim_id) Long name of the signals.
!    signal_type_id (nsignal_dim_id) Type of the signals.
!    signal_obser_value_id (maxnsteps_dim_id, nsignal_dim_id) Observed signals 
!                          for each reconstruction step.
!    signal_model_value_id (maxnsteps_dim_id, nsignal_dim_id) Modeled signals 
!                          for each reconstruction step.
!    signal_sigma_id (maxnsteps_dim_id, nsignal_dim_id) Signal sigmas 
!                          for each reconstruction step.
!
!*******************************************************************************
      TYPE result_file
         PRIVATE
         INTEGER :: ncid, status
! Dimension IDs
         INTEGER :: maxnsteps_dim_id, nparam_dim_id, nsignal_dim_id
         INTEGER :: string_len_id
! Variable IDs
         INTEGER :: nsteps_id
         INTEGER :: g2_id
         INTEGER :: param_name_id, param_index_id, param_value_id,             &
     &              param_sigma_id
         INTEGER :: signal_name_id, signal_type_id,                            &
     &              signal_obser_value_id,                                     &
     &              signal_model_value_id, signal_sigma_id
      END TYPE result_file

      PUBLIC  :: result_output_file_create, result_output_file_write_step,     &
     &           result_output_file_close
      PRIVATE :: result_output_check_error
      
      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Create and setup a new netcdf file.
!
!  A discription of the resulting NetCDF file is provided.
!    
!  Dimensions:
!    maxnsteps  Maximum number of steps.
!    nparam     Number of parameters.
!    nsignal    Number of signals.
!
!  Variables:
!    nsteps (INT) Actual number of reconstruction steps. This starts from zero.
!
!    g2 (DOUBLE) g^2 of each reconstruction step.
!
!    param_name (STRING) Name of the reconstruction parameter.
!    param_index (INT) Index of the reconstruction parameter.
!    param_value (DOUBLE) Value of the reconstruction parameter at each 
!                         reconstruction step.
!    param_sigma (DOUBLE) Sigma of the reconstruction parameter at each 
!                         reconstruction step.
!
!    signal_name (STRING) Long name of the signal.
!    signal_type (STRING) Description of the signal.
!    signal_observed_value (DOUBLE) Observed value of the signal at each 
!                                   reconstruction step.
!    signal_model_value (DOUBLE) Modeled value of the signal at each 
!                                reconstruction step.
!    signal_sigma (DOUBLE) Sigma value of the reconstruction step.
!                          sig = Sqrt( Model_sig^2 + Observe_sig^2 )
!
!*******************************************************************************
      FUNCTION result_output_file_create(path, sd_model, sd_observed)
      USE v3f_global, ONLY: nrstep, n_rparam, n_sdata_o
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      CHARACTER (len = *), INTENT(in) :: path
      TYPE (signal_data), DIMENSION(:), INTENT (in) ::                         &
     &   sd_model, sd_observed
      TYPE (result_file) :: result_output_file_create
      INTEGER :: i, beginIndex, endIndex

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Fine the start and end index of the of the "." to get the base name of the
! name list input file.

      beginIndex = INDEX(path, "input") + 5
      endIndex = INDEX(path, "v3fit", .FALSE.)

      result_output_file_create%status =                                       &
     &   nf90_create(path  = "result" //                                       &
     &                       TRIM(path(beginIndex:endIndex)) //                &
     &                       "nc",                                             &
     &               cmode = nf90_noclobber,                                   &
     &               ncid  = result_output_file_create%ncid)        
      CALL result_output_check_error(result_output_file_create%status)

! Define dimensions for data.
      result_output_file_create%status =                                       &
     &   nf90_def_dim(result_output_file_create%ncid,                          &
     &                "maxnsteps",                                             &
     &                nrstep,                                                  &
     &                result_output_file_create%maxnsteps_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_dim(result_output_file_create%ncid,                          &
     &                "nparam",                                                &
     &                n_rparam,                                                &
     &                result_output_file_create%nparam_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_dim(result_output_file_create%ncid,                          &
     &                "nsignal",                                               &
     &                n_sdata_o,                                               &
     &                result_output_file_create%nsignal_dim_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_dim(result_output_file_create%ncid,                          &
     &                "string_len",                                            &
     &                80,                                                      &
     &                result_output_file_create%string_len_id)
      CALL result_output_check_error(result_output_file_create%status)

! Define variables for data.
      result_output_file_create%status =                                       &
     &   nf90_def_var(ncid  = result_output_file_create%ncid,                  &
     &                name  = "nsteps",                                        &
     &                xtype = nf90_int,                                        &
     &                varid = result_output_file_create%nsteps_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "g2",                                                    &
     &                nf90_double,                                             &
     &                result_output_file_create%maxnsteps_dim_id,              &
     &                result_output_file_create%g2_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "param_name",                                            &
     &                nf90_char,                                               &
     &                (/ result_output_file_create%nparam_dim_id,              &
     &                   result_output_file_create%string_len_id /),           & 
     &                result_output_file_create%param_name_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "param_index",                                           &
     &                nf90_int,                                                &
     &                result_output_file_create%nparam_dim_id,                 &
     &                result_output_file_create%param_index_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "param_value",                                           &
     &                nf90_double,                                             &
     &                (/ result_output_file_create%maxnsteps_dim_id,           &
     &                   result_output_file_create%nparam_dim_id /),           &
     &                result_output_file_create%param_value_id)                &
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "param_sigma",                                           &
     &                nf90_double,                                             &
     &                (/ result_output_file_create%maxnsteps_dim_id,           &
     &                   result_output_file_create%nparam_dim_id /),           &
     &                result_output_file_create%param_sigma_id)                &
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "signal_name",                                           &
     &                nf90_char,                                               &
     &                (/ result_output_file_create%nparam_dim_id,              &
     &                   result_output_file_create%string_len_id /),           & 
     &                result_output_file_create%signal_name_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "signal_type",                                            &
     &                nf90_char,                                               &
     &                (/ result_output_file_create%nparam_dim_id,              &
     &                   result_output_file_create%string_len_id /),           & 
     &                result_output_file_create%signal_type_id)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "signal_observed_value",                                 &
     &                nf90_double,                                             &
     &                (/ result_output_file_create%maxnsteps_dim_id,           &
     &                   result_output_file_create%nsignal_dim_id /),          &
     &                result_output_file_create%signal_obser_value_id)         &
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "signal_model_value",                                    &
     &                nf90_double,                                             &
     &                (/ result_output_file_create%maxnsteps_dim_id,           &
     &                   result_output_file_create%nparam_dim_id /),           &
     &                result_output_file_create%signal_model_value_id)         &
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_def_var(result_output_file_create%ncid,                          &
     &                "signal_sigma",                                          &
     &                nf90_double,                                             &
     &                (/ result_output_file_create%maxnsteps_dim_id,           &
     &                   result_output_file_create%nparam_dim_id /),           &
     &                result_output_file_create%signal_sigma_id)               &
      CALL result_output_check_error(result_output_file_create%status)

! Finished defining netcdf file. Exit out of define mode.
      result_output_file_create%status =                                       &
     &   nf90_enddef(result_output_file_create%ncid)
      CALL result_output_check_error(result_output_file_create%status)

! Fill in initial variable values.
      result_output_file_create%status =                                       &
     &   nf90_put_var(result_output_file_create%ncid,                          &
     &                result_output_file_create%nsteps_id,                     &
     &                0)
      CALL result_output_check_error(result_output_file_create%status)

      result_output_file_create%status =                                       &
     &   nf90_put_var(result_output_file_create%ncid,                          &
     &                result_output_file_create%g2_id,                         & 
     &                gsq_g2)
      CALL result_output_check_error(result_output_file_create%status)

! Fill in reconstuction parameters.
      DO i = 1, n_rparam
         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%param_name_id,              &
     &                   rparam(i)%p_type,                                     &
     &                   start = (/ i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%param_index_id,             &
     &                   rparam(i)%index,                                      &
     &                   start = (/ i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%param_value_id,             &
     &                   rparam(i)%value,                                      &
     &                   start = (/ 1, i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%param_index_id,             &
     &                   rparam(i)%vrnc,                                       &
     &                   start = (/ 1, i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)
      END DO

! Fill in signals
      DO i = 1, n_sdata_o
         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%signal_name_id,             &
     &                   s_desc(i)%l_name,                                &
     &                   start = (/ i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%signal_type_id,             &
     &                   s_desc(i)%s_type,                                &
     &                   start = (/ i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(                                                      &
     &         result_output_file_create%ncid,                                 &
     &         result_output_file_create%signal_obser_value_id,                &
     &         sd_observed(i)%data(1),                                         &
     &         start = (/ 1, i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(                                                      &
     &         result_output_file_create%ncid,                                 &
     &         result_output_file_create%signal_model_value_id,                &
     &         sd_model(i)%data(1),                                            &
     &         start = (/ 1, i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)

         result_output_file_create%status =                                    &
     &      nf90_put_var(result_output_file_create%ncid,                       &
     &                   result_output_file_create%signal_sigma_id,            &
     &                   SQRT(sd_model(i)%sigma(1)**2 +                        &
     &                        sd_observed(i)%sigma(1)**2),                     &
     &                   start = (/ 1, i /))
         CALL                                                                  &
     &      result_output_check_error(result_output_file_create%status)
      END DO

      END FUNCTION result_output_file_create

!*******************************************************************************
!  Write reconstruction step to the netcdf file.
!*******************************************************************************
      SUBROUTINE result_output_file_write_step(result_file_ref,                &
     &                                         sd_model,                       &
     &                                         sd_observed)
      USE v3f_global, ONLY: n_rparam, n_sdata_o
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE (result_file) :: result_file_ref
      TYPE (signal_data), DIMENSION(:), INTENT (in) ::                         &
     &   sd_model, sd_observed
      INTEGER :: nsteps, i

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Read and update current step.
      result_file_ref%status =                                                 &
     &   nf90_get_var(result_file_ref%ncid,                                    &
     &                result_file_ref%nsteps_id,                               &
     &                nsteps)
      CALL result_output_check_error(result_file_ref%status)
      nsteps = nsteps + 1
      result_file_ref%status =                                                 &
     &   nf90_put_var(result_file_ref%ncid,                                    &
     &                result_file_ref%nsteps_id,                               &
     &                nsteps)
      CALL result_output_check_error(result_file_ref%status)

! Steps start at zero however array offsets in the netcdf routines start at one. 
! Increment nsteps to reflect the correct array position.
      nsteps = nsteps + 1

      result_file_ref%status =                                                 &
     &   nf90_put_var(result_file_ref%ncid,                                    &
     &                result_file_ref%g2_id,                                   & 
     &                gsq_g2,                                                  &
     &                start = (/ nsteps /))
      CALL result_output_check_error(result_file_ref%status)

! Fill in reconstuction parameters.
      DO i = 1, n_rparam
         result_file_ref%status =                                              &
     &      nf90_put_var(result_file_ref%ncid,                                 &
     &                   result_file_ref%param_value_id,                       &
     &                   rparam(i)%value,                                      &
     &                   start = (/ nsteps, i /))
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf90_put_var(result_file_ref%ncid,                                 &
     &                   result_file_ref%param_index_id,                       &
     &                   rparam(i)%vrnc,                                       &
     &                   start = (/ nsteps, i /))
         CALL result_output_check_error(result_file_ref%status)
      END DO

! Fill in signals
      DO i = 1, n_sdata_o
         result_file_ref%status =                                              &
     &      nf90_put_var(                                                      &
     &         result_file_ref%ncid,                                           &
     &         result_file_ref%signal_obser_value_id,                          &
     &         sd_observed(i)%data(1),                                         &
     &         start = (/ nsteps, i /))
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf90_put_var(                                                      &
     &         result_file_ref%ncid,                                           &
     &         result_file_ref%signal_model_value_id,                          &
     &         sd_model(i)%data(1),                                            &
     &         start = (/ nsteps, i /))
         CALL result_output_check_error(result_file_ref%status)

         result_file_ref%status =                                              &
     &      nf90_put_var(result_file_ref%ncid,                                 &
     &                   result_file_ref%signal_sigma_id,                      &
     &                   SQRT(sd_model(i)%sigma(1)**2 +                        &
     &                        sd_observed(i)%sigma(1)**2),                     &
     &                   start = (/ nsteps, i /))
         CALL result_output_check_error(result_file_ref%status)
      END DO

      END SUBROUTINE

!*******************************************************************************
!  Close the netcdf file.
!*******************************************************************************
      SUBROUTINE result_output_file_close(result_file_ref)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      TYPE (result_file) :: result_file_ref

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      result_file_ref%status = nf90_close(result_file_ref%ncid)
      CALL result_output_check_error(result_file_ref%status)

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
! Check and handle errors.
!*******************************************************************************
      SUBROUTINE result_output_check_error(status)
      USE v3_utilities                                  ! subroutine assert_eq
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: status

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      CALL assert_eq(status, nf90_noerr, nf90_strerror(status))

      END SUBROUTINE

      END MODULE result_output
